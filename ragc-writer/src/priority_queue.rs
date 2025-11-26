// Bounded priority queue for streaming compression pipeline
// Matches C++ AGC's CBoundedPQueue (queue.h:153-346)

use std::collections::BinaryHeap;
use std::sync::{Arc, Condvar, Mutex};

/// Result type for pop operations
///
/// Matches C++ AGC's result_t enum (queue.h:171):
/// ```cpp
/// enum class result_t { empty, completed, normal };
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PopResult {
    /// Queue is empty but producers are still active (wait and retry)
    Empty,
    /// Queue is empty AND no producers remain (exit worker loop)
    Completed,
    /// Successfully popped an item
    Normal,
}

/// Priority queue entry
///
/// Matches C++ AGC's multimap key: `pair<size_t, size_t>` (queue.h:155)
/// - First: priority (higher = processed first)
/// - Second: cost (for capacity limiting)
#[derive(Debug, Clone, Eq, PartialEq)]
struct QueueEntry<T> {
    priority: usize,
    cost: usize,
    data: T,
}

impl<T> PartialOrd for QueueEntry<T>
where
    T: Eq,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<T> Ord for QueueEntry<T>
where
    T: Eq,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Compare by (priority, cost) - matches C++ multimap ordering
        // Note: BinaryHeap is max-heap, which matches PopLarge() behavior
        (self.priority, self.cost).cmp(&(other.priority, other.cost))
    }
}

/// Internal queue state
struct QueueState<T> {
    /// Priority queue (max-heap by priority, then cost)
    queue: BinaryHeap<QueueEntry<T>>,
    /// Number of active producers
    n_producers: usize,
    /// Current total cost in queue
    current_cost: usize,
    /// Maximum allowed cost
    max_cost: usize,
}

/// Bounded priority queue with capacity limiting
///
/// This matches C++ AGC's CBoundedPQueue (queue.h:153-346):
/// ```cpp
/// template<typename T> class CBoundedPQueue {
///     typedef multimap<pair<size_t, size_t>, T> queue_t;
///     // ...
///     void Emplace(T&& data, size_t priority, size_t cost);
///     result_t PopLarge(T& data);
///     void MarkCompleted();
/// }
/// ```
///
/// Key features:
/// - **Priority ordering**: Higher priority items processed first
/// - **Capacity limiting**: Sum of costs must stay below max_cost
/// - **Thread-safe**: Multiple producers and consumers
/// - **Completion signaling**: MarkCompleted() when producer done
///
/// Example:
/// ```no_run
/// use ragc_writer::priority_queue::{BoundedPriorityQueue, PopResult};
///
/// let queue = BoundedPriorityQueue::new(2, 1000); // 2 producers, 1000 max cost
///
/// // Producer thread
/// queue.emplace("task1".to_string(), 100, 50);  // priority 100, cost 50
/// queue.mark_completed();  // This producer is done
///
/// // Consumer thread
/// loop {
///     match queue.pop_large() {
///         (PopResult::Normal, Some(data)) => {
///             // Process data
///         }
///         (PopResult::Empty, None) => {
///             continue;  // Wait for more
///         }
///         (PopResult::Completed, None) => {
///             break;  // All done
///         }
///         _ => {
///             // Handle other cases
///         }
///     }
/// }
/// ```
pub struct BoundedPriorityQueue<T> {
    state: Arc<(Mutex<QueueState<T>>, Condvar, Condvar)>,
}

impl<T> BoundedPriorityQueue<T>
where
    T: Clone + Eq,
{
    /// Create a new bounded priority queue
    ///
    /// Matches C++ AGC's constructor (queue.h:174-180):
    /// ```cpp
    /// CBoundedPQueue(const int _n_producers, const size_t _max_cost);
    /// ```
    ///
    /// # Arguments
    /// * `n_producers` - Number of producer threads
    /// * `max_cost` - Maximum total cost allowed in queue
    pub fn new(n_producers: usize, max_cost: usize) -> Self {
        BoundedPriorityQueue {
            state: Arc::new((
                Mutex::new(QueueState {
                    queue: BinaryHeap::new(),
                    n_producers,
                    current_cost: 0,
                    max_cost,
                }),
                Condvar::new(), // cv_queue_empty
                Condvar::new(), // cv_queue_full
            )),
        }
    }

    /// Add an item to the queue (blocks if at capacity)
    ///
    /// Matches C++ AGC's Emplace (queue.h:238-251):
    /// ```cpp
    /// void Emplace(T&& data, const size_t priority, const size_t cost);
    /// ```
    ///
    /// # Arguments
    /// * `data` - Item to enqueue
    /// * `priority` - Priority (higher = processed first)
    /// * `cost` - Memory cost (for capacity limiting)
    pub fn emplace(&self, data: T, priority: usize, cost: usize) {
        let (mutex, cv_empty, cv_full) = &*self.state;
        let mut state = mutex.lock().unwrap();

        // Wait until there's space (current_cost < max_cost)
        while state.current_cost >= state.max_cost {
            state = cv_full.wait(state).unwrap();
        }

        let was_empty = state.queue.is_empty();
        state.queue.push(QueueEntry {
            priority,
            cost,
            data,
        });
        state.current_cost += cost;

        if was_empty {
            cv_empty.notify_all();
        }
    }

    /// Add multiple copies of an item with zero cost (for sync barriers)
    ///
    /// Matches C++ AGC's EmplaceManyNoCost (queue.h:270-280):
    /// ```cpp
    /// void EmplaceManyNoCost(T&& data, size_t priority, size_t n_items);
    /// ```
    ///
    /// This is used to send synchronization tokens to all workers.
    ///
    /// # Arguments
    /// * `data` - Item to enqueue (will be cloned n_items times)
    /// * `priority` - Priority
    /// * `n_items` - Number of copies to enqueue
    pub fn emplace_many_no_cost(&self, data: T, priority: usize, n_items: usize) {
        let (mutex, cv_empty, _) = &*self.state;
        let mut state = mutex.lock().unwrap();

        for _ in 0..n_items {
            state.queue.push(QueueEntry {
                priority,
                cost: 0,
                data: data.clone(),
            });
        }

        cv_empty.notify_all();
    }

    /// Pop the highest priority item (blocks if empty)
    ///
    /// Matches C++ AGC's PopLarge (queue.h:284-313):
    /// ```cpp
    /// result_t PopLarge(T& data);
    /// ```
    ///
    /// Returns:
    /// - (Normal, Some(data)): Successfully popped highest priority item
    /// - (Empty, None): Queue empty but producers still active
    /// - (Completed, None): Queue empty and no producers (exit)
    pub fn pop_large(&self) -> (PopResult, Option<T>) {
        let (mutex, cv_empty, cv_full) = &*self.state;
        let mut state = mutex.lock().unwrap();

        // Wait until queue has items or all producers are done
        while state.queue.is_empty() && state.n_producers > 0 {
            state = cv_empty.wait(state).unwrap();
        }

        if state.queue.is_empty() {
            // No items and no producers left
            return if state.n_producers > 0 {
                (PopResult::Empty, None)
            } else {
                (PopResult::Completed, None)
            };
        }

        // Pop highest priority item (BinaryHeap is max-heap)
        let entry = state.queue.pop().unwrap();
        state.current_cost -= entry.cost;

        if state.queue.is_empty() {
            cv_empty.notify_all();
        }

        cv_full.notify_all();

        (PopResult::Normal, Some(entry.data))
    }

    /// Signal that a producer is done
    ///
    /// Matches C++ AGC's MarkCompleted (queue.h:212-219):
    /// ```cpp
    /// void MarkCompleted();
    /// ```
    ///
    /// When all producers call this, consumers will receive Completed.
    pub fn mark_completed(&self) {
        let (mutex, cv_empty, _) = &*self.state;
        let mut state = mutex.lock().unwrap();

        state.n_producers -= 1;

        if state.n_producers == 0 {
            cv_empty.notify_all();
        }
    }

    /// Check if queue is empty
    ///
    /// Matches C++ AGC's IsEmpty (queue.h:197-201)
    pub fn is_empty(&self) -> bool {
        let (mutex, _, _) = &*self.state;
        let state = mutex.lock().unwrap();
        state.queue.is_empty()
    }

    /// Check if queue is completed (empty and no producers)
    ///
    /// Matches C++ AGC's IsCompleted (queue.h:204-209)
    pub fn is_completed(&self) -> bool {
        let (mutex, _, _) = &*self.state;
        let state = mutex.lock().unwrap();
        state.queue.is_empty() && state.n_producers == 0
    }

    /// Get current queue size (items, total_cost)
    ///
    /// Matches C++ AGC's GetSize (queue.h:340-345)
    pub fn get_size(&self) -> (usize, usize) {
        let (mutex, _, _) = &*self.state;
        let state = mutex.lock().unwrap();
        (state.queue.len(), state.current_cost)
    }
}

// Make queue cloneable (shares Arc internally)
impl<T> Clone for BoundedPriorityQueue<T> {
    fn clone(&self) -> Self {
        BoundedPriorityQueue {
            state: Arc::clone(&self.state),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;
    use std::time::Duration;

    #[test]
    fn test_basic_operations() {
        let queue: BoundedPriorityQueue<String> = BoundedPriorityQueue::new(1, 1000);

        queue.emplace("task1".to_string(), 100, 50);
        queue.emplace("task2".to_string(), 200, 50);
        queue.emplace("task3".to_string(), 150, 50);

        // Should pop in priority order: 200, 150, 100
        let (result, data) = queue.pop_large();
        assert_eq!(result, PopResult::Normal);
        assert_eq!(data, Some("task2".to_string()));

        let (result, data) = queue.pop_large();
        assert_eq!(result, PopResult::Normal);
        assert_eq!(data, Some("task3".to_string()));

        let (result, data) = queue.pop_large();
        assert_eq!(result, PopResult::Normal);
        assert_eq!(data, Some("task1".to_string()));
    }

    #[test]
    fn test_completion_signaling() {
        let queue: BoundedPriorityQueue<String> = BoundedPriorityQueue::new(1, 1000);

        queue.mark_completed();

        // Queue is empty and producer is done
        let (result, data) = queue.pop_large();
        assert_eq!(result, PopResult::Completed);
        assert_eq!(data, None);
    }

    #[test]
    fn test_emplace_many_no_cost() {
        let queue: BoundedPriorityQueue<String> = BoundedPriorityQueue::new(1, 1000);

        queue.emplace_many_no_cost("sync".to_string(), 500, 3);

        for _ in 0..3 {
            let (result, data) = queue.pop_large();
            assert_eq!(result, PopResult::Normal);
            assert_eq!(data, Some("sync".to_string()));
        }
    }

    #[test]
    fn test_multi_threaded() {
        let queue: BoundedPriorityQueue<String> = BoundedPriorityQueue::new(2, 1000);
        let q1 = queue.clone();
        let q2 = queue.clone();

        // Producer 1
        let p1 = thread::spawn(move || {
            for i in 0..5 {
                q1.emplace(format!("p1-{}", i), 100 + i, 10);
                thread::sleep(Duration::from_millis(1));
            }
            q1.mark_completed();
        });

        // Producer 2
        let p2 = thread::spawn(move || {
            for i in 0..5 {
                q2.emplace(format!("p2-{}", i), 200 + i, 10);
                thread::sleep(Duration::from_millis(1));
            }
            q2.mark_completed();
        });

        // Consumer
        let mut count = 0;
        loop {
            match queue.pop_large() {
                (PopResult::Normal, Some(_)) => {
                    count += 1;
                }
                (PopResult::Empty, None) => {
                    thread::sleep(Duration::from_millis(1));
                    continue;
                }
                (PopResult::Completed, None) => {
                    break;
                }
                _ => panic!("Unexpected queue state"),
            }
        }

        p1.join().unwrap();
        p2.join().unwrap();

        assert_eq!(count, 10);
    }

    #[test]
    fn test_capacity_limiting() {
        let queue: BoundedPriorityQueue<String> = BoundedPriorityQueue::new(1, 100);
        let q = queue.clone();

        // This should fill the queue (50 + 50 = 100)
        queue.emplace("task1".to_string(), 100, 50);
        queue.emplace("task2".to_string(), 100, 50);

        // Try to add one more in a separate thread (should block)
        let producer = thread::spawn(move || {
            q.emplace("task3".to_string(), 100, 50);
        });

        // Give producer time to try to enqueue
        thread::sleep(Duration::from_millis(10));

        // Should still have only 2 items
        assert_eq!(queue.get_size(), (2, 100));

        // Pop one to make space
        queue.pop_large();

        // Now producer should be able to proceed
        producer.join().unwrap();

        // Should now have 2 items again
        assert_eq!(queue.get_size(), (2, 100));
    }
}
