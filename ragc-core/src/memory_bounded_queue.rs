// Memory-bounded queue with backpressure
// Matches C++ AGC's CBoundedPQueue behavior

use std::collections::VecDeque;
use std::sync::{Arc, Condvar, Mutex};

/// A queue bounded by total bytes (not item count)
///
/// Key properties:
/// - `push()` blocks when adding would exceed capacity
/// - `pull()` blocks when queue is empty (returns None when closed)
/// - Provides automatic backpressure for constant memory usage
///
/// This matches C++ AGC's CBoundedPQueue architecture.
pub struct MemoryBoundedQueue<T> {
    inner: Arc<Mutex<QueueInner<T>>>,
    capacity_bytes: usize,
    not_full: Arc<Condvar>,
    not_empty: Arc<Condvar>,
}

struct QueueInner<T> {
    items: VecDeque<(T, usize)>, // (item, size_bytes)
    current_size: usize,         // Total bytes currently in queue
    closed: bool,                // No more pushes allowed
}

impl<T> MemoryBoundedQueue<T> {
    /// Create a new memory-bounded queue
    ///
    /// # Arguments
    /// * `capacity_bytes` - Maximum total bytes allowed in queue
    ///
    /// # Example
    /// ```
    /// use ragc_core::MemoryBoundedQueue;
    /// let queue: MemoryBoundedQueue<Vec<u8>> = MemoryBoundedQueue::new(2 * 1024 * 1024 * 1024); // 2 GB
    /// ```
    pub fn new(capacity_bytes: usize) -> Self {
        Self {
            inner: Arc::new(Mutex::new(QueueInner {
                items: VecDeque::new(),
                current_size: 0,
                closed: false,
            })),
            capacity_bytes,
            not_full: Arc::new(Condvar::new()),
            not_empty: Arc::new(Condvar::new()),
        }
    }

    /// Push an item to the queue with its size
    ///
    /// **BLOCKS** if adding this item would exceed capacity.
    /// Returns error if queue is closed.
    ///
    /// # Arguments
    /// * `item` - The item to push
    /// * `size_bytes` - Size of the item in bytes
    ///
    /// # Example
    /// ```no_run
    /// # use ragc_core::MemoryBoundedQueue;
    /// # let queue: MemoryBoundedQueue<Vec<u8>> = MemoryBoundedQueue::new(1024);
    /// let contig_data = vec![b'A'; 1000];
    /// queue.push(contig_data.clone(), contig_data.len()).unwrap(); // Blocks if queue is full!
    /// ```
    pub fn push(&self, item: T, size_bytes: usize) -> Result<(), PushError> {
        let mut inner = self.inner.lock().unwrap();

        // Wait while queue would be too full
        while inner.current_size + size_bytes > self.capacity_bytes && !inner.closed {
            inner = self.not_full.wait(inner).unwrap();
        }

        // Check if closed while we were waiting
        if inner.closed {
            return Err(PushError::Closed);
        }

        // Add item
        inner.items.push_back((item, size_bytes));
        inner.current_size += size_bytes;

        // Signal that queue is not empty
        self.not_empty.notify_one();

        Ok(())
    }

    /// Try to push without blocking
    ///
    /// Returns `Err(WouldBlock)` if adding would exceed capacity.
    pub fn try_push(&self, item: T, size_bytes: usize) -> Result<(), TryPushError> {
        let mut inner = self.inner.lock().unwrap();

        if inner.closed {
            return Err(TryPushError::Closed);
        }

        if inner.current_size + size_bytes > self.capacity_bytes {
            return Err(TryPushError::WouldBlock);
        }

        // Add item
        inner.items.push_back((item, size_bytes));
        inner.current_size += size_bytes;

        // Signal that queue is not empty
        self.not_empty.notify_one();

        Ok(())
    }

    /// Pull an item from the queue
    ///
    /// **BLOCKS** if queue is empty.
    /// Returns `None` if queue is closed and empty.
    ///
    /// # Example
    /// ```no_run
    /// # use ragc_core::MemoryBoundedQueue;
    /// # let queue: MemoryBoundedQueue<Vec<u8>> = MemoryBoundedQueue::new(1024);
    /// while let Some(item) = queue.pull() {
    ///     // process(item);
    /// }
    /// // Queue is closed and empty - we're done!
    /// ```
    pub fn pull(&self) -> Option<T> {
        let mut inner = self.inner.lock().unwrap();

        // Wait while queue is empty and not closed
        while inner.items.is_empty() && !inner.closed {
            inner = self.not_empty.wait(inner).unwrap();
        }

        // If closed and empty, return None
        if inner.items.is_empty() {
            return None;
        }

        // Remove item
        let (item, size) = inner.items.pop_front().unwrap();
        inner.current_size -= size;

        // Signal that queue has space
        self.not_full.notify_one();

        Some(item)
    }

    /// Try to pull without blocking
    ///
    /// Returns `None` if queue is empty (even if not closed).
    pub fn try_pull(&self) -> Option<T> {
        let mut inner = self.inner.lock().unwrap();

        if inner.items.is_empty() {
            return None;
        }

        // Remove item
        let (item, size) = inner.items.pop_front().unwrap();
        inner.current_size -= size;

        // Signal that queue has space
        self.not_full.notify_one();

        Some(item)
    }

    /// Close the queue
    ///
    /// After closing:
    /// - No more pushes allowed (returns error)
    /// - Pulls will drain remaining items, then return None
    /// - Workers can detect completion via `pull()` returning None
    pub fn close(&self) {
        let mut inner = self.inner.lock().unwrap();
        inner.closed = true;

        // Wake up all waiting threads
        self.not_full.notify_all();
        self.not_empty.notify_all();
    }

    /// Check if queue is closed
    pub fn is_closed(&self) -> bool {
        self.inner.lock().unwrap().closed
    }

    /// Get current size in bytes
    pub fn current_size(&self) -> usize {
        self.inner.lock().unwrap().current_size
    }

    /// Get current item count
    pub fn len(&self) -> usize {
        self.inner.lock().unwrap().items.len()
    }

    /// Check if queue is empty
    pub fn is_empty(&self) -> bool {
        self.inner.lock().unwrap().items.is_empty()
    }

    /// Get capacity in bytes
    pub fn capacity(&self) -> usize {
        self.capacity_bytes
    }
}

// Make queue cloneable (clones share the same underlying queue)
impl<T> Clone for MemoryBoundedQueue<T> {
    fn clone(&self) -> Self {
        Self {
            inner: Arc::clone(&self.inner),
            capacity_bytes: self.capacity_bytes,
            not_full: Arc::clone(&self.not_full),
            not_empty: Arc::clone(&self.not_empty),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PushError {
    Closed,
}

impl std::fmt::Display for PushError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PushError::Closed => write!(f, "Queue is closed"),
        }
    }
}

impl std::error::Error for PushError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TryPushError {
    Closed,
    WouldBlock,
}

impl std::fmt::Display for TryPushError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TryPushError::Closed => write!(f, "Queue is closed"),
            TryPushError::WouldBlock => write!(f, "Queue is full - would block"),
        }
    }
}

impl std::error::Error for TryPushError {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::thread;
    use std::time::Duration;

    #[test]
    fn test_basic_push_pull() {
        let queue: MemoryBoundedQueue<Vec<u8>> = MemoryBoundedQueue::new(1024);

        // Push some data
        let data = vec![0u8; 100];
        queue.push(data.clone(), 100).unwrap();

        // Pull it back
        let pulled = queue.pull().unwrap();
        assert_eq!(pulled, data);
    }

    #[test]
    fn test_backpressure() {
        let queue: MemoryBoundedQueue<Vec<u8>> = MemoryBoundedQueue::new(1024);

        // Fill queue to capacity
        queue.push(vec![0u8; 512], 512).unwrap();
        queue.push(vec![0u8; 512], 512).unwrap();

        // Try to push more - should block!
        let blocked = Arc::new(AtomicBool::new(false));
        let blocked_clone = Arc::clone(&blocked);
        let queue_clone = queue.clone();

        let handle = thread::spawn(move || {
            blocked_clone.store(true, Ordering::SeqCst);
            queue_clone.push(vec![0u8; 100], 100).unwrap();
            blocked_clone.store(false, Ordering::SeqCst);
        });

        // Wait a bit - thread should still be blocked
        thread::sleep(Duration::from_millis(100));
        assert!(blocked.load(Ordering::SeqCst), "Push should be blocked!");

        // Pull an item - should unblock the push
        queue.pull().unwrap();

        // Wait for thread to finish
        handle.join().unwrap();
        assert!(
            !blocked.load(Ordering::SeqCst),
            "Push should have completed!"
        );
    }

    #[test]
    fn test_close_queue() {
        let queue: MemoryBoundedQueue<Vec<u8>> = MemoryBoundedQueue::new(1024);

        // Push some items
        queue.push(vec![0u8; 100], 100).unwrap();
        queue.push(vec![0u8; 100], 100).unwrap();

        // Close queue
        queue.close();

        // Can't push anymore
        assert!(queue.push(vec![0u8; 100], 100).is_err());

        // Can still pull existing items
        assert!(queue.pull().is_some());
        assert!(queue.pull().is_some());

        // Now empty - returns None
        assert!(queue.pull().is_none());
    }

    #[test]
    fn test_try_operations() {
        let queue: MemoryBoundedQueue<Vec<u8>> = MemoryBoundedQueue::new(100);

        // try_push succeeds when there's space
        assert!(queue.try_push(vec![0u8; 50], 50).is_ok());

        // try_push fails when would exceed capacity
        assert_eq!(
            queue.try_push(vec![0u8; 60], 60),
            Err(TryPushError::WouldBlock)
        );

        // try_pull succeeds when there's an item
        assert!(queue.try_pull().is_some());

        // try_pull returns None when empty
        assert!(queue.try_pull().is_none());
    }

    #[test]
    fn test_multiple_producers_consumers() {
        let queue: MemoryBoundedQueue<usize> = MemoryBoundedQueue::new(1000);

        // Spawn 3 producers
        let mut producers = vec![];
        for i in 0..3 {
            let q = queue.clone();
            producers.push(thread::spawn(move || {
                for j in 0..100 {
                    q.push(i * 100 + j, 10).unwrap();
                }
            }));
        }

        // Spawn 2 consumers
        let mut consumers = vec![];
        for _ in 0..2 {
            let q = queue.clone();
            consumers.push(thread::spawn(move || {
                let mut count = 0;
                while let Some(_) = q.pull() {
                    count += 1;
                    if count == 150 {
                        break; // Each consumer gets 150 items
                    }
                }
                count
            }));
        }

        // Wait for producers
        for p in producers {
            p.join().unwrap();
        }

        // Close queue
        queue.close();

        // Wait for consumers
        let mut total = 0;
        for c in consumers {
            total += c.join().unwrap();
        }

        // Should have consumed all 300 items
        assert_eq!(total, 300);
    }

    #[test]
    fn test_size_tracking() {
        let queue: MemoryBoundedQueue<Vec<u8>> = MemoryBoundedQueue::new(1024);

        assert_eq!(queue.current_size(), 0);

        queue.push(vec![0u8; 100], 100).unwrap();
        assert_eq!(queue.current_size(), 100);

        queue.push(vec![0u8; 200], 200).unwrap();
        assert_eq!(queue.current_size(), 300);

        queue.pull().unwrap();
        assert_eq!(queue.current_size(), 200);

        queue.pull().unwrap();
        assert_eq!(queue.current_size(), 0);
    }
}
