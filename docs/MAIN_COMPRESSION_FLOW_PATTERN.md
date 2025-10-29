# Main Compression Flow Integration Pattern (C++ AGC)

This document analyzes C++ AGC's main compression flow to guide Rust implementation.

**Purpose**: Understand how C++ AGC orchestrates the entire compression pipeline from sample files to archive completion.

**Key Files**:
- `agc_compressor.cpp` lines 2121-2270 (AddSampleFiles)
- `agc_compressor.cpp` lines 2134-2141 (thread management)
- `agc_compressor.cpp` lines 1093-1270 (start_compressing_threads)

---

## Overview

The main compression flow in C++ AGC follows a producer-consumer pattern:

1. **Main Thread (Producer)**:
   - Opens sample files sequentially
   - Reads contigs
   - Enqueues tasks to priority queue
   - Sends synchronization tokens after each sample

2. **Worker Threads (Consumers)**:
   - Pull tasks from priority queue
   - Compress contigs (inline segmentation)
   - Handle synchronization barriers
   - Store segments to archive

**Critical Pattern**: The main thread streams contigs into the queue while workers process them in parallel. Priority ensures earlier samples are processed first, establishing them as reference.

---

## 1. AddSampleFiles() - Main Entry Point

**Source**: `agc_compressor.cpp` lines 2121-2270

```cpp
bool CAGCCompressor::AddSampleFiles(vector<pair<string, string>> _v_sample_file_name, const uint32_t no_threads)
{
    if (_v_sample_file_name.empty())
        return true;

    processed_bases = 0;

    // Queue setup
    size_t queue_capacity = max(2ull << 30, no_threads * (192ull << 20));
    pq_contigs_desc = make_shared<CBoundedPQueue<task_t>>(1, queue_capacity);
    pq_contigs_desc_aux = make_shared<CBoundedPQueue<task_t>>(1, ~0ull);
    pq_contigs_desc_working = pq_contigs_desc;

    // Worker thread count
    uint32_t no_workers = (no_threads < 8) ? no_threads : no_threads - 1;

    // Spawn worker threads
    vector<thread> v_threads;
    v_threads.reserve((size_t)no_workers);
    my_barrier bar(no_workers);
    start_compressing_threads(v_threads, bar, no_workers);

    // Main processing loop
    pair<string, string> sf;
    CGenomeIO gio;
    string id;
    contig_t contig;
    size_t sample_priority = ~0ull;
    size_t cnt_contigs_in_sample = 0;
    const size_t max_no_contigs_before_synchronization = pack_cardinality;
    const size_t min_size_before_synchronization = 1ull << 30;

    if (archive_version >= 3000 && in_archive != nullptr)
        processed_samples = (uint32_t) dynamic_pointer_cast<CCollection_V3>(collection_desc)->get_no_samples();
    else
        processed_samples = 0;

    if (concatenated_genomes)
        cnt_contigs_in_sample = processed_samples % pack_cardinality;

    size_t num_empty_input = 0;

    for(auto sf : _v_sample_file_name)
    {
        if (archive_version >= 3000)
            dynamic_pointer_cast<CCollection_V3>(collection_desc)->reset_prev_sample_name();

        if (!gio.Open(sf.second, false))
        {
            cerr << "Cannot open file: " << sf.second << endl;
            continue;
        }

        bool any_contigs_read = false;
        bool any_contigs_added = false;

        while (gio.ReadContigRaw(id, contig))
        {
            if (concatenated_genomes)
            {
                // Concatenated mode: treat all contigs as one "sample"
                if (!collection_desc->register_sample_contig("", id))
                    cerr << "Error: Pair sample_name:contig_name " << id << ":" << id << " is already in the archive!\n";
                else
                {
                    auto cost = contig.size();
                    pq_contigs_desc->Emplace(make_tuple(contig_processing_stage_t::all_contigs, "", id, move(contig)), sample_priority, cost);
                    contig.clear();

                    if (++cnt_contigs_in_sample >= max_no_contigs_before_synchronization)
                    {
                        // Send synchronization tokens
                        pq_contigs_desc->EmplaceManyNoCost(make_tuple(
                            adaptive_compression ? contig_processing_stage_t::new_splitters : contig_processing_stage_t::registration,
                            "", "", contig_t()), sample_priority, no_workers);

                        cnt_contigs_in_sample = 0;
                        --sample_priority;
                    }

                    any_contigs_added = true;
                }
            }
            else
            {
                // Normal mode: one sample per file
                if (collection_desc->register_sample_contig(sf.first, id))
                {
                    auto cost = contig.size();
                    pq_contigs_desc->Emplace(make_tuple(contig_processing_stage_t::all_contigs, sf.first, id, move(contig)), sample_priority, cost);
                    contig.clear();
                    any_contigs_added = true;
                }
                else
                    cerr << "Error: Pair sample_name:contig_name " << sf.first << ":" << id << " is already in the archive!\n";
            }

            any_contigs_read = true;
        }

        if (!any_contigs_read)
            cerr << "Warning: Pair sample_name:file_path " << sf.first << ":" << sf.second << " contains no contigs and will not be included in the archive!\n";

        if (!any_contigs_added)
        {
            cerr << "Warning: Pair sample_name:file_path " << sf.first << ":" << sf.second << " contains only contigs already present in the archive!\n";
            ++num_empty_input;
        }

        if (!concatenated_genomes && any_contigs_added)
        {
            // Send synchronization tokens after each sample
            pq_contigs_desc->EmplaceManyNoCost(make_tuple(
                adaptive_compression ? contig_processing_stage_t::new_splitters : contig_processing_stage_t::registration,
                "", "", contig_t()), sample_priority, no_workers);

            --sample_priority;
        }

        gio.Close();
    }

    // Final synchronization for concatenated mode
    if (concatenated_genomes)
    {
        pq_contigs_desc->EmplaceManyNoCost(make_tuple(
            adaptive_compression ? contig_processing_stage_t::new_splitters : contig_processing_stage_t::registration,
            "", "", contig_t()), sample_priority, no_workers);

        cnt_contigs_in_sample = 0;
        --sample_priority;
    }

    // Signal completion and wait for workers
    pq_contigs_desc->MarkCompleted();
    join_threads(v_threads);

    // Final cleanup
    if(concatenated_genomes)
        processed_samples = (uint32_t) dynamic_pointer_cast<CCollection_V3>(collection_desc)->get_no_samples();

    if (archive_version >= 3000 && processed_samples % pack_cardinality != 0)
        dynamic_pointer_cast<CCollection_V3>(collection_desc)->store_contig_batch((processed_samples / pack_cardinality) * pack_cardinality, processed_samples);

    out_archive->FlushOutBuffers();

    pq_contigs_desc.reset();
    pq_contigs_desc_aux.reset();
    pq_contigs_desc_working.reset();

    no_samples_in_archive += _v_sample_file_name.size() - num_empty_input;

    return true;
}
```

---

## 2. Step-by-Step Analysis

### Step 1: Queue Initialization (lines 2128-2132)

```cpp
size_t queue_capacity = max(2ull << 30, no_threads * (192ull << 20));
pq_contigs_desc = make_shared<CBoundedPQueue<task_t>>(1, queue_capacity);
pq_contigs_desc_aux = make_shared<CBoundedPQueue<task_t>>(1, ~0ull);
pq_contigs_desc_working = pq_contigs_desc;
```

**Main Queue (`pq_contigs_desc`)**:
- Capacity: `max(2GB, num_threads * 192MB)`
- Used for normal contig processing (all_contigs stage)
- Bounded by memory to prevent main thread from reading too far ahead

**Auxiliary Queue (`pq_contigs_desc_aux`)**:
- Capacity: unlimited (`~0ull`)
- Used for hard contigs reprocessing (hard_contigs stage)
- Only used in adaptive mode

**Working Queue Pointer (`pq_contigs_desc_working`)**:
- Points to currently active queue
- Starts with main queue
- Switches to aux queue in new_splitters stage (adaptive mode only)

**Rationale**: The bounded main queue prevents memory explosion when reading large genomes. Workers pull from the queue, creating backpressure on the main thread if queue fills up.

### Step 2: Worker Thread Count (line 2134)

```cpp
uint32_t no_workers = (no_threads < 8) ? no_threads : no_threads - 1;
```

**Logic**:
- If `num_threads < 8`: Use ALL threads as workers
- If `num_threads >= 8`: Reserve 1 thread for main, use `num_threads - 1` as workers

**Rationale**: For low thread counts (1-7), every thread is precious. For high thread counts (8+), reserving 1 thread for the main loop ensures smooth I/O without blocking workers.

**Example**:
- 1 thread: 1 worker (main thread blocks during I/O)
- 6 threads: 6 workers
- 8 threads: 7 workers (1 reserved for main)
- 16 threads: 15 workers (1 reserved for main)

### Step 3: Thread Spawning (lines 2136-2141)

```cpp
vector<thread> v_threads;
v_threads.reserve((size_t)no_workers);

my_barrier bar(no_workers);

start_compressing_threads(v_threads, bar, no_workers);
```

**Barrier Initialization**:
- Created with `no_workers` threads
- All workers must arrive at barrier before any continue
- Used for synchronization at registration and new_splitters stages

**Worker Thread Spawning**:
- `start_compressing_threads()` creates `no_workers` threads
- Each runs the worker loop (see Phase 2 documentation)
- All share the same barrier

**Critical**: Worker threads start immediately and begin pulling from the queue. Main thread will start enqueuing tasks next.

### Step 4: Variables Initialization (lines 2143-2161)

```cpp
pair<string, string> sf;
CGenomeIO gio;
string id;
contig_t contig;
size_t sample_priority = ~0ull;
size_t cnt_contigs_in_sample = 0;
const size_t max_no_contigs_before_synchronization = pack_cardinality;
const size_t min_size_before_synchronization = 1ull << 30;

if (archive_version >= 3000 && in_archive != nullptr)
    processed_samples = (uint32_t) dynamic_pointer_cast<CCollection_V3>(collection_desc)->get_no_samples();
else
    processed_samples = 0;

if (concatenated_genomes)
    cnt_contigs_in_sample = processed_samples % pack_cardinality;

size_t num_empty_input = 0;
```

**Key Variables**:
- `sample_priority`: Starts at `~0ull` (max value), decrements after each sample
  - **Why**: Ensures earlier samples have higher priority (processed first)
  - Priority queue pops highest priority first
  - Establishes earlier samples as reference for later samples
- `cnt_contigs_in_sample`: Tracks contigs processed in concatenated mode
- `max_no_contigs_before_synchronization`: Sync every `pack_cardinality` contigs (concatenated mode)
- `processed_samples`: Tracks total samples processed (for appending to existing archive)

**Concatenated Genomes Mode**:
- Treats ALL contigs as one logical "sample"
- Syncs every `pack_cardinality` contigs (e.g., every 128 contigs)
- Used for pangenomes where individual sample identity doesn't matter

### Step 5: Main Processing Loop (lines 2163-2242)

```cpp
for(auto sf : _v_sample_file_name)
{
    if (archive_version >= 3000)
        dynamic_pointer_cast<CCollection_V3>(collection_desc)->reset_prev_sample_name();

    if (!gio.Open(sf.second, false))
    {
        cerr << "Cannot open file: " << sf.second << endl;
        continue;
    }

    bool any_contigs_read = false;
    bool any_contigs_added = false;

    while (gio.ReadContigRaw(id, contig))
    {
        // Mode-specific processing...
    }

    // Post-sample sync...
}
```

**Flow**:
1. Reset previous sample name (archive version 3+)
2. Open FASTA file
3. Read contigs one-by-one
4. Enqueue tasks for workers
5. Send sync tokens after sample
6. Close file
7. Move to next sample

#### Normal Mode (not concatenated_genomes)

```cpp
if (collection_desc->register_sample_contig(sf.first, id))
{
    auto cost = contig.size();
    pq_contigs_desc->Emplace(make_tuple(contig_processing_stage_t::all_contigs, sf.first, id, move(contig)), sample_priority, cost);
    contig.clear();
    any_contigs_added = true;
}
else
    cerr << "Error: Pair sample_name:contig_name " << sf.first << ":" << id << " is already in the archive!\n";
```

**Steps**:
1. **Register contig**: `register_sample_contig(sample_name, contig_name)`
   - Returns false if (sample, contig) already exists in archive
   - Prevents duplicate contigs
2. **Calculate cost**: `contig.size()` (number of bases)
   - Used for queue capacity management
3. **Enqueue task**: `Emplace(task, priority, cost)`
   - Task: `(all_contigs, sample_name, contig_name, sequence)`
   - Priority: `sample_priority` (decrements after each sample)
   - Cost: contig size
4. **Clear contig**: `move(contig)` transfers ownership, `contig.clear()` frees local memory

**Post-Sample Sync**:
```cpp
if (!concatenated_genomes && any_contigs_added)
{
    // Send synchronization tokens
    pq_contigs_desc->EmplaceManyNoCost(make_tuple(
        adaptive_compression ? contig_processing_stage_t::new_splitters : contig_processing_stage_t::registration,
        "", "", contig_t()), sample_priority, no_workers);

    --sample_priority;
}
```
- After each sample, enqueue `no_workers` sync tokens
- **Sync stage**: `new_splitters` (adaptive) or `registration` (normal)
- **Why `no_workers` tokens**: Each worker must receive one token to trigger barrier
- Decrement `sample_priority` so next sample has lower priority

#### Concatenated Mode (concatenated_genomes)

```cpp
if (!collection_desc->register_sample_contig("", id))
    cerr << "Error: Pair sample_name:contig_name " << id << ":" << id << " is already in the archive!\n";
else
{
    auto cost = contig.size();
    pq_contigs_desc->Emplace(make_tuple(contig_processing_stage_t::all_contigs, "", id, move(contig)), sample_priority, cost);
    contig.clear();

    if (++cnt_contigs_in_sample >= max_no_contigs_before_synchronization)
    {
        // Send synchronization tokens
        pq_contigs_desc->EmplaceManyNoCost(make_tuple(
            adaptive_compression ? contig_processing_stage_t::new_splitters : contig_processing_stage_t::registration,
            "", "", contig_t()), sample_priority, no_workers);

        cnt_contigs_in_sample = 0;
        --sample_priority;
    }

    any_contigs_added = true;
}
```

**Differences from Normal Mode**:
1. Sample name is empty string `""`
2. Sync every `pack_cardinality` contigs (e.g., 128), not every file
3. Reset `cnt_contigs_in_sample` after sync

**Rationale**: In pangenome mode, individual sample identity doesn't matter. All contigs are treated as one logical collection. Syncing every 128 contigs (or similar) balances memory usage with synchronization overhead.

### Step 6: Completion (lines 2244-2270)

```cpp
// Final synchronization for concatenated mode
if (concatenated_genomes)
{
    pq_contigs_desc->EmplaceManyNoCost(make_tuple(
        adaptive_compression ? contig_processing_stage_t::new_splitters : contig_processing_stage_t::registration,
        "", "", contig_t()), sample_priority, no_workers);

    cnt_contigs_in_sample = 0;
    --sample_priority;
}

// Signal completion and wait for workers
pq_contigs_desc->MarkCompleted();
join_threads(v_threads);

// Final cleanup
if(concatenated_genomes)
    processed_samples = (uint32_t) dynamic_pointer_cast<CCollection_V3>(collection_desc)->get_no_samples();

if (archive_version >= 3000 && processed_samples % pack_cardinality != 0)
    dynamic_pointer_cast<CCollection_V3>(collection_desc)->store_contig_batch((processed_samples / pack_cardinality) * pack_cardinality, processed_samples);

out_archive->FlushOutBuffers();

pq_contigs_desc.reset();
pq_contigs_desc_aux.reset();
pq_contigs_desc_working.reset();

no_samples_in_archive += _v_sample_file_name.size() - num_empty_input;

return true;
```

**Final Sync (Concatenated Mode)**:
- Send one last sync token batch to flush remaining buffered segments

**Mark Completed**:
- `pq_contigs_desc->MarkCompleted()`: Signal no more tasks coming
- Workers will exit their loops after processing all tasks

**Join Threads**:
- `join_threads(v_threads)`: Wait for all workers to complete
- Blocks until all workers exit

**Batch Storage** (Archive Version 3+):
- Store any remaining contigs in the last incomplete batch
- Batches are `pack_cardinality` contigs (e.g., 128)

**Flush Buffers**:
- Write all pending archive data to disk

**Cleanup**:
- Reset queue shared pointers (free memory)
- Update archive sample count

---

## 3. Priority and Ordering

### Sample Priority

```cpp
size_t sample_priority = ~0ull;  // Start at max

// After each sample:
--sample_priority;
```

**Ordering**:
```
Sample 1: priority = 18446744073709551615
Sample 2: priority = 18446744073709551614
Sample 3: priority = 18446744073709551613
...
```

**Why High-to-Low**:
- Priority queue pops **highest priority** first
- Earlier samples processed before later samples
- Establishes Sample 1 as "reference" for differential encoding

### Task Cost

```cpp
auto cost = contig.size();
pq_contigs_desc->Emplace(..., sample_priority, cost);
```

**Purpose**: Queue capacity management
- Bounded queue has capacity in bytes (e.g., 2GB)
- Each task has cost (contig size)
- When total cost exceeds capacity, `Emplace()` blocks until workers free space

**Example**:
```
Queue capacity: 2GB
Contig A: 100MB → enqueued, total cost = 100MB
Contig B: 200MB → enqueued, total cost = 300MB
...
Contig X: 50MB → would exceed 2GB, blocks until workers process some contigs
```

**Benefit**: Prevents main thread from reading entire genome into memory before workers start processing.

---

## 4. Synchronization Token Pattern

### Normal Mode

```cpp
// After each sample file:
pq_contigs_desc->EmplaceManyNoCost(
    make_tuple(
        adaptive_compression ? new_splitters : registration,
        "", "", contig_t()
    ),
    sample_priority,
    no_workers  // Send N tokens for N workers
);
```

**Flow**:
```
Sample 1 contigs → [Worker 1, Worker 2, Worker 3 process]
Sample 1 sync tokens (3x) → [Each worker receives one, triggers barrier]
                          → [All arrive at barrier]
                          → [Leader: register_segments]
                          → [All: store_segments]
                          → [Continue to next sample]
Sample 2 contigs → ...
```

**Critical**: `no_workers` sync tokens ensures each worker receives exactly one. All workers hit the barrier simultaneously.

### Concatenated Mode

```cpp
// Every pack_cardinality contigs:
if (++cnt_contigs_in_sample >= max_no_contigs_before_synchronization)
{
    pq_contigs_desc->EmplaceManyNoCost(..., no_workers);
    cnt_contigs_in_sample = 0;
    --sample_priority;
}
```

**Flow**:
```
Contigs 1-128 → [Workers process]
Sync tokens → [Barrier, register, store]
Contigs 129-256 → [Workers process]
Sync tokens → [Barrier, register, store]
...
```

---

## 5. Adaptive Mode vs Normal Mode

### Normal Mode

```cpp
pq_contigs_desc->EmplaceManyNoCost(
    make_tuple(registration, "", "", contig_t()),
    sample_priority, no_workers
);
```
- Sync stage: **registration**
- Workers execute registration stage handler (4 barriers)
- No adaptive splitter finding

### Adaptive Mode

```cpp
pq_contigs_desc->EmplaceManyNoCost(
    make_tuple(new_splitters, "", "", contig_t()),
    sample_priority, no_workers
);
```
- Sync stage: **new_splitters**
- Workers execute new_splitters stage handler (2 barriers)
- Hard contigs re-enqueued as **hard_contigs** stage
- After new_splitters, additional **registration** sync tokens sent
- Workers process hard contigs with new splitters
- Then hit final registration barrier

**Full Adaptive Flow**:
```
Sample 1:
  Contigs → workers compress (some fail, saved to v_raw_contigs)
  new_splitters sync → find new splitters, re-enqueue hard contigs
  registration sync → register/store segments from hard contigs
Sample 2:
  ... (repeat)
```

---

## 6. Queue Switching (Adaptive Mode Only)

**In `new_splitters` handler** (thread 0):
```cpp
// Re-enqueue hard contigs to aux queue
for (auto& x : v_raw_contigs)
{
    pq_contigs_desc_aux->EmplaceNoLock(
        make_tuple(hard_contigs, get<0>(x), get<1>(x), move(get<2>(x))),
        1, cost
    );
}

// Enqueue registration sync to aux queue
pq_contigs_desc_aux->EmplaceManyNoCost(
    make_tuple(registration, "", "", contig_t()),
    0, n_t
);

// Switch working queue
pq_contigs_desc_working = pq_contigs_desc_aux;
```

**Rationale**:
- Main queue (`pq_contigs_desc`) still has unprocessed tasks from upcoming samples
- Hard contigs need immediate reprocessing with new splitters
- Aux queue (`pq_contigs_desc_aux`) provides isolated processing
- After hard contigs complete, switch back to main queue

**Why Unlimited Capacity for Aux Queue**:
- Hard contigs are already in memory (in `v_raw_contigs`)
- No new I/O, so no need to limit memory

---

## 7. Rust Implementation Strategy

### Main Compressor Structure

```rust
pub struct StreamingCompressor {
    config: CompressorConfig,
    archive: Arc<Mutex<Archive>>,
    collection: Arc<Mutex<Collection>>,
    shared: Arc<SharedCompressorState>,

    // Not created until add_samples_streaming() called
    task_queue: Option<Arc<BoundedPriorityQueue<Task>>>,
    task_queue_aux: Option<Arc<BoundedPriorityQueue<Task>>>,
}

impl StreamingCompressor {
    pub fn add_samples_streaming(
        &mut self,
        sample_files: &[(String, PathBuf)],
    ) -> Result<()> {
        // Step 1: Initialize queues
        let queue_capacity = (2usize << 30).max(self.config.num_threads * (192usize << 20));
        let main_queue = Arc::new(BoundedPriorityQueue::new(queue_capacity));
        let aux_queue = Arc::new(BoundedPriorityQueue::new(usize::MAX));

        self.task_queue = Some(Arc::clone(&main_queue));
        self.task_queue_aux = Some(Arc::clone(&aux_queue));

        // Update shared state with queue reference
        self.shared.set_working_queue(Arc::clone(&main_queue));

        // Step 2: Determine worker count
        let num_workers = if self.config.num_threads < 8 {
            self.config.num_threads
        } else {
            self.config.num_threads - 1
        };

        // Step 3: Spawn worker threads
        let barrier = Arc::new(Barrier::new(num_workers));
        let handles = self.spawn_workers(num_workers, Arc::clone(&barrier))?;

        // Step 4: Main processing loop
        let mut sample_priority = usize::MAX;

        for (sample_name, path) in sample_files {
            // Reset prev sample name (archive v3+)
            if self.archive_version >= 3000 {
                self.collection.lock().unwrap().reset_prev_sample_name();
            }

            // Open file
            let mut reader = GenomeIO::open(path)?;

            let mut any_contigs_added = false;

            // Read contigs
            while let Some((contig_name, sequence)) = reader.read_contig()? {
                // Register contig
                if !self.collection.lock().unwrap()
                    .register_sample_contig(sample_name, &contig_name)?
                {
                    eprintln!("Error: Pair {}:{} already exists", sample_name, contig_name);
                    continue;
                }

                // Calculate cost
                let cost = sequence.len();

                // Enqueue task
                main_queue.emplace(
                    Task {
                        stage: ContigProcessingStage::AllContigs,
                        sample_name: sample_name.clone(),
                        contig_name,
                        sequence,
                    },
                    sample_priority,
                    cost,
                )?;

                any_contigs_added = true;
            }

            // Send sync tokens after sample
            if any_contigs_added {
                let sync_stage = if self.config.adaptive_mode {
                    ContigProcessingStage::NewSplitters
                } else {
                    ContigProcessingStage::Registration
                };

                main_queue.emplace_many_nocost(
                    Task {
                        stage: sync_stage,
                        sample_name: String::new(),
                        contig_name: String::new(),
                        sequence: Vec::new(),
                    },
                    sample_priority,
                    num_workers,
                )?;

                sample_priority -= 1;
            }
        }

        // Step 5: Mark completed and wait
        main_queue.mark_completed();

        for handle in handles {
            handle.join().unwrap();
        }

        // Step 6: Final cleanup
        self.archive.lock().unwrap().flush_out_buffers()?;

        self.task_queue = None;
        self.task_queue_aux = None;

        Ok(())
    }

    fn spawn_workers(
        &self,
        num_workers: usize,
        barrier: Arc<Barrier>,
    ) -> Result<Vec<JoinHandle<()>>> {
        let mut handles = Vec::new();

        for worker_id in 0..num_workers {
            let queue = Arc::clone(self.task_queue.as_ref().unwrap());
            let barrier = Arc::clone(&barrier);
            let shared = Arc::clone(&self.shared);

            let handle = thread::spawn(move || {
                worker_thread(worker_id, num_workers, queue, barrier, shared)
            });

            handles.push(handle);
        }

        Ok(handles)
    }
}
```

### Concatenated Genomes Mode (Future)

```rust
// Inside main loop:
let mut cnt_contigs_in_sample = 0;
let max_no_contigs_before_sync = self.config.pack_cardinality;

for (sample_name, path) in sample_files {
    while let Some((contig_name, sequence)) = reader.read_contig()? {
        // Register with empty sample name
        if self.config.concatenated_genomes {
            self.collection.lock().unwrap()
                .register_sample_contig("", &contig_name)?;
        } else {
            // ... normal mode
        }

        // Enqueue task
        main_queue.emplace(...)?;

        // Sync every pack_cardinality contigs (concatenated mode)
        if self.config.concatenated_genomes {
            cnt_contigs_in_sample += 1;

            if cnt_contigs_in_sample >= max_no_contigs_before_sync {
                main_queue.emplace_many_nocost(..., num_workers)?;
                cnt_contigs_in_sample = 0;
                sample_priority -= 1;
            }
        }
    }

    // Sync after sample (normal mode)
    if !self.config.concatenated_genomes {
        main_queue.emplace_many_nocost(..., num_workers)?;
        sample_priority -= 1;
    }
}

// Final sync for concatenated mode
if self.config.concatenated_genomes && cnt_contigs_in_sample > 0 {
    main_queue.emplace_many_nocost(..., num_workers)?;
}
```

---

## 8. Key Implementation Notes

### Thread Count Decision

- **C++ AGC**: `no_workers = (no_threads < 8) ? no_threads : no_threads - 1`
- **Rust**: Match exactly for compatibility
- **Testing**: Verify memory usage scales correctly with thread count

### Queue Capacity

- **Main queue**: `max(2GB, num_threads * 192MB)`
- **Aux queue**: Unlimited
- **Critical**: Bounded capacity prevents memory explosion

### Priority Ordering

- Start at `usize::MAX`, decrement after each sample
- Earlier samples have higher priority
- Ensures deterministic reference selection

### Sync Token Count

- MUST send exactly `num_workers` sync tokens
- Each worker receives one, triggers barrier
- Miscount causes deadlock

### Queue Switching (Adaptive Mode)

- Switch to aux queue in `new_splitters` handler
- Re-enqueue hard contigs to aux
- Send registration sync to aux
- Workers process from aux until exhausted
- (C++ AGC doesn't explicitly switch back, but aux exhausts first)

---

## 9. Testing Strategy

### Unit Tests

1. **Queue initialization**:
   - Verify capacity calculation
   - Verify aux queue unlimited

2. **Worker count**:
   - Test boundary (threads < 8, threads >= 8)

3. **Priority ordering**:
   - Enqueue tasks with different priorities
   - Verify pop order

### Integration Tests

1. **Single sample**:
   - Verify sync tokens sent after sample
   - Verify worker completion

2. **Multiple samples**:
   - Verify priority ordering (Sample 1 before Sample 2)
   - Verify deterministic output

3. **Adaptive mode** (future):
   - Verify queue switching
   - Verify hard contig reprocessing

---

## 10. Summary

**AddSampleFiles** is the orchestration layer:
- Spawns worker threads
- Streams contigs from files
- Enqueues tasks with priority
- Sends synchronization tokens
- Waits for completion
- Flushes and cleans up

**Critical for correctness**:
- Priority ensures earlier samples are reference
- Bounded queue prevents memory explosion
- Sync token count must match worker count
- Queue switching enables adaptive reprocessing

**Next Steps**:
- Phase 4.2: Worker thread spawning details
- Phase 4.3: Queue completion and cleanup
- Phase 5: Adaptive mode integration
