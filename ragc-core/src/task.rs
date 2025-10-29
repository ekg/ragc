// Task types for streaming compression pipeline
// Matches C++ AGC's contig_processing_stage_t and task_t (agc_compressor.h:654-656)

use ragc_common::Contig;

/// Processing stage for a contig task
///
/// This matches C++ AGC's contig_processing_stage_t enum (agc_compressor.h:654):
/// ```cpp
/// enum class contig_processing_stage_t {
///     unknown, all_contigs, new_splitters, hard_contigs, registration
/// };
/// ```
///
/// Stages:
/// - **AllContigs**: Initial contig processing (normal compression flow)
/// - **Registration**: Synchronization barrier for segment registration
///   - Sent after each sample when adaptive_mode is OFF
///   - Workers arrive at barrier, register segments, then continue
/// - **NewSplitters**: Synchronization barrier for adaptive splitter finding
///   - Sent after each sample when adaptive_mode is ON
///   - Workers arrive at barrier, find new splitters, then continue
/// - **HardContigs**: Re-process contigs that need new splitters
///   - Contigs that didn't segment well are re-enqueued with this stage
///   - Uses expanded splitter set from adaptive mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContigProcessingStage {
    /// Initial contig processing (compression)
    AllContigs,

    /// Synchronization barrier: register segments and update terminators
    Registration,

    /// Synchronization barrier: find adaptive splitters for hard contigs
    NewSplitters,

    /// Re-process contigs with new splitters (adaptive mode)
    HardContigs,
}

/// A task for the worker queue
///
/// This matches C++ AGC's task_t typedef (agc_compressor.h:656):
/// ```cpp
/// using task_t = tuple<contig_processing_stage_t, string, string, contig_t>;
/// ```
///
/// Fields:
/// - `stage`: What type of processing this task needs
/// - `sample_name`: Which sample this contig belongs to
/// - `contig_name`: Contig identifier (chromosome name, etc.)
/// - `sequence`: The actual sequence data
///
/// For synchronization tasks (Registration, NewSplitters):
/// - `sample_name`, `contig_name`, and `sequence` are empty
/// - Workers detect sync tasks and enter barrier synchronization
#[derive(Debug, Clone)]
pub struct Task {
    pub stage: ContigProcessingStage,
    pub sample_name: String,
    pub contig_name: String,
    pub sequence: Contig,
}

impl Task {
    /// Create a new contig processing task
    pub fn new_contig(
        sample_name: String,
        contig_name: String,
        sequence: Contig,
        stage: ContigProcessingStage,
    ) -> Self {
        Task {
            stage,
            sample_name,
            contig_name,
            sequence,
        }
    }

    /// Create a synchronization barrier task
    ///
    /// This is sent N times (once per worker) to trigger barrier synchronization.
    /// Matches C++ AGC's EmplaceManyNoCost() pattern (agc_compressor.cpp:2197).
    pub fn new_sync(stage: ContigProcessingStage) -> Self {
        assert!(
            stage == ContigProcessingStage::Registration
                || stage == ContigProcessingStage::NewSplitters,
            "Sync tasks must be Registration or NewSplitters"
        );

        Task {
            stage,
            sample_name: String::new(),
            contig_name: String::new(),
            sequence: Vec::new(),
        }
    }

    /// Check if this is a synchronization task
    pub fn is_sync(&self) -> bool {
        (self.stage == ContigProcessingStage::Registration
            || self.stage == ContigProcessingStage::NewSplitters)
            && self.sample_name.is_empty()
            && self.contig_name.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contig_task_creation() {
        let task = Task::new_contig(
            "sample1".to_string(),
            "chr1".to_string(),
            vec![0, 1, 2, 3],
            ContigProcessingStage::AllContigs,
        );

        assert_eq!(task.stage, ContigProcessingStage::AllContigs);
        assert_eq!(task.sample_name, "sample1");
        assert_eq!(task.contig_name, "chr1");
        assert_eq!(task.sequence, vec![0, 1, 2, 3]);
        assert!(!task.is_sync());
    }

    #[test]
    fn test_registration_sync_task() {
        let task = Task::new_sync(ContigProcessingStage::Registration);

        assert_eq!(task.stage, ContigProcessingStage::Registration);
        assert!(task.sample_name.is_empty());
        assert!(task.contig_name.is_empty());
        assert!(task.sequence.is_empty());
        assert!(task.is_sync());
    }

    #[test]
    fn test_new_splitters_sync_task() {
        let task = Task::new_sync(ContigProcessingStage::NewSplitters);

        assert_eq!(task.stage, ContigProcessingStage::NewSplitters);
        assert!(task.is_sync());
    }

    #[test]
    #[should_panic(expected = "Sync tasks must be Registration or NewSplitters")]
    fn test_invalid_sync_task() {
        Task::new_sync(ContigProcessingStage::AllContigs);
    }

    #[test]
    fn test_hard_contig_task() {
        let task = Task::new_contig(
            "sample1".to_string(),
            "chr1".to_string(),
            vec![0, 1, 2, 3],
            ContigProcessingStage::HardContigs,
        );

        assert_eq!(task.stage, ContigProcessingStage::HardContigs);
        assert!(!task.is_sync());
    }
}
