use std::path::PathBuf;

use serde::{Serialize, Deserialize};

pub trait PipelineStep {
    fn output(&self) -> PathBuf;
    //fn input()
}

/// The current status of a [`PipelineStep`].
#[derive(Debug, Clone, Copy, PartialEq, Hash, Serialize, Deserialize)]
pub enum PipelineStepStatus {
    /// The pipeline step is currently executed.
    Running,
    /// The pipeline step is waiting for prerequisites to finish.
    Pending,
    /// The pipeline step was finished successfully.
    Success,
    /// The pipeline step failed.
    Failed,
}

impl std::fmt::Display for PipelineStepStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PipelineStepStatus::Running => write!(f, "Running"),
            PipelineStepStatus::Pending => write!(f, "Pending"),
            PipelineStepStatus::Success => write!(f, "Success"),
            PipelineStepStatus::Failed => write!(f, "Failed"),
        }
    }
}