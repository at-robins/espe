use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::application::error::{SeqError, SeqErrorType};

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

const STATUS_FAILED: &str = "Failed";
const STATUS_PENDING: &str = "Pending";
const STATUS_RUNNING: &str = "Running";
const STATUS_SUCCESS: &str = "Success";

impl std::fmt::Display for PipelineStepStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PipelineStepStatus::Failed => write!(f, "{}", STATUS_FAILED),
            PipelineStepStatus::Pending => write!(f, "{}", STATUS_PENDING),
            PipelineStepStatus::Running => write!(f, "{}", STATUS_RUNNING),
            PipelineStepStatus::Success => write!(f, "{}", STATUS_SUCCESS),
        }
    }
}

impl TryFrom<&str> for PipelineStepStatus {
    type Error = SeqError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            STATUS_FAILED => Ok(PipelineStepStatus::Failed),
            STATUS_PENDING => Ok(PipelineStepStatus::Pending),
            STATUS_RUNNING => Ok(PipelineStepStatus::Running),
            STATUS_SUCCESS => Ok(PipelineStepStatus::Success),
            _ => Err(SeqError::new(
                "PipelineStepStatus",
                SeqErrorType::InternalServerError,
                format!("{} cannot be converted into a status.", value),
                "An invalid pipeline step was supplied.",
            )),
        }
    }
}
