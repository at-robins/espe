use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
/// The full logs of an experiment step.
pub struct ExperimentStepLogs {
    /// The build logs.
    pub build: ExperimentStepLog,
    /// The run logs.
    pub run: ExperimentStepLog,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
/// The logs of a specific process of an experiment step.
pub struct ExperimentStepLog {
    /// The stdout output.
    pub stdout: Option<String>,
    /// The stderr output.
    pub stderr: Option<String>,
    /// The process exit code.
    pub exit_code: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
/// Information required to request log file
/// for a specific pipeline step.
pub struct ExperimentStepLogRequest {
    /// The ID of the pipeline.
    pub pipeline_id: String,
    /// The ID of the pipeline step.
    pub step_id: String,
}
