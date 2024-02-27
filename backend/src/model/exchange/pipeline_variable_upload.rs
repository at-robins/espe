use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PipelineStepVariableUpload {
    pub pipeline_id: String,
    pub pipeline_step_id: String,
    pub variable_id: String,
    pub variable_value: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PipelineGlobalVariableUpload {
    pub pipeline_id: String,
    pub variable_id: String,
    pub variable_value: Option<String>,
}

