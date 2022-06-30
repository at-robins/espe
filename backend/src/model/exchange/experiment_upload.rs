use serde::{Serialize, Deserialize};

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ExperimentUpload {
    pub name: String,
    pub mail: String,
    pub comment: String,
    pub pipeline_id: i32,
}
