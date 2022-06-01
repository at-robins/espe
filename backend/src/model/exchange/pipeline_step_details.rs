use chrono::{DateTime, Utc};
use serde::{Serialize, Deserialize};

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PipelineStepDetails {
    pub id: i64,
    pub name: String,
    pub status: String,
    pub creation_time: DateTime<Utc>,
}
