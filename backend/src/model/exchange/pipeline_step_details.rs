use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};

use crate::model::internal::pipeline_blueprint::PipelineStepVariableCategory;

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PipelineStepDetails {
    pub id: i64,
    pub name: String,
    pub status: String,
    pub creation_time: DateTime<Utc>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PipelineStepVariableInstance {
    pub id: String,
    pub value: String,
    pub category: PipelineStepVariableCategory,
}

impl PipelineStepVariableInstance {
    /// Returns `true` if the variable is an instance of a reference to global data.
    pub fn is_global_data_reference(&self) -> bool {
        self.category.eq(&PipelineStepVariableCategory::Global)
    }
}
