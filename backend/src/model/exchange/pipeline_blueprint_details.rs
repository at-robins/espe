use serde::{Deserialize, Serialize};

use crate::model::db::pipeline::Pipeline;

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PipelineBlueprintDetails {
    pub id: i32,
    pub name: String,
    pub comment: String,
}

impl From<&Pipeline> for PipelineBlueprintDetails {
    fn from(value: &Pipeline) -> Self {
        PipelineBlueprintDetails {
            id: value.id,
            name: value.pipeline_name.clone(),
            comment: value.comment.clone(),
        }
    }
}
