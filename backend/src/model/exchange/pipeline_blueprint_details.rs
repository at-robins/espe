use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::model::internal::pipeline_blueprint::ContextualisedPipelineBlueprint;

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PipelineBlueprintDetails {
    pub id: String,
    pub name: String,
    pub comment: String,
}

impl From<&Arc<ContextualisedPipelineBlueprint>> for PipelineBlueprintDetails {
    fn from(value: &Arc<ContextualisedPipelineBlueprint>) -> Self {
        let inner = value.pipeline();
        PipelineBlueprintDetails {
            id: inner.id().clone(),
            name: inner.name().clone(),
            comment: inner.description().clone(),
        }
    }
}
