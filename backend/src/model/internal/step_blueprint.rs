use std::path::PathBuf;

use getset::Getters;
use serde::{Deserialize, Serialize};

/// The definition of a pipeline step.
#[derive(Debug, Clone, Getters, PartialEq, Serialize, Deserialize)]
pub struct PipelineStepBlueprint {
    /// The unique ID of the pipeline step.
    #[getset(get = "pub")]
    id: String,
    /// The name for display.
    #[getset(get = "pub")]
    name: String,
    /// A description of the pipeline step.
    #[getset(get = "pub")]
    description: String,
    /// The container used to execute the pipeline step.
    #[getset(get = "pub")]
    container: String,
    /// The IDs of pipeline steps this step depends on.
    #[getset(get = "pub")]
    dependencies: Vec<String>,
}

/// The definition of a pipeline.
#[derive(Debug, Clone, Getters, PartialEq, Serialize, Deserialize)]
pub struct PipelineBlueprint {
    /// The unique ID of the pipeline.
    #[getset(get = "pub")]
    id: String,
    /// The name for display.
    #[getset(get = "pub")]
    name: String,
    /// A description of the pipeline.
    #[getset(get = "pub")]
    description: String,
    /// The [`PipelineStepBlueprint`] that make up the pipeline.
    #[getset(get = "pub")]
    steps: Vec<PipelineStepBlueprint>,
}

/// The definition of a pipeline in the context of its containing directory.
#[derive(Debug, Clone, Getters, PartialEq)]
pub struct ContextualisedPipelineBlueprint {
    /// The pipeline.
    #[getset(get = "pub")]
    pipeline: PipelineBlueprint,
    /// The context folder of the pipeline.
    #[getset(get = "pub")]
    context: PathBuf,
}

impl ContextualisedPipelineBlueprint {
    /// Creates a new `ContextualisedPipelineBlueprint`
    /// 
    /// # Parameters
    /// 
    /// * `pipeline` - the underlying [`PipelineBlueprint`] definition
    /// * `context` - the directory in which the pipeline definition resides
    pub fn new<P: Into<PathBuf>>(pipeline: PipelineBlueprint, context: P) -> Self {
        ContextualisedPipelineBlueprint {
            pipeline,
            context: context.into(),
        }
    }
}
