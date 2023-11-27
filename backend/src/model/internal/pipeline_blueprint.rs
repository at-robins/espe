use std::path::PathBuf;

use getset::Getters;
use serde::{Deserialize, Serialize};

use crate::application::error::SeqError;

/// The pipeline sub-directory that imports are stored in.
const IMPORT_DIRECOTRY: &str = "imports";
/// The tag marking the start of an import.
const IMPORT_START_TAG: &str = "<espe-import>";
/// The tag marking the end of an import.
const IMPORT_END_TAG: &str = "</espe-import>";

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
    /// The variables that can be specified for the pipeline step.
    #[getset(get = "pub")]
    variables: Vec<PipelineStepVariable>,
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

    /// Returns the context path of the imports directory.
    pub fn imports_path(&self) -> PathBuf {
        self.context.join(IMPORT_DIRECOTRY)
    }

    /// Resolves all imports specified within the pipeline blueprint.
    pub fn resolve_imports(&mut self) -> Result<(), SeqError> {
        let imports_path = self.imports_path();
        if imports_path.is_dir() {
            log::info!("Resolving pipeline imports for {}", imports_path.display());
            // Resolve pipeline description imports.
            if let Some(import_content) =
                Self::load_description_import(self.pipeline.description(), &imports_path)?
            {
                self.pipeline.description = import_content;
            }
            // Resolve pipeline step description imports.
            for step_index in 0..self.pipeline().steps().len() {
                if let Some(import_content) = Self::load_description_import(
                    self.pipeline().steps()[step_index].description(),
                    &imports_path,
                )? {
                    self.pipeline.steps[step_index].description = import_content;
                }
                // Resolve pipeline step variable description imports.
                for variable_index in 0..self.pipeline().steps()[step_index].variables().len() {
                    if let Some(import_content) = Self::load_description_import(
                        self.pipeline().steps()[step_index].variables()[variable_index]
                            .description(),
                        &imports_path,
                    )? {
                        self.pipeline.steps[step_index].variables[variable_index].description =
                            import_content;
                    }
                }
            }
        }
        Ok(())
    }

    /// Checks the specified description for an import statement and
    /// loads the import file's content if specified.
    ///
    /// # Parameters
    ///
    /// * `description` - the description to check for an import statement
    /// * `imports_path` - the context path for pipeline imports
    fn load_description_import(
        description: &String,
        imports_path: &PathBuf,
    ) -> Result<Option<String>, SeqError> {
        let trimmed = description.trim();
        let trimmed_lowercase = trimmed.to_lowercase();
        if trimmed_lowercase.starts_with(&IMPORT_START_TAG.to_lowercase())
            && trimmed_lowercase.ends_with(&IMPORT_END_TAG.to_lowercase())
        {
            // Appends the import file specified in the import tag
            // to the pipeline specific import path.
            let import_path = imports_path.join(
                trimmed[IMPORT_START_TAG.len()..trimmed.len() - IMPORT_END_TAG.len()]
                    .trim()
                    .to_string(),
            );
            // If the import is valid, return the import file's content.
            if import_path.is_file() {
                return Ok(Some(std::fs::read_to_string(import_path)?));
            }
        }
        Ok(None)
    }
}

/// The definition of a pipeline step variable.
#[derive(Debug, Clone, Getters, PartialEq, Serialize, Deserialize)]
pub struct PipelineStepVariable {
    /// The unique ID of the variable.
    #[getset(get = "pub")]
    id: String,
    /// The name for display.
    #[getset(get = "pub")]
    name: String,
    /// A description of the variable.
    #[getset(get = "pub")]
    description: String,
    /// The type of variable.
    #[getset(get = "pub")]
    category: PipelineStepVariableCategory,
    /// If the variable must be set by its instance.
    #[getset(get = "pub")]
    required: Option<bool>,
}

/// The definition of a pipeline step variable type.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "tag", content = "content")]
pub enum PipelineStepVariableCategory {
    /// A boolean checkbox.
    Boolean,
    /// A reference to global data.
    Global,
    /// A number field.
    Number,
    /// An option dropdown.
    Option(Vec<PipelineStepVariableCategoryOption>),
    /// A text field.
    String,
}

/// A pipeline step variable option type.
#[derive(Debug, Clone, Getters, PartialEq, Serialize, Deserialize)]
pub struct PipelineStepVariableCategoryOption {
    /// The name for display.
    #[getset(get = "pub")]
    name: String,
    /// The actual variable name.
    #[getset(get = "pub")]
    value: String,
}
