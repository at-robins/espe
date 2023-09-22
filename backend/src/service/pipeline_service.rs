use std::{collections::HashMap, path::PathBuf, sync::Arc};

use actix_web::web;
use parking_lot::Mutex;

use crate::{
    application::{
        config::{Configuration, PIPELINE_DEFINITION_FILE},
        error::{SeqError, SeqErrorType},
    },
    model::internal::pipeline_blueprint::{ContextualisedPipelineBlueprint, PipelineBlueprint},
};

#[derive(Debug)]
/// All pipelines currently loaded by the application.
pub struct LoadedPipelines {
    pipeline_map: Mutex<HashMap<String, Arc<ContextualisedPipelineBlueprint>>>,
}

impl LoadedPipelines {
    /// Loads and stores all pipelines based on the supplied application [`Configuration`].
    ///
    /// # Parameters
    ///
    /// * `app_config` - the app [`Configuration`]
    pub fn new(app_config: web::Data<Configuration>) -> Result<Self, SeqError> {
        Ok(Self {
            pipeline_map: Mutex::new(Self::load_pipeline_map(app_config)?),
        })
    }

    /// Returns the pipeline with the specified ID if loaded.
    ///
    /// # Parameters
    ///
    /// * `id` - the pipeline ID
    pub fn get<T: AsRef<str>>(&self, id: T) -> Option<Arc<ContextualisedPipelineBlueprint>> {
        self.pipeline_map
            .lock()
            .get(id.as_ref())
            .map(|value| Arc::clone(value))
    }

    /// Returns ```true``` if the specified variable exists in the loaded pipelines.
    ///
    /// # Parameters
    ///
    /// * `pipeline_id` - the ID of the pipeline the variable belongs to
    /// * `pipeline_step_id` - the ID of the pipeline step the variable belongs to
    /// * `variable_id` - the ID of the variable
    pub fn has_variable<T: AsRef<str>, R: AsRef<str>, S: AsRef<str>>(
        &self,
        pipeline_id: T,
        pipeline_step_id: R,
        variable_id: S,
    ) -> bool {
        if let Some(pipeline) = self.get(pipeline_id) {
            pipeline
                .pipeline()
                .steps()
                .iter()
                .filter(|step| step.id() == pipeline_step_id.as_ref())
                .any(|step| {
                    step.variables()
                        .iter()
                        .any(|variable| variable.id() == variable_id.as_ref())
                })
        } else {
            false
        }
    }

    /// Returns all loaded pipelines.
    pub fn pipelines(&self) -> Vec<Arc<ContextualisedPipelineBlueprint>> {
        self.pipeline_map
            .lock()
            .values()
            .map(|value| Arc::clone(value))
            .collect()
    }

    /// Returns ```true``` if the pipeline with the specified ID is loaded.
    ///
    /// # Parameters
    ///
    /// * `id` - the pipeline ID
    pub fn is_loaded<T: AsRef<str>>(&self, id: T) -> bool {
        self.pipeline_map.lock().contains_key(id.as_ref())
    }

    /// Updates the currently loaded pipelines.
    ///
    /// # Parameters
    ///
    /// * `app_config` - the app [`Configuration`]
    pub fn update_loaded_pipelines(
        &self,
        app_config: web::Data<Configuration>,
    ) -> Result<(), SeqError> {
        *(self.pipeline_map.lock()) = Self::load_pipeline_map(app_config)?;
        Ok(())
    }

    /// Creates a map of loaded pipelines by their respective ID.
    ///
    /// # Parameters
    ///
    /// * `app_config` - the app [`Configuration`]
    fn load_pipeline_map(
        app_config: web::Data<Configuration>,
    ) -> Result<HashMap<String, Arc<ContextualisedPipelineBlueprint>>, SeqError> {
        let pipelines = load_pipelines(Arc::clone(&app_config))?;
        let mut pipeline_map = HashMap::new();
        for pipeline in pipelines {
            let duplicate =
                pipeline_map.insert(pipeline.pipeline().id().clone(), Arc::new(pipeline));
            if let Some(duplicate_pipeline) = duplicate {
                log::warn!(
                    "The pipeline {:?} was overwritten due to pipeline ID {} not being unique.",
                    duplicate_pipeline,
                    duplicate_pipeline.pipeline().id()
                );
            }
        }
        Ok(pipeline_map)
    }
}

/// Returns all pipelines defined in the respective directory.
///
/// # Parameters
///
/// * `app_config` - the app [`Configuration`]
pub fn load_pipelines(
    app_config: Arc<Configuration>,
) -> Result<Vec<ContextualisedPipelineBlueprint>, SeqError> {
    let pipeline_path: PathBuf = app_config.pipeline_folder().into();
    let list = std::fs::read_dir(pipeline_path)?;
    let (dirs, errors): (Vec<std::fs::DirEntry>, Vec<std::io::Error>) =
        list.fold((Vec::new(), Vec::new()), |mut acc, entry| {
            if entry.is_ok() {
                acc.0.push(entry.unwrap());
            } else {
                acc.1.push(entry.unwrap_err());
            }
            acc
        });
    if errors.is_empty() {
        let pipeline_definition_files: Vec<PathBuf> = dirs
            .into_iter()
            .map(|dir| dir.path())
            .filter(|dir| dir.is_dir().into())
            .map(|mut dir| {
                dir.push(PIPELINE_DEFINITION_FILE);
                dir
            })
            .collect();
        let mut pipelines = Vec::with_capacity(pipeline_definition_files.len());
        for pipeline_def_file in pipeline_definition_files {
            let pipeline_def: PipelineBlueprint =
                serde_json::from_reader(std::fs::File::open(&pipeline_def_file)?)?;
            pipelines.push(ContextualisedPipelineBlueprint::new(
                pipeline_def,
                pipeline_def_file
                    .parent()
                    .expect("This unwrap of the parent path must work since we just appended the definition file."),
            ));
        }
        Ok(pipelines)
    } else {
        let combined_error = errors.into_iter().fold(String::new(), |mut acc, error| {
            acc.push_str(&error.to_string());
            acc.push('\n');
            acc
        });
        Err(SeqError::new(
            "std::io::Error",
            SeqErrorType::InternalServerError,
            combined_error,
            "An unforseen error occured during pipeline parsing. Please consult the logs.",
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use serial_test::serial;

    use crate::{
        application::config::Configuration,
        model::internal::pipeline_blueprint::PipelineStepVariableCategory,
        test_utility::{TestContext, TEST_RESOURCES_PATH},
    };

    use super::*;

    // The tests need to be serial to prevent the same context UUID to be issued
    // to different tests at the same time.

    #[test]
    #[serial]
    fn test_load_pipelines_folder_not_exist() {
        let context = TestContext::new();
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let app_config: Arc<Configuration> = Arc::new((&context).into());
        // Remove the folder that was automatically create with the test context.
        if PathBuf::from(context.pipeline_folder()).exists() {
            std::fs::remove_dir_all(context.pipeline_folder()).unwrap();
        }
        let pipelines = load_pipelines(app_config);
        // The pipeline folder does not exist.
        assert!(pipelines.is_err());
    }

    #[test]
    #[serial]
    fn test_load_pipelines_empty() {
        let context = TestContext::new();
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let app_config: Arc<Configuration> = Arc::new((&context).into());
        // Remove the folder that was automatically create with the test context.
        if PathBuf::from(context.pipeline_folder()).exists() {
            std::fs::remove_dir_all(context.pipeline_folder()).unwrap();
        }
        std::fs::create_dir_all(context.pipeline_folder()).unwrap();
        let pipelines = load_pipelines(app_config).unwrap();
        assert!(pipelines.is_empty());
    }

    #[test]
    #[serial]
    fn test_load_pipelines() {
        let context = TestContext::new();
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let app_config: Arc<Configuration> = Arc::new(Configuration::new(
            context.database_url(),
            "info",
            "127.0.0.1",
            "8080",
            context.context_folder(),
            format!("{}/pipelines", TEST_RESOURCES_PATH),
        ));
        let pipelines = load_pipelines(app_config).unwrap();
        assert_eq!(pipelines.len(), 1);
        let test_pipeline = &pipelines[0];
        let pipeline_context: PathBuf =
            format!("{}/pipelines/testing_pipeline", TEST_RESOURCES_PATH,).into();
        assert_eq!(test_pipeline.context(), &pipeline_context);
        let pipeline: &PipelineBlueprint = test_pipeline.pipeline();
        assert_eq!(pipeline.id(), "testing_pipeline");
        assert_eq!(pipeline.name(), "Testing pipeline");
        assert_eq!(pipeline.description(), "This pipeline is for testing purposes.");
        assert_eq!(pipeline.steps().len(), 1);
        let step = &pipeline.steps()[0];
        assert_eq!(step.id(), "fastqc");
        assert_eq!(step.name(), "FastQC");
        assert_eq!(step.description(), "Performs a quality control.");
        assert_eq!(step.container(), "fastqc");
        assert_eq!(step.dependencies(), &vec!["123", "456"]);
        assert_eq!(step.variables().len(), 5);
        assert_eq!(step.variables()[0].id(), "bool");
        assert_eq!(step.variables()[0].name(), "Boolean");
        assert_eq!(step.variables()[0].description(), "A boolean checkbox.");
        assert_eq!(step.variables()[0].category(), &PipelineStepVariableCategory::Boolean);
        assert_eq!(step.variables()[0].required(), &Some(true));
        assert_eq!(step.variables()[1].id(), "global");
        assert_eq!(step.variables()[1].name(), "Global");
        assert_eq!(step.variables()[1].description(), "A global data reference.");
        assert_eq!(step.variables()[1].category(), &PipelineStepVariableCategory::Global);
        assert_eq!(step.variables()[1].required(), &Some(false));
        assert_eq!(step.variables()[2].id(), "number");
        assert_eq!(step.variables()[2].name(), "Number");
        assert_eq!(step.variables()[2].description(), "A number field.");
        assert_eq!(step.variables()[2].category(), &PipelineStepVariableCategory::Number);
        assert_eq!(step.variables()[2].required(), &None);
        assert_eq!(step.variables()[3].id(), "option");
        assert_eq!(step.variables()[3].name(), "Option");
        assert_eq!(step.variables()[3].description(), "An option dropdown.");
        if let PipelineStepVariableCategory::Option(options) = step.variables()[3].category() {
            assert_eq!(options.len(), 2);
            assert_eq!(options[0].name(), "Option 1");
            assert_eq!(options[0].value(), "option1");
            assert_eq!(options[1].name(), "Option 2");
            assert_eq!(options[1].value(), "option2");
        } else {
            panic!("Not an option variable!");
        }
        assert_eq!(step.variables()[3].required(), &None);
        assert_eq!(step.variables()[4].id(), "string");
        assert_eq!(step.variables()[4].name(), "String");
        assert_eq!(step.variables()[4].description(), "A string text field.");
        assert_eq!(step.variables()[4].category(), &PipelineStepVariableCategory::String);
        assert_eq!(step.variables()[4].required(), &None);
    }

    #[test]
    #[serial]
    fn test_loaded_pipelines_has_varaible() {
        let context = TestContext::new();
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let app_config: web::Data<Configuration> = web::Data::new(Configuration::new(
            context.database_url(),
            "info",
            "127.0.0.1",
            "8080",
            context.context_folder(),
            format!("{}/pipelines", TEST_RESOURCES_PATH),
        ));
        let pipelines = LoadedPipelines::new(app_config).unwrap();
        assert!(pipelines.has_variable("testing_pipeline", "fastqc", "bool"));
        assert!(pipelines.has_variable("testing_pipeline", "fastqc", "global"));
        assert!(pipelines.has_variable("testing_pipeline", "fastqc", "number"));
        assert!(pipelines.has_variable("testing_pipeline", "fastqc", "option"));
        assert!(pipelines.has_variable("testing_pipeline", "fastqc", "string"));
        assert!(!pipelines.has_variable("invalid_pipeline", "fastqc", "string"));
        assert!(!pipelines.has_variable("testing_pipeline", "invalid_step", "string"));
        assert!(!pipelines.has_variable("testing_pipeline", "fastqc", "invalid_variable"));
    }
}
