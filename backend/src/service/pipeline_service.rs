use std::{
    collections::HashMap,
    path::{Path, PathBuf},
    sync::Arc,
};

use actix_web::web;
use parking_lot::Mutex;

use crate::{
    application::{
        config::{Configuration, PIPELINE_DEFINITION_FILE},
        error::{SeqError, SeqErrorType, DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE},
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
    pub fn new(app_config: web::Data<Configuration>) -> Result<Self, (Vec<SeqError>, Self)> {
        let (pipeline_map, errors) = match Self::load_pipeline_map(app_config) {
            Ok(pipeline_map) => (Mutex::new(pipeline_map), Vec::new()),
            Err((errors, pipeline_map)) => (
                Mutex::new(pipeline_map),
                errors
                    .into_iter()
                    .map(|err| {
                        err.chain(
                            "A pipeline could not be loaded during creation of all loaded pipeline data.",
                        )
                    })
                    .collect(),
            ),
        };
        if errors.is_empty() {
            Ok(Self { pipeline_map })
        } else {
            Err((errors, Self { pipeline_map }))
        }
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

    /// Sets the pipeline with the specified ID and returns the replaced value if any.
    ///
    /// # Parameters
    ///
    /// * `id` - the pipeline ID
    /// * `blueprint` - the new pipeline blueprint to set
    pub fn set<T: AsRef<str>>(
        &self,
        id: T,
        blueprint: ContextualisedPipelineBlueprint,
    ) -> Option<Arc<ContextualisedPipelineBlueprint>> {
        self.pipeline_map
            .lock()
            .insert(id.as_ref().to_string(), Arc::new(blueprint))
    }

    /// The number of pipelines loaded.
    pub fn size(&self) -> usize {
        self.pipeline_map.lock().len()
    }

    /// Returns ```true``` if the specified global variable exists in the loaded pipelines.
    ///
    /// # Parameters
    ///
    /// * `pipeline_id` - the ID of the pipeline the variable belongs to
    /// * `pipeline_step_id` - the ID of the pipeline step the variable belongs to
    /// * `variable_id` - the ID of the variable
    pub fn has_global_variable<T: AsRef<str>, S: AsRef<str>>(
        &self,
        pipeline_id: T,
        variable_id: S,
    ) -> bool {
        if let Some(pipeline) = self.get(pipeline_id) {
            pipeline
                .pipeline()
                .global_variables()
                .iter()
                .any(|variable| variable.id() == variable_id.as_ref())
        } else {
            false
        }
    }

    /// Returns ```true``` if the specified step variable exists in the loaded pipelines.
    ///
    /// # Parameters
    ///
    /// * `pipeline_id` - the ID of the pipeline the variable belongs to
    /// * `pipeline_step_id` - the ID of the pipeline step the variable belongs to
    /// * `variable_id` - the ID of the variable
    pub fn has_step_variable<T: AsRef<str>, R: AsRef<str>, S: AsRef<str>>(
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
    /// If errors occur during loading only successfully loaded pipelines
    /// are available and all remaining errors are logged.
    /// An error is only returned if no pipeline could be loaded.
    ///
    /// # Parameters
    ///
    /// * `app_config` - the app [`Configuration`]
    pub fn update_loaded_pipelines(
        &self,
        app_config: web::Data<Configuration>,
    ) -> Result<(), SeqError> {
        let (pipeline_map, errors) = match Self::load_pipeline_map(app_config) {
            Ok(pipeline_map) => (pipeline_map, Vec::new()),
            Err((errors, pipeline_map)) => (pipeline_map, errors.into_iter().map(|err| {
                err.chain("Updating some of the loaded pipelines failed. The pipeline map could not be loaded completely.")
            }).collect()),
        };

        // Log all errors.
        for error in &errors {
            error.log_default();
        }

        let is_pipeline_map_empty = pipeline_map.is_empty();

        if !is_pipeline_map_empty || errors.is_empty() {
            // Do not update the pipelines if there is clearly something wrong (empty pipeline map and errors),
            // but allow updating on partial errors.
            *(self.pipeline_map.lock()) = pipeline_map;
            log::warn!("Not all pipelines could be successfully updated.")
        }

        if is_pipeline_map_empty && !errors.is_empty() {
            // Only return an error if nothing could be loaded.
            Err(SeqError::new(
                "Pipeline update error",
                SeqErrorType::InternalServerError,
                "During updating all pipelines a severe error occured. Please check the individual pipeline loading logs for further information.",
                DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
            ))
        } else {
            Ok(())
        }
    }

    /// Updates a single currently loaded pipeline.
    ///
    /// # Parameters
    ///
    /// * `pipeline_id` - the ID of the pipeline
    pub fn update_loaded_pipeline<T: AsRef<str>>(&self, pipeline_id: T) -> Result<(), SeqError> {
        match self.get(&pipeline_id) {
            Some(loaded_pipeline) => {
                self.set(
                    pipeline_id.as_ref(),
                    load_pipeline(loaded_pipeline.context()).map_err(|err| {
                        err.chain(format!(
                            "The pipeline with ID {} could not be updated.",
                            pipeline_id.as_ref()
                        ))
                    })?,
                );
                Ok(())
            },
            None => Err(SeqError::new(
                "Pipeline not found",
                SeqErrorType::NotFoundError,
                format!(
                    "The pipeline with ID {} is not loaded and can thus not be updated.",
                    pipeline_id.as_ref()
                ),
                "The pipeline is not loaded. For further information please check the logs.",
            )),
        }
    }

    /// Creates a map of loaded pipelines by their respective ID.
    /// If loading of specific pipelines fails the map of pipelines
    /// that could be loaded is returned together with all errors
    /// that occurred during loading.
    ///
    /// # Parameters
    ///
    /// * `app_config` - the app [`Configuration`]
    fn load_pipeline_map(
        app_config: web::Data<Configuration>,
    ) -> Result<
        HashMap<String, Arc<ContextualisedPipelineBlueprint>>,
        (Vec<SeqError>, HashMap<String, Arc<ContextualisedPipelineBlueprint>>),
    > {
        let (pipelines, errors) = match load_pipelines(Arc::clone(&app_config)) {
            Ok(pipelines) => (pipelines, Vec::new()),
            Err((errors, pipelines)) => (
                pipelines,
                errors
                    .into_iter()
                    .map(|err| {
                        err.chain("A pipeline could not be loaded during pipeline map generation.")
                    })
                    .collect(),
            ),
        };
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
        if errors.is_empty() {
            Ok(pipeline_map)
        } else {
            Err((errors, pipeline_map))
        }
    }
}

/// Loads and returns a ['ContextualisedPipelineBlueprint'] from the specified file.
///
/// # Parameters
///
/// * `pipeline_definition_path` - the path to the pipeline definition file
pub fn load_pipeline_definition<P: AsRef<Path>>(
    pipeline_definition_path: P,
) -> Result<ContextualisedPipelineBlueprint, SeqError> {
    let pipeline_def: PipelineBlueprint =
        serde_json::from_reader(std::fs::File::open(&pipeline_definition_path)?)?;
    let mut contextualised_pipeline = ContextualisedPipelineBlueprint::new(
        pipeline_def,
        pipeline_definition_path.as_ref().parent().expect(
            "This unwrap of the parent path must work since we just deserialised the definition file.",
        ),
    );
    // Loads potential data that is not directly present in the pipeline definition.
    contextualised_pipeline.resolve_imports()?;

    Ok(contextualised_pipeline)
}

/// Loads and returns a ['ContextualisedPipelineBlueprint'] from the specified directory.
/// This will return an error if the specified path is not a directory.
///
/// # Parameters
///
/// * `pipeline_path` - the path to the pipeline directory
pub fn load_pipeline<P: AsRef<Path>>(
    pipeline_path: P,
) -> Result<ContextualisedPipelineBlueprint, SeqError> {
    if pipeline_path.as_ref().is_dir() {
        load_pipeline_definition(pipeline_path.as_ref().join(PIPELINE_DEFINITION_FILE)).map_err(
            |err| {
                err.chain(format!(
                    "The pipeline directory {} could not be loaded.",
                    pipeline_path.as_ref().display()
                ))
            },
        )
    } else {
        Err(SeqError::new(
            "Pipeline loading error",
            SeqErrorType::InternalServerError,
            format!("The pipeline path {} is not a directory.", pipeline_path.as_ref().display()),
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        ))
    }
}

/// Returns all pipelines defined in the respective directory.
/// In case of an error it returns all errors as well as all pipelines
/// that could be read successfully.
///
/// # Parameters
///
/// * `app_config` - the app [`Configuration`]
pub fn load_pipelines(
    app_config: Arc<Configuration>,
) -> Result<
    Vec<ContextualisedPipelineBlueprint>,
    (Vec<SeqError>, Vec<ContextualisedPipelineBlueprint>),
> {
    let pipeline_path: PathBuf = app_config.pipeline_folder().into();
    let list = match std::fs::read_dir(&pipeline_path) {
        Ok(list) => list,
        Err(err) => {
            return Err((
                vec![SeqError::from(err).chain(format!(
                    "Could not read pipeline directory {}.",
                    pipeline_path.display()
                ))],
                Vec::new(),
            ))
        },
    };

    let (dirs, mut errors): (Vec<std::fs::DirEntry>, Vec<SeqError>) =
        list.fold((Vec::new(), Vec::new()), |mut acc, entry| {
            if entry.is_ok() {
                acc.0.push(entry.unwrap());
            } else {
                acc.1.push(SeqError::from(entry.unwrap_err()).chain(format!(
                    "Error while reading the sub-directories of the pipeline directory {}.",
                    pipeline_path.display(),
                )));
            }
            acc
        });
    // Only loads folders that actually contain a pipeline definition file.
    let pipeline_definition_files: Vec<PathBuf> = dirs
        .into_iter()
        .map(|dir| dir.path())
        .filter(|dir| dir.is_dir())
        .map(|mut dir| {
            dir.push(PIPELINE_DEFINITION_FILE);
            dir
        })
        .filter(|file| file.is_file())
        .collect();
    let mut pipelines = Vec::with_capacity(pipeline_definition_files.len());
    for pipeline_def_file in pipeline_definition_files {
        match load_pipeline_definition(&pipeline_def_file) {
            Ok(pipeline) => pipelines.push(pipeline),
            Err(err) => errors.push(err.chain(format!(
                "Pipeline definition file {} could not be loaded.",
                pipeline_def_file.display()
            ))),
        }
    }
    if errors.is_empty() {
        Ok(pipelines)
    } else {
        Err((errors, pipelines))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use serial_test::serial;

    use crate::{
        application::config::{ApplicationMode, Configuration},
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
            ApplicationMode::Release,
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
    fn test_loaded_pipelines_has_step_varaible() {
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
            ApplicationMode::Release,
        ));
        let pipelines = LoadedPipelines::new(app_config).unwrap();
        assert!(pipelines.has_step_variable("testing_pipeline", "fastqc", "bool"));
        assert!(pipelines.has_step_variable("testing_pipeline", "fastqc", "global"));
        assert!(pipelines.has_step_variable("testing_pipeline", "fastqc", "number"));
        assert!(pipelines.has_step_variable("testing_pipeline", "fastqc", "option"));
        assert!(pipelines.has_step_variable("testing_pipeline", "fastqc", "string"));
        assert!(!pipelines.has_step_variable("invalid_pipeline", "fastqc", "string"));
        assert!(!pipelines.has_step_variable("testing_pipeline", "invalid_step", "string"));
        assert!(!pipelines.has_step_variable("testing_pipeline", "fastqc", "invalid_variable"));
    }

    #[test]
    #[serial]
    fn test_loaded_pipelines_has_global_varaible() {
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
            ApplicationMode::Release,
        ));
        let pipelines = LoadedPipelines::new(app_config).unwrap();
        assert!(pipelines.has_global_variable("testing_pipeline", "global_number"));
        assert!(!pipelines.has_global_variable("invalid_pipeline", "string"));
        assert!(!pipelines.has_global_variable("testing_pipeline", "invalid_variable"));
    }
}
