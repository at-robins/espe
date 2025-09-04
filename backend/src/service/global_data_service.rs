use std::collections::HashSet;

use actix_web::web;
use diesel::SqliteConnection;

use crate::{
    application::error::{SeqError, SeqErrorType},
    model::{
        db::{
            experiment_execution::ExperimentExecution,
            pipeline_global_variable::PipelineGlobalVariable,
            pipeline_step_variable::PipelineStepVariable,
        },
        internal::pipeline_blueprint::{PipelineBlueprint, PipelineStepVariableCategory},
    },
    service::pipeline_service::LoadedPipelines,
};

#[derive(PartialEq, Eq, Hash)]
/// A helper struct to manage return values for lock checking.
/// Contains an experiment ID and an associated pipeline ID.
struct PipelineExperiment {
    pub pipeline_id: String,
    pub experiment_id: i32,
}

/// A helper struct to manage return values for lock checking.
/// Contains a global variable ID and an associated pipeline ID.
/// Represents a global variable that references a global data repository.
struct GlobalGlobalVariableID<'a> {
    pub variable_id: &'a String,
    pub pipeline_id: &'a String,
}

impl<'a> GlobalGlobalVariableID<'a> {
    /// Extracts all global variables that reference a global data repository
    /// from a pipeline blueprint.
    ///
    /// # Parameters
    ///
    /// * `blueprint` - the pipeline blueprint to extract the variables from
    pub fn generate(blueprint: &'a PipelineBlueprint) -> Vec<Self> {
        blueprint
            .global_variables()
            .iter()
            .filter(|global_var| global_var.category() == &PipelineStepVariableCategory::Global)
            .into_iter()
            .map(|global_data_var| GlobalGlobalVariableID {
                variable_id: global_data_var.id(),
                pipeline_id: blueprint.id(),
            })
            .collect()
    }

    /// Returns the experiment IDs and associated pipelines of all experiments that have a global pipeline variable set,
    /// which references the specified global data repository.
    ///
    /// # Parameters
    /// * `global_data_id` - the ID of the global data repository
    /// * `pipelines` - all loaded pipelines
    /// * `connection` - a connection to the database
    pub fn global_experiment_ids_with_pipeline(
        global_data_id: i32,
        pipelines: web::Data<LoadedPipelines>,
        connection: &mut SqliteConnection,
    ) -> Result<HashSet<PipelineExperiment>, SeqError> {
        let pipelines = pipelines.pipelines();
        let global_variable_global_ids: Vec<GlobalGlobalVariableID> = pipelines
            .iter()
            .flat_map(|c_pipeline| GlobalGlobalVariableID::generate(c_pipeline.pipeline()))
            .collect();
        let mut experiment_id_set = HashSet::new();
        for i in global_variable_global_ids {
            experiment_id_set.extend(
                PipelineGlobalVariable::get_by_pipeline_and_variable_id(
                    i.pipeline_id,
                    i.variable_id,
                    connection,
                )
                .map(|all_variables| {
                    all_variables
                        .into_iter()
                        .filter(|single_variable| {
                            single_variable.variable_value == Some(global_data_id.to_string())
                        })
                        .map(|single_variable| PipelineExperiment {
                            pipeline_id: single_variable.pipeline_id,
                            experiment_id: single_variable.experiment_id,
                        })
                })
                .map_err(|err| {
                    SeqError::from(err).chain(format!(
                        "Could not load global variables of \
                        the pipeline {} with ID {} from the database \
                        while generating an experiment ID set.",
                        i.pipeline_id, i.variable_id
                    ))
                })?,
            );
        }
        Ok(experiment_id_set)
    }
}

/// A helper struct to manage return values for lock checking.
/// Contains a step variable ID, an associated pipeline ID and a pipeline step ID.
/// Represents a step variable that references a global data repository.
struct GlobalStepVariableID<'a> {
    pub variable_id: &'a String,
    pub pipeline_id: &'a String,
    pub pipeline_step_id: &'a String,
}

impl<'a> GlobalStepVariableID<'a> {
    /// Extracts all step variables that reference a global data repository
    /// from a pipeline blueprint.
    ///
    /// # Parameters
    ///
    /// * `blueprint` - the pipeline blueprint to extract the variables from
    pub fn generate(blueprint: &'a PipelineBlueprint) -> Vec<Self> {
        blueprint
            .steps()
            .iter()
            .flat_map(|step| {
                step.variables()
                    .into_iter()
                    .map(|step_var| (step.id(), step_var))
            })
            .filter(|(_, global_var)| {
                global_var.category() == &PipelineStepVariableCategory::Global
            })
            .map(|(step_id, global_data_var)| GlobalStepVariableID {
                variable_id: global_data_var.id(),
                pipeline_id: blueprint.id(),
                pipeline_step_id: step_id,
            })
            .collect()
    }

    /// Returns the experiment IDs and associated pipelines of all experiments that have a pipeline step variable set,
    /// which reference the specified global data repository.
    ///
    /// # Parameters
    /// * `global_data_id` - the ID of the global data repository
    /// * `pipelines` - all loaded pipelines
    /// * `connection` - a connection to the database
    pub fn step_experiment_ids_with_pipeline(
        global_data_id: i32,
        pipelines: web::Data<LoadedPipelines>,
        connection: &mut SqliteConnection,
    ) -> Result<HashSet<PipelineExperiment>, SeqError> {
        let pipelines = pipelines.pipelines();
        let step_variable_global_ids: Vec<GlobalStepVariableID> = pipelines
            .iter()
            .flat_map(|c_pipeline| GlobalStepVariableID::generate(c_pipeline.pipeline()))
            .collect();
        let mut experiment_id_set = HashSet::new();
        for i in step_variable_global_ids {
            experiment_id_set.extend(
                PipelineStepVariable::get_by_pipeline_step_and_variable_id(
                    i.pipeline_id,
                    i.pipeline_step_id,
                    i.variable_id,
                    connection,
                )
                .map(|all_variables| {
                    all_variables
                        .into_iter()
                        .filter(|single_variable| {
                            single_variable.variable_value == Some(global_data_id.to_string())
                        })
                        .map(|single_variable| PipelineExperiment {
                            pipeline_id: single_variable.pipeline_id,
                            experiment_id: single_variable.experiment_id,
                        })
                })
                .map_err(|err| {
                    SeqError::from(err).chain(format!(
                        "Could not load step variables of \
                        the pipeline {} with ID {} from the database \
                        while generating an experiment ID set.",
                        i.pipeline_id, i.variable_id
                    ))
                })?,
            );
        }
        Ok(experiment_id_set)
    }
}

/// Returns `true` if the global data repo is currently locked.
///
/// # Parameters
///
/// * `global_data_id` - the ID of the global data repo
/// * `pipelines` - all loaded pipelines
/// * `connection` - a connection to the database
pub fn is_global_data_locked(
    global_data_id: i32,
    pipelines: web::Data<LoadedPipelines>,
    connection: &mut SqliteConnection,
) -> Result<bool, SeqError> {
    let global_experiments = GlobalGlobalVariableID::global_experiment_ids_with_pipeline(
        global_data_id,
        web::Data::clone(&pipelines),
        connection,
    )
    .map_err(|err| {
        err.chain(format!(
            "Loading the experiments and pipelines \
            containing global references to data \
            repository {} failed.",
            global_data_id
        ))
    })?;
    let step_experiment = GlobalStepVariableID::step_experiment_ids_with_pipeline(
        global_data_id,
        pipelines,
        connection,
    )
    .map_err(|err| {
        err.chain(format!(
            "Loading the experiments and pipelines \
            containing step references to data \
            repository {} failed.",
            global_data_id
        ))
    })?;
    for exp in global_experiments.union(&step_experiment) {
        if ExperimentExecution::is_executed_with_pipeline(
            exp.experiment_id,
            &exp.pipeline_id,
            connection,
        )
        .map_err(|err| {
            SeqError::from(err).chain(format!(
                "Obtaining lock state for global data repository {} failed \
                when checking if experiment {} is currently executed with pipeline {}.",
                global_data_id, exp.experiment_id, exp.pipeline_id
            ))
        })? {
            return Ok(true);
        }
    }
    return Ok(false);
}

/// Returns [`Ok`] if the global data repo is currently locked
/// and an error otherwise.
///
/// # Parameters
///
/// * `global_data_id` - the ID of the global data repo
/// * `pipelines` - all loaded pipelines
/// * `connection` - a connection to the database
pub fn is_global_data_locked_err(
    global_data_id: i32,
    pipelines: web::Data<LoadedPipelines>,
    connection: &mut SqliteConnection,
) -> Result<(), SeqError> {
    is_global_data_locked(global_data_id, pipelines, connection).and_then(|locked| {
        if locked {
            Err(SeqError::new(
                "Global data repository locked",
                SeqErrorType::PreconditionFailed,
                format!("Global data repository {} is locked.", global_data_id),
                "The global data repository is locked.",
            ))
        } else {
            Ok(())
        }
    })
}

#[cfg(test)]
mod tests {
    use actix_web::{http::StatusCode, ResponseError};
    use diesel::RunQueryDsl;

    use crate::{
        application::config::Configuration,
        model::db::{
            experiment_execution::ExecutionStatus, pipeline_step_variable::NewPipelineStepVariable,
        },
        test_utility::{
            create_default_experiment, create_default_global_data, TestContext,
            DEFAULT_EXPERIMENT_ID, DEFAULT_GLOBAL_DATA_ID, DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID, TEST_RESOURCES_PATH,
        },
    };

    use super::*;

    #[test]
    fn test_is_global_data_locked_aborted() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let mut context = TestContext::new();
        context.set_pipeline_folder(format!("{}/pipelines", TEST_RESOURCES_PATH));
        let mut connection = context.get_connection();
        let loaded_pipelines = web::Data::new(
            LoadedPipelines::new(web::Data::new(Configuration::from(&context))).unwrap(),
        );
        // Create an experiment containing all different stati.
        create_default_experiment(&mut connection);
        let new_records_all: Vec<ExperimentExecution> = vec![ExperimentExecution {
            id: 0,
            experiment_id: DEFAULT_EXPERIMENT_ID,
            pipeline_id: DEFAULT_PIPELINE_ID.to_string(),
            pipeline_step_id: DEFAULT_PIPELINE_STEP_ID.to_string(),
            execution_status: ExecutionStatus::Aborted.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        }];
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();

        create_default_global_data(&mut connection);
        // Associates the running pipeline with the test global data repository.
        let association_variable = NewPipelineStepVariable::new(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
            "global",
            Some(DEFAULT_GLOBAL_DATA_ID.to_string()),
        );
        diesel::insert_into(crate::schema::pipeline_step_variable::table)
            .values(&association_variable)
            .execute(&mut connection)
            .unwrap();
        assert!(!is_global_data_locked(
            DEFAULT_EXPERIMENT_ID,
            loaded_pipelines.clone(),
            &mut connection
        )
        .unwrap());
    }

    #[test]
    fn test_is_global_data_locked_failed() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let mut context = TestContext::new();
        context.set_pipeline_folder(format!("{}/pipelines", TEST_RESOURCES_PATH));
        let mut connection = context.get_connection();
        let loaded_pipelines = web::Data::new(
            LoadedPipelines::new(web::Data::new(Configuration::from(&context))).unwrap(),
        );
        // Create an experiment containing all different stati.
        create_default_experiment(&mut connection);
        let new_records_all: Vec<ExperimentExecution> = vec![ExperimentExecution {
            id: 0,
            experiment_id: DEFAULT_EXPERIMENT_ID,
            pipeline_id: DEFAULT_PIPELINE_ID.to_string(),
            pipeline_step_id: DEFAULT_PIPELINE_STEP_ID.to_string(),
            execution_status: ExecutionStatus::Failed.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        }];
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();

        create_default_global_data(&mut connection);
        // Associates the running pipeline with the test global data repository.
        let association_variable = NewPipelineStepVariable::new(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
            "global",
            Some(DEFAULT_GLOBAL_DATA_ID.to_string()),
        );
        diesel::insert_into(crate::schema::pipeline_step_variable::table)
            .values(&association_variable)
            .execute(&mut connection)
            .unwrap();
        assert!(!is_global_data_locked(
            DEFAULT_EXPERIMENT_ID,
            loaded_pipelines.clone(),
            &mut connection
        )
        .unwrap());
    }

    #[test]
    fn test_is_global_data_locked_finished() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let mut context = TestContext::new();
        context.set_pipeline_folder(format!("{}/pipelines", TEST_RESOURCES_PATH));
        let mut connection = context.get_connection();
        let loaded_pipelines = web::Data::new(
            LoadedPipelines::new(web::Data::new(Configuration::from(&context))).unwrap(),
        );
        // Create an experiment containing all different stati.
        create_default_experiment(&mut connection);
        let new_records_all: Vec<ExperimentExecution> = vec![ExperimentExecution {
            id: 0,
            experiment_id: DEFAULT_EXPERIMENT_ID,
            pipeline_id: DEFAULT_PIPELINE_ID.to_string(),
            pipeline_step_id: DEFAULT_PIPELINE_STEP_ID.to_string(),
            execution_status: ExecutionStatus::Finished.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        }];
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();

        create_default_global_data(&mut connection);
        // Associates the running pipeline with the test global data repository.
        let association_variable = NewPipelineStepVariable::new(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
            "global",
            Some(DEFAULT_GLOBAL_DATA_ID.to_string()),
        );
        diesel::insert_into(crate::schema::pipeline_step_variable::table)
            .values(&association_variable)
            .execute(&mut connection)
            .unwrap();
        assert!(!is_global_data_locked(
            DEFAULT_EXPERIMENT_ID,
            loaded_pipelines.clone(),
            &mut connection
        )
        .unwrap());
    }

    #[test]
    fn test_is_global_data_locked_running() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let mut context = TestContext::new();
        context.set_pipeline_folder(format!("{}/pipelines", TEST_RESOURCES_PATH));
        let mut connection = context.get_connection();
        let loaded_pipelines = web::Data::new(
            LoadedPipelines::new(web::Data::new(Configuration::from(&context))).unwrap(),
        );
        // Create an experiment containing all different stati.
        create_default_experiment(&mut connection);
        let new_records_all: Vec<ExperimentExecution> = vec![ExperimentExecution {
            id: 0,
            experiment_id: DEFAULT_EXPERIMENT_ID,
            pipeline_id: DEFAULT_PIPELINE_ID.to_string(),
            pipeline_step_id: DEFAULT_PIPELINE_STEP_ID.to_string(),
            execution_status: ExecutionStatus::Running.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        }];
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();

        create_default_global_data(&mut connection);
        // Associates the running pipeline with the test global data repository.
        let association_variable = NewPipelineStepVariable::new(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
            "global",
            Some(DEFAULT_GLOBAL_DATA_ID.to_string()),
        );
        diesel::insert_into(crate::schema::pipeline_step_variable::table)
            .values(&association_variable)
            .execute(&mut connection)
            .unwrap();
        assert!(is_global_data_locked(
            DEFAULT_EXPERIMENT_ID,
            loaded_pipelines.clone(),
            &mut connection
        )
        .unwrap());
    }

    #[test]
    fn test_is_global_data_locked_waiting() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let mut context = TestContext::new();
        context.set_pipeline_folder(format!("{}/pipelines", TEST_RESOURCES_PATH));
        let mut connection = context.get_connection();
        let loaded_pipelines = web::Data::new(
            LoadedPipelines::new(web::Data::new(Configuration::from(&context))).unwrap(),
        );
        // Create an experiment containing all different stati.
        create_default_experiment(&mut connection);
        let new_records_all: Vec<ExperimentExecution> = vec![ExperimentExecution {
            id: 0,
            experiment_id: DEFAULT_EXPERIMENT_ID,
            pipeline_id: DEFAULT_PIPELINE_ID.to_string(),
            pipeline_step_id: DEFAULT_PIPELINE_STEP_ID.to_string(),
            execution_status: ExecutionStatus::Waiting.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        }];
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();

        create_default_global_data(&mut connection);
        // Associates the running pipeline with the test global data repository.
        let association_variable = NewPipelineStepVariable::new(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
            "global",
            Some(DEFAULT_GLOBAL_DATA_ID.to_string()),
        );
        diesel::insert_into(crate::schema::pipeline_step_variable::table)
            .values(&association_variable)
            .execute(&mut connection)
            .unwrap();
        assert!(is_global_data_locked(
            DEFAULT_EXPERIMENT_ID,
            loaded_pipelines.clone(),
            &mut connection
        )
        .unwrap());
    }

    #[test]
    fn test_is_global_data_locked_err_aborted() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let mut context = TestContext::new();
        context.set_pipeline_folder(format!("{}/pipelines", TEST_RESOURCES_PATH));
        let mut connection = context.get_connection();
        let loaded_pipelines = web::Data::new(
            LoadedPipelines::new(web::Data::new(Configuration::from(&context))).unwrap(),
        );
        // Create an experiment containing all different stati.
        create_default_experiment(&mut connection);
        let new_records_all: Vec<ExperimentExecution> = vec![ExperimentExecution {
            id: 0,
            experiment_id: DEFAULT_EXPERIMENT_ID,
            pipeline_id: DEFAULT_PIPELINE_ID.to_string(),
            pipeline_step_id: DEFAULT_PIPELINE_STEP_ID.to_string(),
            execution_status: ExecutionStatus::Aborted.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        }];
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();

        create_default_global_data(&mut connection);
        // Associates the running pipeline with the test global data repository.
        let association_variable = NewPipelineStepVariable::new(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
            "global",
            Some(DEFAULT_GLOBAL_DATA_ID.to_string()),
        );
        diesel::insert_into(crate::schema::pipeline_step_variable::table)
            .values(&association_variable)
            .execute(&mut connection)
            .unwrap();
        assert!(is_global_data_locked_err(
            DEFAULT_EXPERIMENT_ID,
            loaded_pipelines.clone(),
            &mut connection
        )
        .is_ok());
    }

    #[test]
    fn test_is_global_data_locked_err_failed() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let mut context = TestContext::new();
        context.set_pipeline_folder(format!("{}/pipelines", TEST_RESOURCES_PATH));
        let mut connection = context.get_connection();
        let loaded_pipelines = web::Data::new(
            LoadedPipelines::new(web::Data::new(Configuration::from(&context))).unwrap(),
        );
        // Create an experiment containing all different stati.
        create_default_experiment(&mut connection);
        let new_records_all: Vec<ExperimentExecution> = vec![ExperimentExecution {
            id: 0,
            experiment_id: DEFAULT_EXPERIMENT_ID,
            pipeline_id: DEFAULT_PIPELINE_ID.to_string(),
            pipeline_step_id: DEFAULT_PIPELINE_STEP_ID.to_string(),
            execution_status: ExecutionStatus::Failed.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        }];
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();

        create_default_global_data(&mut connection);
        // Associates the running pipeline with the test global data repository.
        let association_variable = NewPipelineStepVariable::new(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
            "global",
            Some(DEFAULT_GLOBAL_DATA_ID.to_string()),
        );
        diesel::insert_into(crate::schema::pipeline_step_variable::table)
            .values(&association_variable)
            .execute(&mut connection)
            .unwrap();
        assert!(is_global_data_locked_err(
            DEFAULT_EXPERIMENT_ID,
            loaded_pipelines.clone(),
            &mut connection
        )
        .is_ok());
    }

    #[test]
    fn test_is_global_data_locked_err_finished() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let mut context = TestContext::new();
        context.set_pipeline_folder(format!("{}/pipelines", TEST_RESOURCES_PATH));
        let mut connection = context.get_connection();
        let loaded_pipelines = web::Data::new(
            LoadedPipelines::new(web::Data::new(Configuration::from(&context))).unwrap(),
        );
        // Create an experiment containing all different stati.
        create_default_experiment(&mut connection);
        let new_records_all: Vec<ExperimentExecution> = vec![ExperimentExecution {
            id: 0,
            experiment_id: DEFAULT_EXPERIMENT_ID,
            pipeline_id: DEFAULT_PIPELINE_ID.to_string(),
            pipeline_step_id: DEFAULT_PIPELINE_STEP_ID.to_string(),
            execution_status: ExecutionStatus::Finished.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        }];
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();

        create_default_global_data(&mut connection);
        // Associates the running pipeline with the test global data repository.
        let association_variable = NewPipelineStepVariable::new(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
            "global",
            Some(DEFAULT_GLOBAL_DATA_ID.to_string()),
        );
        diesel::insert_into(crate::schema::pipeline_step_variable::table)
            .values(&association_variable)
            .execute(&mut connection)
            .unwrap();
        assert!(is_global_data_locked_err(
            DEFAULT_EXPERIMENT_ID,
            loaded_pipelines.clone(),
            &mut connection
        )
        .is_ok());
    }

    #[test]
    fn test_is_global_data_locked_err_running() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let mut context = TestContext::new();
        context.set_pipeline_folder(format!("{}/pipelines", TEST_RESOURCES_PATH));
        let mut connection = context.get_connection();
        let loaded_pipelines = web::Data::new(
            LoadedPipelines::new(web::Data::new(Configuration::from(&context))).unwrap(),
        );
        // Create an experiment containing all different stati.
        create_default_experiment(&mut connection);
        let new_records_all: Vec<ExperimentExecution> = vec![ExperimentExecution {
            id: 0,
            experiment_id: DEFAULT_EXPERIMENT_ID,
            pipeline_id: DEFAULT_PIPELINE_ID.to_string(),
            pipeline_step_id: DEFAULT_PIPELINE_STEP_ID.to_string(),
            execution_status: ExecutionStatus::Running.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        }];
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();

        create_default_global_data(&mut connection);
        // Associates the running pipeline with the test global data repository.
        let association_variable = NewPipelineStepVariable::new(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
            "global",
            Some(DEFAULT_GLOBAL_DATA_ID.to_string()),
        );
        diesel::insert_into(crate::schema::pipeline_step_variable::table)
            .values(&association_variable)
            .execute(&mut connection)
            .unwrap();
        assert_eq!(
            is_global_data_locked_err(
                DEFAULT_EXPERIMENT_ID,
                loaded_pipelines.clone(),
                &mut connection
            )
            .unwrap_err()
            .status_code(),
            StatusCode::PRECONDITION_FAILED
        );
    }

    #[test]
    fn test_is_global_data_locked_err_waiting() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let mut context = TestContext::new();
        context.set_pipeline_folder(format!("{}/pipelines", TEST_RESOURCES_PATH));
        let mut connection = context.get_connection();
        let loaded_pipelines = web::Data::new(
            LoadedPipelines::new(web::Data::new(Configuration::from(&context))).unwrap(),
        );
        // Create an experiment containing all different stati.
        create_default_experiment(&mut connection);
        let new_records_all: Vec<ExperimentExecution> = vec![ExperimentExecution {
            id: 0,
            experiment_id: DEFAULT_EXPERIMENT_ID,
            pipeline_id: DEFAULT_PIPELINE_ID.to_string(),
            pipeline_step_id: DEFAULT_PIPELINE_STEP_ID.to_string(),
            execution_status: ExecutionStatus::Waiting.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        }];
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();

        create_default_global_data(&mut connection);
        // Associates the running pipeline with the test global data repository.
        let association_variable = NewPipelineStepVariable::new(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
            "global",
            Some(DEFAULT_GLOBAL_DATA_ID.to_string()),
        );
        diesel::insert_into(crate::schema::pipeline_step_variable::table)
            .values(&association_variable)
            .execute(&mut connection)
            .unwrap();
        assert_eq!(
            is_global_data_locked_err(
                DEFAULT_EXPERIMENT_ID,
                loaded_pipelines.clone(),
                &mut connection
            )
            .unwrap_err()
            .status_code(),
            StatusCode::PRECONDITION_FAILED
        );
    }
}
