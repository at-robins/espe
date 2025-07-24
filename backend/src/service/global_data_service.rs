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
            return Ok(false);
        }
    }
    return Ok(true);
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
    use std::path::Path;

    use actix_web::{http::StatusCode, ResponseError};
    use diesel::RunQueryDsl;

    use crate::{
        model::db::{experiment::NewExperiment, experiment_execution::ExecutionStatus},
        test_utility::TestContext,
    };

    use super::*;

    // #[test]
    // fn test_is_experiment_locked() {
    //     // Use a reference to the context, so the context is not dropped early
    //     // and messes up test context folder deletion.
    //     let context = TestContext::new();
    //     let mut connection = context.get_connection();
    //     // Create an experiment containing all different stati.
    //     let experiment_all = NewExperiment::new("all".to_string());
    //     let experiment_waiting = NewExperiment::new("waiting".to_string());
    //     let experiment_running = NewExperiment::new("running".to_string());
    //     let experiment_not_executed = NewExperiment::new("not executed".to_string());
    //     let experiment_empty = NewExperiment::new("empty".to_string());
    //     let experiment_id_all: i32 = diesel::insert_into(crate::schema::experiment::table)
    //         .values(&experiment_all)
    //         .returning(crate::schema::experiment::id)
    //         .get_result(&mut connection)
    //         .unwrap();
    //     let experiment_id_waiting: i32 = diesel::insert_into(crate::schema::experiment::table)
    //         .values(&experiment_waiting)
    //         .returning(crate::schema::experiment::id)
    //         .get_result(&mut connection)
    //         .unwrap();
    //     let experiment_id_running: i32 = diesel::insert_into(crate::schema::experiment::table)
    //         .values(&experiment_running)
    //         .returning(crate::schema::experiment::id)
    //         .get_result(&mut connection)
    //         .unwrap();
    //     let experiment_id_not_executed: i32 = diesel::insert_into(crate::schema::experiment::table)
    //         .values(&experiment_not_executed)
    //         .returning(crate::schema::experiment::id)
    //         .get_result(&mut connection)
    //         .unwrap();
    //     let experiment_id_empty: i32 = diesel::insert_into(crate::schema::experiment::table)
    //         .values(&experiment_empty)
    //         .returning(crate::schema::experiment::id)
    //         .get_result(&mut connection)
    //         .unwrap();
    //     let new_records_all: Vec<ExperimentExecution> = vec![
    //         (ExecutionStatus::Aborted, experiment_id_all),
    //         (ExecutionStatus::Failed, experiment_id_all),
    //         (ExecutionStatus::Finished, experiment_id_all),
    //         (ExecutionStatus::Running, experiment_id_all),
    //         (ExecutionStatus::Finished, experiment_id_all),
    //         (ExecutionStatus::Waiting, experiment_id_all),
    //         (ExecutionStatus::Finished, experiment_id_waiting),
    //         (ExecutionStatus::Waiting, experiment_id_waiting),
    //         (ExecutionStatus::Aborted, experiment_id_waiting),
    //         (ExecutionStatus::Finished, experiment_id_waiting),
    //         (ExecutionStatus::Finished, experiment_id_running),
    //         (ExecutionStatus::Running, experiment_id_running),
    //         (ExecutionStatus::Aborted, experiment_id_running),
    //         (ExecutionStatus::Finished, experiment_id_running),
    //         (ExecutionStatus::Finished, experiment_id_not_executed),
    //         (ExecutionStatus::Failed, experiment_id_not_executed),
    //         (ExecutionStatus::Aborted, experiment_id_not_executed),
    //         (ExecutionStatus::Finished, experiment_id_not_executed),
    //     ]
    //     .into_iter()
    //     .enumerate()
    //     .map(|(id, (status, experiment_id))| ExperimentExecution {
    //         id: id as i32,
    //         experiment_id,
    //         pipeline_id: id.to_string(),
    //         pipeline_step_id: id.to_string(),
    //         execution_status: status.into(),
    //         start_time: None,
    //         end_time: None,
    //         creation_time: chrono::Utc::now().naive_local(),
    //     })
    //     .collect();
    //     diesel::insert_into(crate::schema::experiment_execution::table)
    //         .values(&new_records_all)
    //         .execute(&mut connection)
    //         .unwrap();
    //     assert!(is_experiment_locked(experiment_id_all, &mut connection).unwrap());
    //     assert!(is_experiment_locked(experiment_id_waiting, &mut connection).unwrap());
    //     assert!(is_experiment_locked(experiment_id_running, &mut connection).unwrap());
    //     assert!(!is_experiment_locked(experiment_id_not_executed, &mut connection).unwrap());
    //     assert!(!is_experiment_locked(experiment_id_empty, &mut connection).unwrap());
    // }

    // #[test]
    // fn test_is_experiment_locked_err() {
    //     // Use a reference to the context, so the context is not dropped early
    //     // and messes up test context folder deletion.
    //     let context = TestContext::new();
    //     let mut connection = context.get_connection();
    //     // Create an experiment containing all different stati.
    //     let experiment_all = NewExperiment::new("all".to_string());
    //     let experiment_waiting = NewExperiment::new("waiting".to_string());
    //     let experiment_running = NewExperiment::new("running".to_string());
    //     let experiment_not_executed = NewExperiment::new("not executed".to_string());
    //     let experiment_empty = NewExperiment::new("empty".to_string());
    //     let experiment_id_all: i32 = diesel::insert_into(crate::schema::experiment::table)
    //         .values(&experiment_all)
    //         .returning(crate::schema::experiment::id)
    //         .get_result(&mut connection)
    //         .unwrap();
    //     let experiment_id_waiting: i32 = diesel::insert_into(crate::schema::experiment::table)
    //         .values(&experiment_waiting)
    //         .returning(crate::schema::experiment::id)
    //         .get_result(&mut connection)
    //         .unwrap();
    //     let experiment_id_running: i32 = diesel::insert_into(crate::schema::experiment::table)
    //         .values(&experiment_running)
    //         .returning(crate::schema::experiment::id)
    //         .get_result(&mut connection)
    //         .unwrap();
    //     let experiment_id_not_executed: i32 = diesel::insert_into(crate::schema::experiment::table)
    //         .values(&experiment_not_executed)
    //         .returning(crate::schema::experiment::id)
    //         .get_result(&mut connection)
    //         .unwrap();
    //     let experiment_id_empty: i32 = diesel::insert_into(crate::schema::experiment::table)
    //         .values(&experiment_empty)
    //         .returning(crate::schema::experiment::id)
    //         .get_result(&mut connection)
    //         .unwrap();
    //     let new_records_all: Vec<ExperimentExecution> = vec![
    //         (ExecutionStatus::Aborted, experiment_id_all),
    //         (ExecutionStatus::Failed, experiment_id_all),
    //         (ExecutionStatus::Finished, experiment_id_all),
    //         (ExecutionStatus::Running, experiment_id_all),
    //         (ExecutionStatus::Finished, experiment_id_all),
    //         (ExecutionStatus::Waiting, experiment_id_all),
    //         (ExecutionStatus::Finished, experiment_id_waiting),
    //         (ExecutionStatus::Waiting, experiment_id_waiting),
    //         (ExecutionStatus::Aborted, experiment_id_waiting),
    //         (ExecutionStatus::Finished, experiment_id_waiting),
    //         (ExecutionStatus::Finished, experiment_id_running),
    //         (ExecutionStatus::Running, experiment_id_running),
    //         (ExecutionStatus::Aborted, experiment_id_running),
    //         (ExecutionStatus::Finished, experiment_id_running),
    //         (ExecutionStatus::Finished, experiment_id_not_executed),
    //         (ExecutionStatus::Failed, experiment_id_not_executed),
    //         (ExecutionStatus::Aborted, experiment_id_not_executed),
    //         (ExecutionStatus::Finished, experiment_id_not_executed),
    //     ]
    //     .into_iter()
    //     .enumerate()
    //     .map(|(id, (status, experiment_id))| ExperimentExecution {
    //         id: id as i32,
    //         experiment_id,
    //         pipeline_id: id.to_string(),
    //         pipeline_step_id: id.to_string(),
    //         execution_status: status.into(),
    //         start_time: None,
    //         end_time: None,
    //         creation_time: chrono::Utc::now().naive_local(),
    //     })
    //     .collect();
    //     diesel::insert_into(crate::schema::experiment_execution::table)
    //         .values(&new_records_all)
    //         .execute(&mut connection)
    //         .unwrap();
    //     assert_eq!(
    //         is_experiment_locked_err(experiment_id_all, &mut connection)
    //             .unwrap_err()
    //             .status_code(),
    //         StatusCode::PRECONDITION_FAILED
    //     );
    //     assert_eq!(
    //         is_experiment_locked_err(experiment_id_waiting, &mut connection)
    //             .unwrap_err()
    //             .status_code(),
    //         StatusCode::PRECONDITION_FAILED
    //     );
    //     assert_eq!(
    //         is_experiment_locked_err(experiment_id_running, &mut connection)
    //             .unwrap_err()
    //             .status_code(),
    //         StatusCode::PRECONDITION_FAILED
    //     );
    //     assert!(is_experiment_locked_err(experiment_id_not_executed, &mut connection).is_ok());
    //     assert!(is_experiment_locked_err(experiment_id_empty, &mut connection).is_ok());
    // }
}
