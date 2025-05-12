use std::collections::HashMap;

use crate::{
    application::{
        config::{Configuration, LogProcessType},
        database::DatabaseManager,
        error::{SeqError, SeqErrorType},
    },
    diesel::RunQueryDsl,
    model::{
        db::{
            experiment::{Experiment, NewExperiment},
            experiment_execution::{ExecutionStatus, ExperimentExecution, NewExperimentExecution},
            pipeline_global_variable::{NewPipelineGlobalVariable, PipelineGlobalVariable},
            pipeline_step_variable::{NewPipelineStepVariable, PipelineStepVariable},
        },
        exchange::{
            experiment_details::ExperimentDetails,
            experiment_pipeline::ExperimentPipelineBlueprint,
            pipeline_variable_upload::{PipelineGlobalVariableUpload, PipelineStepVariableUpload},
        },
    },
    service::{
        execution_service::ExecutionScheduler,
        experiment_service::{delete_step_logs, delete_step_output},
        pipeline_service::LoadedPipelines,
        validation_service::{validate_comment, validate_entity_name, validate_mail},
    },
};
use actix_web::{web, HttpResponse};
use chrono::NaiveDateTime;
use diesel::{ExpressionMethods, QueryDsl};
use parking_lot::Mutex;

pub async fn create_experiment(
    database_manager: web::Data<DatabaseManager>,
    name: actix_web::web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let name: String = name.into_inner();
    validate_entity_name(&name)?;
    let mut connection = database_manager.database_connection()?;
    log::info!("Creating experiment with name {}.", &name);
    let new_record = NewExperiment::new(name);
    let inserted_id: i32 = diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .returning(crate::schema::experiment::id)
        .get_result(&mut connection)?;
    // Return the ID of the created experiment.
    Ok(HttpResponse::Created().json(inserted_id))
}

pub async fn delete_experiment(
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Retrieve the app config.
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    log::info!("Deleting experiment with ID {}.", id);
    // Remove all files belonging to the experiment.
    let experiment_path = app_config.experiment_path(id.to_string());
    if experiment_path.exists() {
        std::fs::remove_dir_all(experiment_path)?;
    }
    // Delete the experiment from the database.
    connection.immediate_transaction(|connection| {
        diesel::delete(crate::schema::experiment::table)
            .filter(crate::schema::experiment::id.eq(id))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn get_experiment(
    database_manager: web::Data<DatabaseManager>,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    let experiment_details: ExperimentDetails = crate::schema::experiment::table
        .find(id)
        .first::<Experiment>(&mut connection)?
        .into();
    Ok(HttpResponse::Ok().json(experiment_details))
}

pub async fn get_experiment_execution_status(
    database_manager: web::Data<DatabaseManager>,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    let execution_steps = ExperimentExecution::get_by_experiment(id, &mut connection)?;
    let result = if execution_steps.is_empty() {
        "None".to_string()
    } else if execution_steps
        .iter()
        .any(|execution| execution.execution_status == ExecutionStatus::Running.to_string())
    {
        ExecutionStatus::Running.to_string()
    } else if execution_steps
        .iter()
        .any(|execution| execution.execution_status == ExecutionStatus::Failed.to_string())
    {
        ExecutionStatus::Failed.to_string()
    } else if execution_steps
        .iter()
        .any(|execution| execution.execution_status == ExecutionStatus::Aborted.to_string())
    {
        ExecutionStatus::Aborted.to_string()
    } else if execution_steps
        .iter()
        .all(|execution| execution.execution_status == ExecutionStatus::Finished.to_string())
    {
        ExecutionStatus::Finished.to_string()
    } else {
        ExecutionStatus::Waiting.to_string()
    };
    Ok(HttpResponse::Ok().json(result))
}

pub async fn patch_experiment_name(
    database_manager: web::Data<DatabaseManager>,
    id: web::Path<i32>,
    new_name: web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_name = new_name.into_inner();
    validate_entity_name(&new_name)?;
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::experiment_name.eq(new_name))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn patch_experiment_mail(
    database_manager: web::Data<DatabaseManager>,
    id: web::Path<i32>,
    new_name: web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_mail = new_name.into_inner();
    validate_mail(&new_mail)?;
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::mail.eq(new_mail))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn patch_experiment_comment(
    database_manager: web::Data<DatabaseManager>,
    id: web::Path<i32>,
    new_comment: web::Json<Option<String>>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Sanitise the HTML and validate.
    let new_comment = new_comment.into_inner().map(|inner| ammonia::clean(&inner));
    if let Some(inner) = &new_comment {
        validate_comment(inner)?;
    }
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::comment.eq(new_comment))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn patch_experiment_pipeline(
    database_manager: web::Data<DatabaseManager>,
    pipelines: web::Data<LoadedPipelines>,
    id: web::Path<i32>,
    new_pipeline: web::Json<Option<String>>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_pipeline = new_pipeline.into_inner();
    if let Some(new_pipeline_id) = &new_pipeline {
        if !pipelines.is_loaded(new_pipeline_id) {
            return Err(SeqError::new(
                "Not Found",
                SeqErrorType::NotFoundError,
                format!("No pipeline with ID {} is currently loaded.", new_pipeline_id),
                "The pipeline ID is invalid.",
            ));
        }
    }
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::pipeline_id.eq(new_pipeline))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn list_experiment(
    database_manager: web::Data<DatabaseManager>,
) -> Result<web::Json<Vec<ExperimentDetails>>, SeqError> {
    let mut connection = database_manager.database_connection()?;
    let experiments: Vec<ExperimentDetails> = Experiment::get_all(&mut connection)?
        .into_iter()
        .map(|val| val.into())
        .collect();
    Ok(web::Json(experiments))
}

pub async fn get_experiment_pipelines(
    database_manager: web::Data<DatabaseManager>,
    pipelines: web::Data<LoadedPipelines>,
    id: web::Path<i32>,
) -> Result<web::Json<Vec<ExperimentPipelineBlueprint>>, SeqError> {
    let experiment_id: i32 = id.into_inner();
    let mut connection = database_manager.database_connection()?;
    let all_experiment_stati =
        ExperimentExecution::get_by_experiment(experiment_id, &mut connection)?;
    let mut experiment_pipelines = Vec::new();
    for pipeline in pipelines.pipelines() {
        let values_global = crate::model::db::pipeline_global_variable::PipelineGlobalVariable::get_values_by_experiment_and_pipeline(experiment_id, pipeline.pipeline().id(), &mut connection)?;
        let values_step = crate::model::db::pipeline_step_variable::PipelineStepVariable::get_values_by_experiment_and_pipeline(experiment_id, pipeline.pipeline().id(), &mut connection)?;
        let stati: HashMap<String, String> = all_experiment_stati
            .iter()
            .filter(|execution| &execution.pipeline_id == pipeline.pipeline().id())
            .map(|execution| {
                (execution.pipeline_step_id.clone(), execution.execution_status.clone())
            })
            .collect();
        experiment_pipelines.push(ExperimentPipelineBlueprint::from_internal(
            pipeline.pipeline(),
            values_global,
            values_step,
            stati,
        ));
    }
    Ok(web::Json(experiment_pipelines))
}

pub async fn get_experiment_pipeline_run(
    database_manager: web::Data<DatabaseManager>,
    pipelines: web::Data<LoadedPipelines>,
    id: web::Path<i32>,
) -> Result<web::Json<Option<ExperimentPipelineBlueprint>>, SeqError> {
    let experiment_id: i32 = id.into_inner();
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(experiment_id, &mut connection)?;
    let experiment = Experiment::get(experiment_id, &mut connection)?;
    let experiment_pipeline = if let Some(pipeline_id) = experiment.pipeline_id {
        if let Some(pipeline) = pipelines.get(&pipeline_id) {
            let values_global = crate::model::db::pipeline_global_variable::PipelineGlobalVariable::get_values_by_experiment_and_pipeline(experiment_id, pipeline.pipeline().id(), &mut connection)?;
            let values_step = crate::model::db::pipeline_step_variable::PipelineStepVariable::get_values_by_experiment_and_pipeline(experiment_id, pipeline.pipeline().id(), &mut connection)?;
            let stati: HashMap<String, String> =
                ExperimentExecution::get_by_experiment(experiment_id, &mut connection)?
                    .into_iter()
                    .filter(|execution| &execution.pipeline_id == &pipeline_id)
                    .map(|execution| (execution.pipeline_step_id, execution.execution_status))
                    .collect();
            Some(ExperimentPipelineBlueprint::from_internal(
                pipeline.pipeline(),
                values_global,
                values_step,
                stati,
            ))
        } else {
            return Err(SeqError::new(
                "Not Found",
                SeqErrorType::NotFoundError,
                format!("No pipeline with ID {} is currently loaded.", pipeline_id),
                "The pipeline ID is invalid.",
            ));
        }
    } else {
        None
    };
    Ok(web::Json(experiment_pipeline))
}

pub async fn post_experiment_pipeline_step_variable(
    database_manager: web::Data<DatabaseManager>,
    pipelines: web::Data<LoadedPipelines>,
    experiment_id: web::Path<i32>,
    new_variable: web::Json<PipelineStepVariableUpload>,
) -> Result<HttpResponse, SeqError> {
    let experiment_id: i32 = experiment_id.into_inner();
    let new_variable: PipelineStepVariableUpload = new_variable.into_inner();
    if !pipelines.has_step_variable(
        &new_variable.pipeline_id,
        &new_variable.pipeline_step_id,
        &new_variable.variable_id,
    ) {
        return Err(SeqError::new(
            "Not Found",
            SeqErrorType::NotFoundError,
            format!(
                "No pipeline step variable with corresponding properties is currently loaded, thus variable {:?} cannot be inserted.",
                new_variable
            ),
            "The variable is invalid.",
        ));
    }
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(experiment_id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        if let Some(existing_variable) = PipelineStepVariable::get(
            experiment_id,
            &new_variable.pipeline_id,
            &new_variable.pipeline_step_id,
            &new_variable.variable_id,
            connection,
        )? {
            // Update if the variable already exists.
            diesel::update(
                crate::schema::pipeline_step_variable::table
                    .filter(crate::schema::pipeline_step_variable::id.eq(existing_variable.id)),
            )
            .set(
                crate::schema::pipeline_step_variable::variable_value
                    .eq(new_variable.variable_value),
            )
            .execute(connection)
        } else {
            // Insert if the variable does not exist.
            let new_variable = NewPipelineStepVariable::new(
                experiment_id,
                new_variable.pipeline_id,
                new_variable.pipeline_step_id,
                new_variable.variable_id,
                new_variable.variable_value,
            );
            diesel::insert_into(crate::schema::pipeline_step_variable::table)
                .values(&new_variable)
                .execute(connection)
        }
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn post_experiment_pipeline_global_variable(
    database_manager: web::Data<DatabaseManager>,
    pipelines: web::Data<LoadedPipelines>,
    experiment_id: web::Path<i32>,
    new_variable: web::Json<PipelineGlobalVariableUpload>,
) -> Result<HttpResponse, SeqError> {
    let experiment_id: i32 = experiment_id.into_inner();
    let new_variable: PipelineGlobalVariableUpload = new_variable.into_inner();
    if !pipelines.has_global_variable(&new_variable.pipeline_id, &new_variable.variable_id) {
        return Err(SeqError::new(
            "Not Found",
            SeqErrorType::NotFoundError,
            format!(
                "No global pipeline variable with corresponding properties is currently loaded, thus variable {:?} cannot be inserted.",
                new_variable
            ),
            "The variable is invalid.",
        ));
    }
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(experiment_id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        if let Some(existing_variable) = PipelineGlobalVariable::get(
            experiment_id,
            &new_variable.pipeline_id,
            &new_variable.variable_id,
            connection,
        )? {
            // Update if the variable already exists.
            diesel::update(
                crate::schema::pipeline_global_variable::table
                    .filter(crate::schema::pipeline_global_variable::id.eq(existing_variable.id)),
            )
            .set(
                crate::schema::pipeline_global_variable::variable_value
                    .eq(new_variable.variable_value),
            )
            .execute(connection)
        } else {
            // Insert if the variable does not exist.
            let new_variable = NewPipelineGlobalVariable::new(
                experiment_id,
                new_variable.pipeline_id,
                new_variable.variable_id,
                new_variable.variable_value,
            );
            diesel::insert_into(crate::schema::pipeline_global_variable::table)
                .values(&new_variable)
                .execute(connection)
        }
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn post_execute_experiment(
    database_manager: web::Data<DatabaseManager>,
    experiment_id: web::Path<i32>,
    pipelines: web::Data<LoadedPipelines>,
) -> Result<HttpResponse, SeqError> {
    let experiment_id: i32 = experiment_id.into_inner();
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(experiment_id, &mut connection)?;
    // Error if the experiment was already submitted for execution.
    if ExperimentExecution::has_experiment_execution_entries(experiment_id, &mut connection)? {
        return Err(SeqError::new(
            "Invalid run",
            SeqErrorType::BadRequestError,
            format!("The experiment {} was already executed.", experiment_id),
            "The requested run is invalid.",
        ));
    }
    if let Some(pipeline_id) = Experiment::get(experiment_id, &mut connection)?.pipeline_id {
        if let Some(pipeline) = pipelines.get(&pipeline_id) {
            let pipeline = pipeline.pipeline();
            // Checks if the required variables are set.
            PipelineGlobalVariable::validate_global_variables(
                pipeline,
                experiment_id,
                &mut connection,
            )?;
            let mut execution_steps: Vec<NewExperimentExecution> =
                Vec::with_capacity(pipeline.steps().len());
            for step in pipeline.steps() {
                PipelineStepVariable::validate_step_variables(
                    step,
                    experiment_id,
                    &pipeline_id,
                    &mut connection,
                )?;
                execution_steps.push(NewExperimentExecution::new(
                    experiment_id,
                    pipeline.id(),
                    step.id(),
                ));
            }
            connection.immediate_transaction(|connection| {
                diesel::insert_into(crate::schema::experiment_execution::table)
                    .values(&execution_steps)
                    .execute(connection)
            })?;
            log::info!(
                "Submitted experiment {} with pipeline {} for execution.",
                experiment_id,
                pipeline.id()
            );
            Ok(HttpResponse::Ok().finish())
        } else {
            // Error if the pipeline is not loaded.
            Err(SeqError::new(
                "Invalid run",
                SeqErrorType::BadRequestError,
                format!(
                    "The selected pipeline {} for experiment {} is not loaded.",
                    pipeline_id, experiment_id
                ),
                "The requested run parameters are invalid.",
            ))
        }
    } else {
        // Error if no pipeline was selected.
        Err(SeqError::new(
            "Invalid run",
            SeqErrorType::BadRequestError,
            format!(
                "No pipeline was selected for experiment {}, so it cannot be run.",
                experiment_id
            ),
            "The requested run parameters are invalid.",
        ))
    }
}

pub async fn post_execute_experiment_step(
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    experiment_id: web::Path<i32>,
    step_id: web::Json<String>,
    pipelines: web::Data<LoadedPipelines>,
) -> Result<HttpResponse, SeqError> {
    let experiment_id: i32 = experiment_id.into_inner();
    let step_id: String = step_id.into_inner();
    let mut connection = database_manager.database_connection().map_err(|err| {
        err.chain(format!(
            "Connection to the database could not be aquired when executing pipeline step {} of experiment {}.",
            step_id, experiment_id
        ))
    })?;
    // Error if the experiment does not exist.
    Experiment::exists_err(experiment_id, &mut connection).map_err(|err| {
        err.chain(format!(
            "Experiment {} does not exist. Pipeline step {} cannot be executed.",
            experiment_id, step_id
        ))
    })?;
    if let Some(pipeline_id) = Experiment::get(experiment_id, &mut connection)?.pipeline_id {
        if let Some(pipeline) = pipelines.get(&pipeline_id) {
            let pipeline = pipeline.pipeline();
            // Checks if the step exists within the currently selected pipeline
            if let Some(step) = pipeline.steps().iter().find(|step| step.id() == &step_id) {
                // Checks if the dependencies are satisfied.
                let executions: Vec<ExperimentExecution> =
                    ExperimentExecution::get_by_experiment(experiment_id, &mut connection)?
                        .into_iter()
                        .filter(|s| &s.pipeline_id == &pipeline_id)
                        .collect();
                let satisfied_dependencies: Vec<&String> = executions
                    .iter()
                    .filter(|s| {
                        s.execution_status == ExecutionStatus::Finished.to_string()
                            || s.execution_status == ExecutionStatus::Running.to_string()
                            || s.execution_status == ExecutionStatus::Waiting.to_string()
                    })
                    .map(|s| &s.pipeline_step_id)
                    .collect();
                let are_dependencies_satisfied: bool = step
                    .dependencies()
                    .iter()
                    .all(|dependency| satisfied_dependencies.contains(&dependency));
                if !are_dependencies_satisfied {
                    return Err(SeqError::new(
                        "Invalid run",
                        SeqErrorType::BadRequestError,
                        format!("The experiment {} is missing dependencies for execution of pipeline {} step {}.\nRequired dependencies: {:?}\nSatisfied dependencies: {:?}", experiment_id, pipeline.id(), step.id(), step.dependencies(), satisfied_dependencies),
                        "The requested run parameters are invalid.",
                    ));
                }
                // Checks if the required variables are set.
                PipelineGlobalVariable::validate_global_variables(
                    pipeline,
                    experiment_id,
                    &mut connection,
                )?;
                PipelineStepVariable::validate_step_variables(
                    step,
                    experiment_id,
                    &pipeline_id,
                    &mut connection,
                )?;
                // Submits execution step.
                if let Some(existing_execution) =
                    executions.iter().find(|s| s.pipeline_step_id == step_id)
                {
                    // Restart an existing pipeline step.
                    if existing_execution.execution_status == ExecutionStatus::Running.to_string()
                        || existing_execution.execution_status
                            == ExecutionStatus::Waiting.to_string()
                    {
                        // Error if the step is currently scheduled for execution.
                        return Err(SeqError::new(
                            "Invalid run",
                            SeqErrorType::BadRequestError,
                            format!("The experiment {} pipeline {} step {} is already scheduled for execution and can thus not be restarted.", experiment_id, pipeline.id(), step.id()),
                            "The requested run parameters are invalid.",
                        ));
                    }
                    // Deletes the output folder and run logs, but keeps the build logs
                    // as the build process might be cached / skipped.
                    delete_step_output(&step_id, experiment_id, web::Data::clone(&app_config))
                        .map_err(|err| {
                            err.chain(format!(
                                "Deletion of pipeline step output ({}/{}) of experiment {} failed during restart request.",
                                pipeline_id, step_id, experiment_id
                            ))
                        })?;
                    delete_step_logs(
                        &pipeline_id,
                        &step_id,
                        experiment_id,
                        &vec![LogProcessType::Run],
                        app_config,
                    ).map_err(|err| {
                        err.chain(format!(
                            "Deletion of pipeline step run logs ({}/{}) of experiment {} failed during restart request.",
                            pipeline_id, step_id, experiment_id
                        ))
                    })?;
                    connection.immediate_transaction(|connection| {
                        let clear_time: Option<NaiveDateTime> = None;
                        diesel::update(
                            crate::schema::experiment_execution::table.find(existing_execution.id),
                        )
                        .set((
                            crate::schema::experiment_execution::execution_status
                                .eq(ExecutionStatus::Waiting.to_string()),
                            crate::schema::experiment_execution::start_time.eq(clear_time.clone()),
                            crate::schema::experiment_execution::end_time.eq(clear_time),
                        ))
                        .execute(connection)
                    })?;
                } else {
                    // Create a newly added pipeline step.
                    let execution_step: NewExperimentExecution =
                        NewExperimentExecution::new(experiment_id, pipeline.id(), step.id());

                    connection.immediate_transaction(|connection| {
                        diesel::insert_into(crate::schema::experiment_execution::table)
                            .values(&execution_step)
                            .execute(connection)
                    })?;
                }
                log::info!(
                    "Submitted experiment {} with pipeline {} step {} for execution.",
                    experiment_id,
                    pipeline.id(),
                    step_id
                );
                Ok(HttpResponse::Ok().finish())
            } else {
                // Error the pipeline does not contain the step.
                Err(SeqError::new(
                    "Invalid run",
                    SeqErrorType::BadRequestError,
                    format!(
                        "The selected pipeline {} for experiment {} is does not contain step {}.",
                        pipeline_id, experiment_id, step_id
                    ),
                    "The requested run parameters are invalid.",
                ))
            }
        } else {
            // Error if the pipeline is not loaded.
            Err(SeqError::new(
                "Invalid run",
                SeqErrorType::BadRequestError,
                format!(
                    "The selected pipeline {} for experiment {} is not loaded.",
                    pipeline_id, experiment_id
                ),
                "The requested run parameters are invalid.",
            ))
        }
    } else {
        // Error if no pipeline was selected.
        Err(SeqError::new(
            "Invalid run",
            SeqErrorType::BadRequestError,
            format!(
                "No pipeline was selected for experiment {}, so it cannot be run.",
                experiment_id
            ),
            "The requested run parameters are invalid.",
        ))
    }
}

pub async fn post_experiment_execution_abort(
    database_manager: web::Data<DatabaseManager>,
    scheduler: web::Data<Mutex<ExecutionScheduler>>,
    experiment_id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let experiment_id: i32 = experiment_id.into_inner();
    {
        // The block drops the connection since its no longer needed and a
        // new connection is retrieved in the scheduler.
        let mut connection = database_manager.database_connection()?;
        Experiment::exists_err(experiment_id, &mut connection)?;
    }
    scheduler.lock().abort(experiment_id)?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn post_experiment_execution_reset(
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    experiment_id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let experiment_id: i32 = experiment_id.into_inner();
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(experiment_id, &mut connection)?;
    // Delete output files and logs.
    log::info!("Deleting output from experiment with ID {}.", experiment_id);
    let experiment_steps_path = app_config.experiment_steps_path(experiment_id.to_string());
    if experiment_steps_path.exists() {
        std::fs::remove_dir_all(experiment_steps_path)?;
    }
    let experiment_logs_path = app_config.experiment_logs_path(experiment_id.to_string());
    if experiment_logs_path.exists() {
        std::fs::remove_dir_all(experiment_logs_path)?;
    }
    // Delete execution steps from the database.
    connection.immediate_transaction(|connection| {
        diesel::delete(crate::schema::experiment_execution::table)
            .filter(crate::schema::experiment_execution::experiment_id.eq(experiment_id))
            .execute(connection)
    })?;

    Ok(HttpResponse::Ok().finish())
}

#[cfg(test)]
mod tests;
