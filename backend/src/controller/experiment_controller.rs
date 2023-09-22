use crate::{
    application::{
        config::Configuration,
        error::{SeqError, SeqErrorType},
    },
    diesel::RunQueryDsl,
    model::{
        db::{
            experiment::{Experiment, NewExperiment},
            experiment_execution::{ExperimentExecution, NewExperimentExecution, ExecutionStatus},
            pipeline_step_variable::{NewPipelineStepVariable, PipelineStepVariable},
        },
        exchange::{
            experiment_details::ExperimentDetails,
            experiment_pipeline::ExperimentPipelineBlueprint,
            pipeline_variable_upload::PipelineStepVariableUpload,
        },
    },
    service::{
        pipeline_service::LoadedPipelines,
        validation_service::{validate_comment, validate_entity_name, validate_mail},
    },
};
use actix_web::{web, HttpResponse};
use diesel::{ExpressionMethods, QueryDsl};

pub async fn create_experiment(
    app_config: web::Data<Configuration>,
    name: actix_web::web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let name: String = name.into_inner();
    validate_entity_name(&name)?;
    let mut connection = app_config.database_connection()?;
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
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Retrieve the app config.
    let mut connection = app_config.database_connection()?;
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
    app_config: web::Data<Configuration>,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    let experiment_details: ExperimentDetails = crate::schema::experiment::table
        .find(id)
        .first::<Experiment>(&mut connection)?
        .into();
    Ok(HttpResponse::Ok().json(experiment_details))
}

pub async fn get_experiment_execution_status(
    app_config: web::Data<Configuration>,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    let execution_steps = ExperimentExecution::get_by_experiment(id, &mut connection)?;
    let result = if execution_steps.is_empty() {
        "None".to_string()
    } else if execution_steps.iter().any(|execution| execution.execution_status == ExecutionStatus::Failed.to_string()) {
        ExecutionStatus::Failed.to_string()
    } else if execution_steps.iter().any(|execution| execution.execution_status == ExecutionStatus::Aborted.to_string()) {
        ExecutionStatus::Aborted.to_string()
    } else if execution_steps.iter().any(|execution| execution.execution_status == ExecutionStatus::Running.to_string()) {
        ExecutionStatus::Running.to_string()
    } else if execution_steps.iter().all(|execution| execution.execution_status == ExecutionStatus::Finished.to_string()) {
        ExecutionStatus::Finished.to_string()
    } else {
        ExecutionStatus::Waiting.to_string()
    };
    Ok(HttpResponse::Ok().json(result))
}

pub async fn patch_experiment_name(
    app_config: web::Data<Configuration>,
    id: web::Path<i32>,
    new_name: web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_name = new_name.into_inner();
    validate_entity_name(&new_name)?;
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::experiment_name.eq(new_name))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn patch_experiment_mail(
    app_config: web::Data<Configuration>,
    id: web::Path<i32>,
    new_name: web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_mail = new_name.into_inner();
    validate_mail(&new_mail)?;
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::mail.eq(new_mail))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn patch_experiment_comment(
    app_config: web::Data<Configuration>,
    id: web::Path<i32>,
    new_comment: web::Json<Option<String>>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Sanitise the HTML and validate.
    let new_comment = new_comment.into_inner().map(|inner| ammonia::clean(&inner));
    if let Some(inner) = &new_comment {
        validate_comment(inner)?;
    }
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::comment.eq(new_comment))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn patch_experiment_pipeline(
    app_config: web::Data<Configuration>,
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
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::pipeline_id.eq(new_pipeline))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn list_experiment(
    app_config: web::Data<Configuration>,
) -> Result<web::Json<Vec<ExperimentDetails>>, SeqError> {
    let mut connection = app_config.database_connection()?;
    let experiments: Vec<ExperimentDetails> = Experiment::get_all(&mut connection)?
        .into_iter()
        .map(|val| val.into())
        .collect();
    Ok(web::Json(experiments))
}

pub async fn get_experiment_pipelines(
    app_config: web::Data<Configuration>,
    pipelines: web::Data<LoadedPipelines>,
    id: web::Path<i32>,
) -> Result<web::Json<Vec<ExperimentPipelineBlueprint>>, SeqError> {
    let experiment_id: i32 = id.into_inner();
    let mut connection = app_config.database_connection()?;
    let mut experiment_pipelines = Vec::new();
    for pipeline in pipelines.pipelines() {
        let values = crate::model::db::pipeline_step_variable::PipelineStepVariable::get_values_by_experiment_and_pipeline(experiment_id, pipeline.pipeline().id(), &mut connection)?;
        experiment_pipelines
            .push(ExperimentPipelineBlueprint::from_internal(pipeline.pipeline(), values));
    }
    Ok(web::Json(experiment_pipelines))
}

pub async fn post_experiment_pipeline_variable(
    app_config: web::Data<Configuration>,
    pipelines: web::Data<LoadedPipelines>,
    experiment_id: web::Path<i32>,
    new_variable: web::Json<PipelineStepVariableUpload>,
) -> Result<HttpResponse, SeqError> {
    let experiment_id: i32 = experiment_id.into_inner();
    let new_variable: PipelineStepVariableUpload = new_variable.into_inner();
    if !pipelines.has_variable(
        &new_variable.pipeline_id,
        &new_variable.pipeline_step_id,
        &new_variable.variable_id,
    ) {
        return Err(SeqError::new(
            "Not Found",
            SeqErrorType::NotFoundError,
            format!(
                "No pipeline variable with corresponding properties is currently loaded, thus variable {:?} cannot be inserted.",
                new_variable
            ),
            "The variable is invalid.",
        ));
    }
    let mut connection = app_config.database_connection()?;
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

pub async fn post_execute_experiment(
    app_config: web::Data<Configuration>,
    experiment_id: web::Path<i32>,
    pipelines: web::Data<LoadedPipelines>,
) -> Result<HttpResponse, SeqError> {
    let experiment_id: i32 = experiment_id.into_inner();
    let mut connection = app_config.database_connection()?;
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
            let experiment_variables = PipelineStepVariable::get_values_by_experiment_and_pipeline(
                experiment_id,
                pipeline_id,
                &mut connection,
            )?;
            let mut execution_steps: Vec<NewExperimentExecution> =
                Vec::with_capacity(pipeline.steps().len());
            for step in pipeline.steps() {
                execution_steps.push(NewExperimentExecution::new(
                    experiment_id,
                    pipeline.id(),
                    step.id(),
                ));
                for variable in step.variables() {
                    if variable.required().unwrap_or(false) {
                        // Error if required variables are not set.
                        if !experiment_variables.contains_key(&format!(
                            "{}{}",
                            step.id(),
                            variable.id()
                        )) {
                            return Err(SeqError::new(
                                "Invalid run",
                                SeqErrorType::BadRequestError,
                                format!("The experiment {} is missing the required variable with pipeline id {} step id {} and variable id {}.", experiment_id, pipeline.id(), step.id(), variable.id()),
                                "The requested run parameters are invalid.",
                            ));
                        }
                    }
                }
            }
            diesel::insert_into(crate::schema::experiment_execution::table)
                .values(&execution_steps)
                .execute(&mut connection)?;
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
                "No pipeline was selected for experiment {}, so i cannot be run.",
                experiment_id
            ),
            "The requested run parameters are invalid.",
        ))
    }
}

#[cfg(test)]
mod tests;
