use std::{
    path::{Path, PathBuf},
    sync::Arc,
};

use crate::{
    application::{
        config::Configuration, config::PATH_FILES_EXPERIMENT_INITIAL_FASTQ, error::SeqError,
        error::SeqErrorType,
    },
    diesel::{ExpressionMethods, RunQueryDsl},
    model::{db::experiment::NewExperiment, exchange::experiment_upload::ExperimentUpload},
    service::multipart_service::{
        create_temporary_file, delete_temporary_file, parse_multipart_file,
    },
};
use actix_multipart::Multipart;
use actix_web::{HttpRequest, HttpResponse};
use diesel::{dsl::exists, QueryDsl};

pub async fn upload_sample(
    request: HttpRequest,
    payload: Multipart,
) -> Result<HttpResponse, SeqError> {
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");

    let (temporary_file_path, temporary_file_id) = create_temporary_file(Arc::clone(&app_config))?;

    let inserted_id =
        persist_multipart(payload, temporary_file_path.as_path(), Arc::clone(app_config))
            .await
            .map_err(|error| {
                // Delete temporary file on error.
                if delete_temporary_file(temporary_file_id, Arc::clone(app_config)).is_err() {
                    log::error!(
                        "Failed to delete temporary file {} upon error {}.",
                        temporary_file_path.display(),
                        error
                    );
                }
                error
            })?;

    // Return the ID of the created attachment.
    Ok(HttpResponse::Created().body(inserted_id.to_string()))
}

async fn persist_multipart<P: AsRef<Path>>(
    payload: Multipart,
    temporary_file_path: P,
    app_config: Arc<Configuration>,
) -> Result<i32, SeqError> {
    let mut connection = app_config.database_connection()?;

    let (upload_info, temp_file_path) =
        parse_multipart_file::<ExperimentUpload, P>(payload, temporary_file_path).await?;

    // Validate the existance of the pipeline.
    let pipeline_exists: bool = diesel::select(exists(
        crate::schema::pipeline::dsl::pipeline
            .filter(crate::schema::pipeline::id.eq(upload_info.pipeline_id)),
    ))
    .get_result(&mut connection)?;
    if !pipeline_exists {
        return Err(SeqError::new(
                "Pipeline invalid",
                SeqErrorType::BadRequestError,
                format!(
                    "Validation of sample information {:?} failed with error: Pipeline with ID {} does not exist.",
                    upload_info, upload_info.pipeline_id
                ),
                "The provided sample information is invalid.",
            ));
    }

    // Write experiment to database.
    let new_experiment: NewExperiment = upload_info.into();
    let inserted_id = diesel::insert_into(crate::schema::experiment::table)
        .values(&new_experiment)
        .returning(crate::schema::experiment::id)
        .get_result(&mut connection)?;

    // Create experiment folder and copy file to destination.
    temp_file_to_experiment(inserted_id, temp_file_path, app_config).map_err(|error| {
        // Roll back database if there is an error while moving the temporary file.
        if let Err(e) = diesel::delete(
            crate::schema::experiment::dsl::experiment
                .filter(crate::schema::experiment::id.eq(inserted_id)),
        )
        .execute(&mut connection)
        {
            log::error!(
                "Roll back of database after insertion of experiment {} failed with error: {}.",
                inserted_id,
                e
            );
        }
        error
    })?;

    // Return the ID of the created attachment.
    Ok(inserted_id)
}

fn temp_file_to_experiment<P: AsRef<Path>>(
    inserted_id: i32,
    temp_file_path: P,
    app_config: Arc<Configuration>,
) -> Result<(), SeqError> {
    let mut final_file_path: PathBuf = app_config.experiments_path();
    final_file_path.push(inserted_id.to_string());
    std::fs::create_dir_all(&final_file_path)?;
    final_file_path.push(PATH_FILES_EXPERIMENT_INITIAL_FASTQ);
    std::fs::rename(temp_file_path, final_file_path)?;
    Ok(())
}

#[cfg(test)]
mod tests;
