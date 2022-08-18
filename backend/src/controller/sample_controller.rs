use std::{
    io::Write,
    path::{Path, PathBuf},
    sync::Arc,
};

use crate::{
    application::config::{PATH_FILES_EXPERIMENTS, PATH_FILES_EXPERIMENT_INITIAL_FASTQ},
    diesel::ExpressionMethods,
};
use crate::{diesel::RunQueryDsl, model::db::experiment::Experiment};
use actix_multipart::Multipart;
use actix_web::{
    web::{self},
    HttpRequest, HttpResponse,
};
use diesel::QueryDsl;
use futures_util::TryStreamExt;
use log::{error, warn};
use uuid::Uuid;

use crate::{
    application::{
        config::{Configuration, PATH_FILES_TEMPORARY},
        error::{InternalError, SeqError},
    },
    model::{db::experiment::NewExperiment, exchange::experiment_upload::ExperimentUpload},
};

const MAX_MULTIPART_FORM_SIZE: usize = 524_288;

pub async fn upload_sample(
    request: HttpRequest,
    payload: Multipart,
) -> Result<HttpResponse, SeqError> {
    upload_sample_internal(request, payload)
        .await
        .map_err(|error| {
            error!("{}", error);
            error
        })
}

async fn upload_sample_internal(
    request: HttpRequest,
    mut payload: Multipart,
) -> Result<HttpResponse, SeqError> {
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");

    let connection = app_config.database_connection()?;

    // Validate that all information is provided.
    let mut is_file_provided = false;
    let mut upload_info: Option<ExperimentUpload> = None;

    // Generate a new UUID for the temporary file.
    let uuid = Configuration::generate_uuid();

    // Create the temporary folder.
    let mut temp_file_path: PathBuf = PATH_FILES_TEMPORARY.into();
    std::fs::create_dir_all(&temp_file_path)?;
    temp_file_path.push(uuid.to_string());
    let temp_file_path: Arc<Path> = temp_file_path.into();

    // Iterate over the multipart stream and save the file.
    while let Some(mut field) = payload.try_next().await? {
        match field.name() {
            "form" => {
                let mut body = web::BytesMut::new();
                while let Some(chunk) = field.try_next().await? {
                    let current_length = body.len() + chunk.len();
                    if (current_length) > MAX_MULTIPART_FORM_SIZE {
                        return Err(SeqError::BadRequestError(InternalError::new(
                            "Multipart overflow", 
                            format!("The maximum length for multipart form data is {} bytes, but the current chunk adds up to {} bytes.", MAX_MULTIPART_FORM_SIZE, current_length), 
                            "The multipart body is to large.")));
                    } else {
                        body.extend_from_slice(&chunk);
                    }
                }
                upload_info = Some(serde_json::from_slice::<ExperimentUpload>(&body)?);
            },
            "file" => {
                // Write the file to disk.
                let file_path_ref = Arc::clone(&temp_file_path);
                let mut file = web::block(|| std::fs::File::create(file_path_ref)).await??;
                while let Some(chunk) = field.try_next().await? {
                    file = web::block(move || file.write_all(&chunk).map(|_| file)).await??;
                }
                is_file_provided = true;
            },
            name => warn!("Unknown content name: {}", name),
        }
    }

    if upload_info.is_some() && is_file_provided {
        let upload_info = upload_info.unwrap();
        upload_info.validate().map_err(|e| {
            SeqError::BadRequestError(InternalError::new(
                "Sample invalid",
                format!(
                    "Validation of sample information {:?} failed with error: {}",
                    upload_info, e
                ),
                "The provided sample information is invalid.",
            ))
        })?;
        let new_experiment: NewExperiment = upload_info.into();
        // Write to database.
        diesel::insert_into(crate::schema::experiment::table)
            .values(&new_experiment)
            .execute(&connection)?;

        // Get ID of the inserted record.
        // The corner case of multiple identical samples being inserted at the same time using the same pipeline is not covered.
        let inserted = crate::schema::experiment::table
            .filter(crate::schema::experiment::experiment_name.eq(new_experiment.experiment_name()))
            .filter(crate::schema::experiment::creation_time.eq(new_experiment.creation_time()))
            .filter(crate::schema::experiment::pipeline_id.eq(new_experiment.pipeline_id()))
            .first::<Experiment>(&connection)?;
        let inserted_id = inserted.id;

        // Create experiment folder and copy file to destination.
        let mut final_file_path: PathBuf = PATH_FILES_EXPERIMENTS.into();
        final_file_path.push(inserted_id.to_string());
        std::fs::create_dir_all(&final_file_path)?;
        final_file_path.push(PATH_FILES_EXPERIMENT_INITIAL_FASTQ);
        std::fs::rename(temp_file_path, final_file_path)?;

        // Return the ID of the created attachment.
        Ok(HttpResponse::Created().body(inserted_id.to_string()))
    } else {
        delete_temporary_file(uuid)?;
        Err(SeqError::BadRequestError(InternalError::new("Missing multipart data",
         format!("For the experiment upload both a file ({}) and the according form data ({:?}) must be present.", is_file_provided, upload_info),
         "For the experiment upload both a file and the according form data must be present.")))
    }
}

/// Deletes the temporary file associated with the specified UUID.
///
/// # Parameters
///
/// `uuid` - the UUID of the temporary file
fn delete_temporary_file(uuid: Uuid) -> Result<(), SeqError> {
    let mut file_path: PathBuf = PATH_FILES_TEMPORARY.into();
    file_path.push(uuid.to_string());
    if file_path.exists() {
        std::fs::remove_file(file_path)?;
    } else {
        warn!("Tried to delete non existing temporary file {:?}.", file_path)
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::test_utility::create_test_app;

    use super::*;
    use actix_web::{
        http::{header::ContentType, StatusCode},
        test,
    };
    use mime;

    #[actix_web::test]
    async fn test_upload_sample_post() {
        let app = test::init_service(create_test_app()).await;
        let payload =
            std::fs::read("../testing_resources/requests/sample_submission_multipart").unwrap();
        let content_type: mime::Mime =
            "multipart/form-data; boundary=---------------------------5851692324164894962235391524"
                .parse()
                .unwrap();
        let req = test::TestRequest::post()
            .uri("/api/experiment")
            .insert_header(ContentType(content_type))
            .set_payload(payload)
            .to_request();
        error!("{:?}", req);
        let resp = test::call_service(&app, req).await;
        assert_eq!(resp.status(), StatusCode::OK);
    }
}
