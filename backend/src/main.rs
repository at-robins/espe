#[macro_use]
extern crate diesel;

use std::path::PathBuf;

use actix_web::{dev::Service, middleware, web::Data, App, HttpServer};
use application::{
    config::Configuration,
    database::DatabaseManager,
    environment::LOG_LEVEL,
    error::{
        SeqError, SeqErrorLogger, SeqErrorType, DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
    },
};
use controller::routing::routing_config;
use diesel_migrations::{
    embed_migrations, EmbeddedMigrations, HarnessWithOutput, MigrationHarness,
};
use dotenv::dotenv;
use futures::FutureExt;
use parking_lot::Mutex;
use service::{
    execution_service::ExecutionScheduler, pipeline_service::LoadedPipelines,
    temp_file_service::TemporaryFileManager,
};

/// The intervall in seconds in which the pipeline execution process is updated.
const PIPELINE_EXECUTION_UPDATE_INTERVALL: u64 = 5;
/// The intervall in seconds in which temporary data is inspected.
const TEMPORARY_DATA_MANAGEMENT_UPDATE_INTERVALL: u64 = 300;
/// The compiled database migrations.
pub const MIGRATIONS: EmbeddedMigrations = embed_migrations!("./migrations");

/// Sets up the environment and the logger.
fn setup_environment() -> Result<(), SeqError> {
    // Setup default enviroment variables.
    let environment_setup_result = dotenv()
        .map_err(|err| SeqError::from(err).chain("Setting up the application environment failed."));
    // Setup the logger.
    env_logger::init_from_env(env_logger::Env::new().filter(LOG_LEVEL));
    environment_setup_result
        .map(|_| ())
        .map_err(|err| SeqError::from(err).chain("Setting up the application environment failed."))
}

/// Sets up the application's [`Configuration`].
fn setup_configuration() -> Result<Data<Configuration>, SeqError> {
    Configuration::create_from_environment()
        .map(|app_config| Data::new(app_config))
        .map_err(|err| {
            SeqError::from(err).chain("Setting up the application's configuration failed.")
        })
}

/// Sets up the database.
///
/// # Parameters
///
/// * `app_config` - the application's [`Configuration`]
fn setup_database(app_config: &Data<Configuration>) -> Result<Data<DatabaseManager>, SeqError> {
    // Setup database and conncetion pool.
    if let Some(database_path) = PathBuf::from(app_config.database_url()).parent() {
        std::fs::create_dir_all(database_path).map_err(|err| {
            SeqError::from(err)
                .chain(format!("Could not create database directory {}", database_path.display()))
        })?;
    }
    let database_manager =
        Data::new(DatabaseManager::new(Data::clone(&app_config)).map_err(|err| {
            err.chain("Could not create database manager during application setup.")
        })?);
    {
        let mut connection = database_manager.database_connection().map_err(|err| {
            err.chain("Could not obtain a database connection during application setup.")
        })?;
        HarnessWithOutput::write_to_stdout(&mut connection)
            .run_pending_migrations(MIGRATIONS)
            .map_err(|error| {
                SeqError::new(
                    "diesel_migrations::MigrationError",
                    SeqErrorType::InternalServerError,
                    error,
                    DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
                )
                .chain("Could not apply database migrations during application setup.")
            })?;
    }
    return Ok(database_manager);
}

/// Loads the pipelines into memory.
///
/// # Parameters
///
/// * `app_config` - the application's [`Configuration`]
fn setup_pipelines(app_config: &Data<Configuration>) -> Result<Data<LoadedPipelines>, SeqError> {
    let (loaded_pipelines, errors) = match LoadedPipelines::new(Data::clone(&app_config)) {
        Ok(loaded_pipelines) => (loaded_pipelines, Vec::new()),
        Err((errors, loaded_pipelines)) => (
            loaded_pipelines,
            errors
                .into_iter()
                .map(|err| err.chain("Loading of a pipeline failed during application setup."))
                .collect(),
        ),
    };

    // Logs all errors during pipeline loading.
    // These might be minor errors that are only relevant for
    // specific pipelines. Others might have been loaded fine.
    for error in &errors {
        error.log_default();
    }

    if loaded_pipelines.size() == 0 && !errors.is_empty() {
        // Something went really wrong. There were errors and not a single pipeline could be loaded.
        Err(SeqError::new(
            "Pipeline setup error",
            SeqErrorType::InternalServerError,
            "Loading the pipelines failed during application setup. Please check the pipeline specific logging statements for more detailed information on the errors.",
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        ))
    } else {
        Ok(Data::new(loaded_pipelines))
    }
}

/// Starts and returns the scheduler that handles pipeline execution.
///
/// # Parameters
///
/// * `app_config` - the application's [`Configuration`]
/// * `database_manager` - the application's [`DatabaseManager`]
/// * `loaded_pipelines` - the application's [`LoadedPipelines`]
fn setup_execution_scheduler(
    app_config: &Data<Configuration>,
    database_manager: &Data<DatabaseManager>,
    loaded_pipelines: &Data<LoadedPipelines>,
) -> Data<Mutex<ExecutionScheduler>> {
    let scheduler = Data::new(Mutex::new(ExecutionScheduler::new(
        Data::clone(&app_config),
        Data::clone(&database_manager),
        Data::clone(&loaded_pipelines),
    )));
    let execution_scheduler = Data::clone(&scheduler);
    std::thread::spawn(move || loop {
        std::thread::sleep(std::time::Duration::new(PIPELINE_EXECUTION_UPDATE_INTERVALL, 0));
        if let Err(err) = execution_scheduler.lock().update_pipeline_execution() {
            err.chain("Updating the pipeline execution failed.")
                .log_default();
        }
    });

    scheduler
}

/// Starts the temporary file management.
///
/// # Parameters
///
/// * `app_config` - the application's [`Configuration`]
fn setup_temporary_file_manager_scheduler(app_config: &Data<Configuration>) {
    let temp_file_manager_config = Data::clone(&app_config);
    std::thread::spawn(move || {
        let temp_file_manager = TemporaryFileManager::new(temp_file_manager_config);
        loop {
            if let Err(err) = temp_file_manager.update() {
                err.chain("Managing temporary data failed.").log_default();
            }
            std::thread::sleep(std::time::Duration::new(
                TEMPORARY_DATA_MANAGEMENT_UPDATE_INTERVALL,
                0,
            ));
        }
    });
}

/// Starts the application.
async fn start_application() -> Result<(), SeqError> {
    let app_config = setup_environment().and(setup_configuration())?;
    let server_address = app_config.server_address_and_port();
    let database_manager = setup_database(&app_config)?;
    let loaded_pipelines = setup_pipelines(&app_config)?;
    let scheduler = setup_execution_scheduler(&app_config, &database_manager, &loaded_pipelines);
    setup_temporary_file_manager_scheduler(&app_config);

    // Setup the application.
    Ok(HttpServer::new(move || {
        App::new()
            .wrap(middleware::Logger::default())
            .wrap_fn(|request, srv| {
                // Logs internal errrors.
                srv.call(request).map(|response_result| {
                    response_result.and_then(|response| {
                        if let Some(seq_error_logger) =
                            response.response().extensions().get::<SeqErrorLogger>()
                        {
                            seq_error_logger.log_default();
                        }
                        Ok(response)
                    })
                })
            })
            .app_data(Data::clone(&app_config))
            .app_data(Data::clone(&loaded_pipelines))
            .app_data(Data::clone(&database_manager))
            .app_data(Data::clone(&scheduler))
            .configure(routing_config)
    })
    .bind(server_address)?
    .run()
    .await?)
}

#[actix_web::main]
async fn main() -> Result<(), SeqError> {
    start_application().await.map_err(|err| {
        let final_err = err.chain("Fatal error while running the application.");
        final_err.log_default();
        final_err
    })
}

mod application;
mod controller;
mod model;
mod schema;
mod service;
#[cfg(test)]
mod test_utility;
