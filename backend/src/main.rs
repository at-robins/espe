#[macro_use]
extern crate diesel;

use std::{collections::HashMap, sync::Arc};

use actix_web::{middleware, web, App, HttpServer};
use application::{
    config::Configuration, database::DatabaseManager, environment::LOG_LEVEL, error::SeqError,
};
use controller::routing::routing_config;
use dotenv::dotenv;
use parking_lot::Mutex;
use service::{
    execution_service::ExecutionScheduler,
    pipeline_service::{load_pipelines, LoadedPipelines},
    temp_file_service::TemporaryFileManager,
};

/// The intervall in seconds in which the pipeline execution process is updated.
const PIPELINE_EXECUTION_UPDATE_INTERVALL: u64 = 10;
/// The intervall in seconds in which temporary data is inspected.
const TEMPORARY_DATA_MANAGEMENT_UPDATE_INTERVALL: u64 = 300;

#[actix_web::main]
async fn main() -> Result<(), SeqError> {
    // Setup default enviroment variables.
    let environment_setup_result = dotenv();
    // Setup the logger.
    env_logger::init_from_env(env_logger::Env::new().filter(LOG_LEVEL));
    // Log potential errors that occurred during environment setup.
    if let Err(enviroment_error) = environment_setup_result {
        log::error!("{}", enviroment_error);
    }
    // Setup the configuration.
    let app_config = web::Data::new(Configuration::create_from_environment()?);
    let server_address = app_config.server_address_and_port();
    // Setup database conncetion pool.
    let database_manager = web::Data::new(DatabaseManager::new(web::Data::clone(&app_config))?);
    // Load all pipelines into memory.
    let pipelines = load_pipelines(Arc::clone(&app_config))?;
    let mut pipeline_map = HashMap::new();
    for pipeline in pipelines {
        let duplicate = pipeline_map.insert(pipeline.pipeline().id().clone(), pipeline);
        if let Some(duplicate_pipeline) = duplicate {
            log::warn!(
                "The pipeline {:?} was overwritten due to pipeline ID {} not being unique.",
                duplicate_pipeline,
                duplicate_pipeline.pipeline().id()
            );
        }
    }
    let loaded_pipelines = web::Data::new(LoadedPipelines::new(web::Data::clone(&app_config))?);
    let scheduler = web::Data::new(Mutex::new(ExecutionScheduler::new(
        web::Data::clone(&app_config),
        web::Data::clone(&database_manager),
        web::Data::clone(&loaded_pipelines),
    )));
    let execution_scheduler = web::Data::clone(&scheduler);
    std::thread::spawn(move || loop {
        std::thread::sleep(std::time::Duration::new(PIPELINE_EXECUTION_UPDATE_INTERVALL, 0));
        if let Err(err) = execution_scheduler.lock().update_pipeline_execution() {
            log::error!("Updating the pipeline execution failed with error: {:?}", err);
        }
    });
    // Setup temporary file management.
    let temp_file_manager_config = web::Data::clone(&app_config);
    std::thread::spawn(move || {
        let temp_file_manager = TemporaryFileManager::new(temp_file_manager_config);
        loop {
            if let Err(err) = temp_file_manager.update() {
                log::error!("Managing temporary data failed with error: {:?}", err);
            }
            std::thread::sleep(std::time::Duration::new(
                TEMPORARY_DATA_MANAGEMENT_UPDATE_INTERVALL,
                0,
            ));
        }
    });
    // Setup the application.
    Ok(HttpServer::new(move || {
        App::new()
            .wrap(middleware::Logger::default())
            .app_data(web::Data::clone(&app_config))
            .app_data(web::Data::clone(&loaded_pipelines))
            .app_data(web::Data::clone(&database_manager))
            .app_data(web::Data::clone(&scheduler))
            .configure(routing_config)
    })
    .bind(server_address)?
    .run()
    .await?)
}

mod application;
mod controller;
mod model;
mod schema;
mod service;
#[cfg(test)]
mod test_utility;
