#[macro_use]
extern crate diesel;

use std::{collections::HashMap, sync::Arc};

use actix_web::{middleware, App, HttpServer, web};
use application::{config::Configuration, environment::LOG_LEVEL, error::SeqError, database::DatabaseManager};
use controller::routing::routing_config;
use dotenv::dotenv;
use service::{pipeline_service::{load_pipelines, LoadedPipelines}, execution_service::ExecutionScheduler};

/// The intervall in seconds in which the pipeline execution process is updated.
const PIPELINE_EXECUTION_UPDATE_INTERVALL: u64 = 10; 

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
    let execution_config = web::Data::clone(&app_config);
    let execution_db_manager = web::Data::clone(&database_manager);
    let execution_pipelines = web::Data::clone(&loaded_pipelines);
    std::thread::spawn(move || {
        let mut scheduler = ExecutionScheduler::new(execution_config, execution_db_manager, execution_pipelines);
        loop {
            std::thread::sleep(std::time::Duration::new(PIPELINE_EXECUTION_UPDATE_INTERVALL, 0));
            if let Err(err) = scheduler.update_pipeline_execution() {
                log::error!("Updating the pipeline execution failed with error: {:?}", err);
            }
        }
    });
    // Setup the application.
    Ok(HttpServer::new(move || {
        App::new()
            .wrap(middleware::Logger::default())
            .app_data(web::Data::clone(&app_config))
            .app_data(web::Data::clone(&loaded_pipelines))
            .app_data(web::Data::clone(&database_manager))
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
