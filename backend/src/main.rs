#[macro_use]
extern crate diesel;
#[macro_use]
extern crate lazy_static;

use std::sync::Arc;

use actix_web::{middleware, App, HttpServer};
use application::{config::Configuration, environment::LOG_LEVEL, error::SeqError};
use controller::routing::routing_config;
use dotenv::dotenv;

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
    let app_config = Arc::new(Configuration::new()?);
    let app_config_internal = Arc::clone(&app_config);
    // Setup the application.
    Ok(HttpServer::new(move || {
        App::new()
            .wrap(middleware::Logger::default())
            .app_data(Arc::clone(&app_config_internal))
            .configure(routing_config)
    })
    .bind(app_config.server_address_and_port())?
    .run()
    .await?)
}

mod application;
mod controller;
mod model;
mod schema;
#[cfg(test)]
mod test_utility;
