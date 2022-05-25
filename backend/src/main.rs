use std::sync::Arc;

use actix_web::{middleware, App, HttpServer};
use application::{config::Configuration, error::SeqError};
use controller::routing::routing_config;
use dotenv::dotenv;

#[actix_web::main]
async fn main() -> Result<(), SeqError> {
    // Setup default enviroment variables.
    dotenv().ok();
    let app_config = Arc::new(Configuration::new()?);
    let app_config_internal = Arc::clone(&app_config);
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
