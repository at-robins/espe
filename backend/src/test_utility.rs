use crate::{application::config::Configuration, controller::routing::routing_config};
use actix_web::{
    body::MessageBody,
    dev::{ServiceFactory, ServiceRequest, ServiceResponse},
    middleware, App, Error,
};
use dotenv::dotenv;
use std::sync::Arc;

/**
 * Creates a fully configured app for testing purposes.
 */
pub fn create_test_app() -> App<
    impl ServiceFactory<
        ServiceRequest,
        Response = ServiceResponse<impl MessageBody>,
        Config = (),
        InitError = (),
        Error = Error,
    >,
> {
    dotenv().unwrap();
    env_logger::try_init_from_env(env_logger::Env::new().filter("debug")).ok();
    App::new()
        .wrap(middleware::Logger::default())
        .app_data(Arc::clone(&Arc::new(Configuration::new().unwrap())))
        .configure(routing_config)
}
