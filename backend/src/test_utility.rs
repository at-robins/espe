use crate::{application::config::Configuration, controller::routing::routing_config};
use actix_web::{
    body::MessageBody,
    dev::{ServiceFactory, ServiceRequest, ServiceResponse},
    middleware, App, Error,
};
use diesel::{Connection, SqliteConnection};
use dotenv::dotenv;
use std::{path::PathBuf, sync::Arc};
use uuid::Uuid;

const TEST_RESOURCES_PATH: &str = "../testing_resources";

/// Creates a fully configured app for testing purposes
/// for testing against the specified test database.
///
/// # Parameters
///
/// * `database_url` - the URL / URI of the test database
pub fn create_test_app(
    context: &TestContext,
) -> App<
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
        .app_data(Arc::clone(&Arc::<Configuration>::new(context.into())))
        .configure(routing_config)
}

/// A test context that provides clean test resources (e.g. test databases) and according
/// initialisation and cleanup on a per test basis.
pub struct TestContext {
    id: Uuid,
}

impl TestContext {
    /// Creates a new `TestContext`.
    pub fn new() -> TestContext {
        let id = Configuration::generate_uuid();
        let context = TestContext { id };
        std::fs::create_dir_all(context.context_folder()).unwrap();
        diesel_migrations::run_pending_migrations(&context.get_connection()).unwrap();
        context
    }

    /// Returns the URL / URI of the test database.
    pub fn database_url(&self) -> String {
        format!("{}/{}.db", self.context_folder(), self.id)
    }

    /// Returns the context folder that stores all context information.
    pub fn context_folder(&self) -> String {
        format!("{}/{}", TEST_RESOURCES_PATH, self.id)
    }

    /// Returns the pipeline definition folder.
    pub fn pipeline_folder(&self) -> String {
        format!("{}/{}_pipelines", self.context_folder(), self.id)
    }

    /// Opens a connection to the test database.
    pub fn get_connection(&self) -> SqliteConnection {
        let connection = SqliteConnection::establish(&self.database_url()).unwrap();
        connection.execute("PRAGMA foreign_keys = ON;").unwrap();
        connection
    }
}

impl Drop for TestContext {
    fn drop(&mut self) {
        let context_path: PathBuf = self.context_folder().into();
        if context_path.exists() {
            if let Err(e) = std::fs::remove_dir_all(&context_path) {
                log::error!("Dropping test context {} failed with error: {}", self.id, e);
            } else {
                log::info!("Removed test context {}.", context_path.display());
            }
        } else {
            log::warn!("Tried to delete non existing temporary context folder {}.", context_path.display());
        }
    }
}

impl From<&TestContext> for Configuration {
    fn from(context: &TestContext) -> Self {
        Configuration::new(
            context.database_url(),
            "info",
            "127.0.0.1",
            "8080",
            context.context_folder(),
            context.pipeline_folder(),
        )
    }
}

impl From<TestContext> for Configuration {
    fn from(context: TestContext) -> Self {
        (&context).into()
    }
}