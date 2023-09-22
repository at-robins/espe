use crate::{
    application::config::Configuration, controller::routing::routing_config,
    service::pipeline_service::LoadedPipelines,
};
use actix_web::{
    body::MessageBody,
    dev::{ServiceFactory, ServiceRequest, ServiceResponse},
    middleware, web, App, Error,
};
use diesel::{connection::SimpleConnection, Connection, SqliteConnection};
use diesel_migrations::MigrationHarness;
use dotenv::dotenv;
use std::path::PathBuf;
use uuid::Uuid;

/// The test resource path.
pub const TEST_RESOURCES_PATH: &str = "../testing_resources";

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
    let app_config = &web::Data::<Configuration>::new(context.into());
    App::new()
        .wrap(middleware::Logger::default())
        .app_data(web::Data::clone(&app_config))
        .app_data(web::Data::new(LoadedPipelines::new(web::Data::clone(&app_config)).unwrap()))
        .configure(routing_config)
}

/// A test context that provides clean test resources (e.g. test databases) and according
/// initialisation and cleanup on a per test basis.
pub struct TestContext {
    id: Uuid,
    pipeline_folder_override: Option<String>,
}

impl TestContext {
    /// Creates a new `TestContext`.
    pub fn new() -> TestContext {
        let id = Uuid::new_v4();
        let context = TestContext { id, pipeline_folder_override: None };
        std::fs::create_dir_all(context.context_folder()).unwrap();
        std::fs::create_dir_all(context.pipeline_folder()).unwrap();
        let con = &mut context.get_connection();
        con.run_pending_migrations(
            diesel_migrations::FileBasedMigrations::find_migrations_directory().unwrap(),
        )
        .unwrap();
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
        if let Some(custom_folder) = &self.pipeline_folder_override {
            custom_folder.clone()
        } else {
            format!("{}/{}_pipelines", self.context_folder(), self.id)
        }
    }

    /// Overrides the default testing pipeline folder with a custom one.
    /// 
    /// # Parameters
    /// 
    /// * `pipeline_folder` - the custom testing pipeline folder
    pub fn set_pipeline_folder<T: Into<String>>(&mut self, pipeline_folder: T) {
        self.pipeline_folder_override = Some(pipeline_folder.into());
    }

    /// Opens a connection to the test database.
    pub fn get_connection(&self) -> SqliteConnection {
        let mut connection = SqliteConnection::establish(&self.database_url()).unwrap();
        connection
            .batch_execute(
                "PRAGMA foreign_keys = ON;
            PRAGMA journal_mode = WAL;
            PRAGMA synchronous = NORMAL;
            PRAGMA busy_timeout = 10000;",
            )
            .unwrap();
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
            log::warn!(
                "Tried to delete non existing temporary context folder {}.",
                context_path.display()
            );
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
