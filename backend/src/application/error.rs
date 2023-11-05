//! The `error` module defines specific error types.

use actix_web::{http::StatusCode, HttpResponse, ResponseError};
use getset::{CopyGetters, Getters};
use log::{error, warn};
use serde::{Deserialize, Serialize};
use uuid::Uuid;

pub const DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE: &str =
    "An unforseen error occurred. Please check the logs for further information.";

/// An application wide error.
#[derive(Debug, Clone, Getters, CopyGetters, Serialize, Deserialize)]
pub struct SeqError {
    /// The error ID.
    #[getset(get_copy = "pub")]
    uuid: Uuid,
    /// The error type.
    #[getset(get_copy = "pub")]
    error_type: SeqErrorType,
    /// The identifier, name or type of the error.
    #[getset(get = "pub")]
    name: String,
    /// The message for internal logging or display.
    /// This may contain sensitive data.
    #[getset(get = "pub")]
    internal_message: String,
    /// The message for external display.
    /// This must not contain sensitive data.
    #[getset(get = "pub")]
    external_message: String,
}

/// An application wide error type.
#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
pub enum SeqErrorType {
    /// A generic error implying an internal problem.
    InternalServerError,
    /// A error representing a missing resource.
    NotFoundError,
    /// A error representing an erroneous request.
    BadRequestError,
    /// A error representing a request conflicting with internal server state.
    Conflict,
}

impl SeqError {
    /// Creates a new error and automatically logs the error.
    pub fn new<T: ToString, U: ToString, V: ToString>(
        name: T,
        error_type: SeqErrorType,
        internal_message: U,
        external_message: V,
    ) -> Self {
        let error = SeqError {
            uuid: Uuid::new_v4(),
            error_type,
            name: name.to_string(),
            internal_message: internal_message.to_string(),
            external_message: external_message.to_string(),
        };
        error.log_default();
        error
    }

    /// Logs the error on its default level.
    fn log_default(&self) {
        match self.error_type() {
            SeqErrorType::InternalServerError => error!("{}", self),
            SeqErrorType::NotFoundError => warn!("{}", self),
            SeqErrorType::BadRequestError => error!("{}", self),
            SeqErrorType::Conflict => error!("{}", self),
        }
    }

    /// Returns a respective error response.
    fn error_response(&self) -> ErrorResponse {
        ErrorResponse {
            code: self.status_code().as_u16(),
            uuid: self.uuid(),
            name: self.status_code().to_string(),
            message: self.external_message().clone(),
        }
    }

    /// Converts a [`VarError`](std::env::VarError) into a [`SeqError`].
    ///
    /// # Parameters
    ///
    /// * `environment_variable` - the environment variable that caused the error
    pub fn from_var_error<T: ToString>(error: std::env::VarError, enviroment_variable: T) -> Self {
        Self::new(
            "std::env::VarError",
            SeqErrorType::InternalServerError,
            format!("{}: {}", enviroment_variable.to_string(), error),
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        )
    }
}

impl std::fmt::Display for SeqError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}: {} - {} - {}",
            self.uuid(),
            self.status_code(),
            self.name(),
            self.internal_message()
        )
    }
}

#[derive(Debug, Serialize)]
/// An informative error response for client side display
/// containing a unique ID to correlate the error with the
/// internally logged error.
/// This response must not contain sensitive data.
struct ErrorResponse {
    code: u16,
    uuid: Uuid,
    name: String,
    message: String,
}

impl ResponseError for SeqError {
    fn error_response(&self) -> HttpResponse {
        HttpResponse::build(self.status_code()).json(self.error_response())
    }

    fn status_code(&self) -> StatusCode {
        match self.error_type() {
            SeqErrorType::InternalServerError => StatusCode::INTERNAL_SERVER_ERROR,
            SeqErrorType::NotFoundError => StatusCode::NOT_FOUND,
            SeqErrorType::BadRequestError => StatusCode::BAD_REQUEST,
            SeqErrorType::Conflict => StatusCode::CONFLICT,
        }
    }
}

impl From<serde_json::Error> for SeqError {
    fn from(error: serde_json::Error) -> Self {
        Self::new(
            "serde_json::Error",
            SeqErrorType::InternalServerError,
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        )
    }
}

impl From<actix_web::error::BlockingError> for SeqError {
    fn from(error: actix_web::error::BlockingError) -> Self {
        Self::new(
            "actix_web::error::BlockingError",
            SeqErrorType::InternalServerError,
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        )
    }
}

impl From<actix_multipart::MultipartError> for SeqError {
    fn from(error: actix_multipart::MultipartError) -> Self {
        Self::new(
            "actix_multipart::MultipartError",
            SeqErrorType::InternalServerError,
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        )
    }
}

impl From<std::io::Error> for SeqError {
    fn from(error: std::io::Error) -> Self {
        Self::new(
            "std::io::Error",
            SeqErrorType::InternalServerError,
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        )
    }
}

impl From<walkdir::Error> for SeqError {
    fn from(error: walkdir::Error) -> Self {
        Self::new(
            "walkdir::Error",
            SeqErrorType::InternalServerError,
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        )
    }
}

impl From<diesel::ConnectionError> for SeqError {
    fn from(error: diesel::ConnectionError) -> Self {
        Self::new(
            "diesel::ConnectionError",
            SeqErrorType::InternalServerError,
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        )
    }
}

impl From<diesel::result::Error> for SeqError {
    fn from(error: diesel::result::Error) -> Self {
        Self::new(
            "diesel::result::Error",
            SeqErrorType::InternalServerError,
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        )
    }
}

impl From<std::time::SystemTimeError> for SeqError {
    fn from(error: std::time::SystemTimeError) -> Self {
        Self::new(
            "std::time::SystemTimeError",
            SeqErrorType::InternalServerError,
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        )
    }
}
