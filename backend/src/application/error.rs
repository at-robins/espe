//! The `error` module defines specific error types.

use actix_web::{http::StatusCode, HttpResponse, ResponseError};
use getset::{CopyGetters, Getters};
use serde::{Deserialize, Serialize};
use uuid::Uuid;

pub const DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE: &str =
    "An unforseen error occurred. Please check the logs for further information.";

/// An application wide error type.
#[derive(Debug)]
pub enum SeqError {
    /// A generic error implying an internal problem.
    InternalServerError(InternalError),
    /// A error representing a missing resource.
    NotFoundError(InternalError),
    /// A error representing an erroneous request.
    BadRequestError(InternalError),
}

#[derive(Debug, Clone, Getters, CopyGetters, Serialize, Deserialize)]
pub struct InternalError {
    /// The error ID.
    #[getset(get_copy = "pub")]
    uuid: Uuid,
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

impl InternalError {
    /// Creates a new internal error.
    pub fn new<T: ToString, U: ToString, V: ToString>(
        name: T,
        internal_message: U,
        external_message: V,
    ) -> Self {
        InternalError {
            uuid: Uuid::new_v4(),
            name: name.to_string(),
            internal_message: internal_message.to_string(),
            external_message: external_message.to_string(),
        }
    }
}

impl std::fmt::Display for InternalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {} - {}", self.uuid(), self.name(), self.internal_message())
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

impl SeqError {
    fn error_response(&self) -> ErrorResponse {
        match self {
            Self::InternalServerError(internal) => ErrorResponse {
                code: self.status_code().as_u16(),
                uuid: internal.uuid(),
                name: self.status_code().to_string(),
                message: internal.external_message().clone(),
            },
            Self::NotFoundError(internal) => ErrorResponse {
                code: self.status_code().as_u16(),
                uuid: internal.uuid(),
                name: self.status_code().to_string(),
                message: internal.external_message().clone(),
            },
            Self::BadRequestError(internal) => ErrorResponse {
                code: self.status_code().as_u16(),
                uuid: internal.uuid(),
                name: self.status_code().to_string(),
                message: internal.external_message().clone(),
            },
        }
    }
}

impl std::fmt::Display for SeqError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InternalServerError(internal) => write!(f, "{}", internal),
            Self::NotFoundError(internal) => write!(f, "{}", internal),
            Self::BadRequestError(internal) => write!(f, "{}", internal),
        }
    }
}

impl ResponseError for SeqError {
    fn error_response(&self) -> HttpResponse {
        HttpResponse::build(self.status_code()).json(self.error_response())
    }

    fn status_code(&self) -> StatusCode {
        match self {
            Self::InternalServerError(_) => StatusCode::INTERNAL_SERVER_ERROR,
            Self::NotFoundError(_) => StatusCode::NOT_FOUND,
            Self::BadRequestError(_) => StatusCode::BAD_REQUEST,
        }
    }
}

impl From<serde_json::Error> for SeqError {
    fn from(error: serde_json::Error) -> Self {
        Self::InternalServerError(InternalError::new(
            "serde_json::Error",
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        ))
    }
}

impl From<actix_web::error::BlockingError> for SeqError {
    fn from(error: actix_web::error::BlockingError) -> Self {
        Self::InternalServerError(InternalError::new(
            "actix_web::error::BlockingError",
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        ))
    }
}

impl From<actix_multipart::MultipartError> for SeqError {
    fn from(error: actix_multipart::MultipartError) -> Self {
        Self::InternalServerError(InternalError::new(
            "actix_multipart::MultipartError",
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        ))
    }
}

impl From<std::io::Error> for SeqError {
    fn from(error: std::io::Error) -> Self {
        Self::InternalServerError(InternalError::new(
            "std::io::Error",
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        ))
    }
}

impl From<std::env::VarError> for SeqError {
    fn from(error: std::env::VarError) -> Self {
        Self::InternalServerError(InternalError::new(
            "std::env::VarError",
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        ))
    }
}

impl From<diesel::ConnectionError> for SeqError {
    fn from(error: diesel::ConnectionError) -> Self {
        Self::InternalServerError(InternalError::new(
            "diesel::ConnectionError",
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        ))
    }
}

impl From<diesel::result::Error> for SeqError {
    fn from(error: diesel::result::Error) -> Self {
        Self::InternalServerError(InternalError::new(
            "diesel::result::Error",
            error,
            DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
        ))
    }
}
