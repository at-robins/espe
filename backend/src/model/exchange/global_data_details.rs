use chrono::NaiveDateTime;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use crate::model::db::global_data::GlobalData;
use crate::service::multipart_service::UploadForm;

const MAX_LENGTH_FILE_PATH: usize = 128;

const ILLEGAL_COMPONENTS: [&str; 3] = [".", "..", "~"];
const ILLEGAL_COMPONENT_CHARACTERS: [&str; 42] = [
    "/", "\\", "<", ">", ":", "*", "?", "|", "\"", "\x00", "\x01", "\x02", "\x03", "\x04", "\x05",
    "\x06", "\x07", "\x08", "\x09", "\x0A", "\x0B", "\x0C", "\x0D", "\x0E", "\x0F", "\x10", "\x11",
    "\x12", "\x13", "\x14", "\x15", "\x16", "\x17", "\x18", "\x19", "\x1A", "\x1B", "\x1C", "\x1D",
    "\x1E", "\x1F", "\x7F",
];

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GlobalDataDetails {
    /// The database ID.
    pub id: i32,
    /// The display name.
    pub name: String,
    /// An optional comment on the global data repository.
    pub comment: Option<String>,
    /// The time of creation.
    pub creation_time: NaiveDateTime,
}

impl From<GlobalData> for GlobalDataDetails {
    fn from(value: GlobalData) -> Self {
        Self {
            id: value.id,
            name: value.global_data_name,
            comment: value.comment,
            creation_time: value.creation_time,
        }
    }
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GlobalDataFileDetails {
    /// The path components.
    pub path_components: Vec<String>,
    pub is_file: bool,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GlobalDataFilePath {
    /// The path components.
    pub path_components: Vec<String>,
}

impl GlobalDataFilePath {
    /// Returns the relative path to the file.
    pub fn file_path(&self) -> PathBuf {
        let mut path = PathBuf::new();
        for path_component in &self.path_components {
            path.push(path_component);
        }
        path
    }

    fn validate_file_path(&self) -> Result<(), String> {
        let file_path = self.file_path().into_os_string();
        if file_path.is_empty() {
            Err("A file path must be set.".to_string())
        } else if file_path.len() > MAX_LENGTH_FILE_PATH {
            Err(format!("The file path may only contain {} letters.", MAX_LENGTH_FILE_PATH))
        } else {
            Ok(())
        }
    }

    fn validate_file_path_component<T: AsRef<str>>(component: T) -> Result<(), String> {
        let component: &str = component.as_ref();
        if component.is_empty() {
            Err("A file path component may not be empty.".to_string())
        } else if component.len() > MAX_LENGTH_FILE_PATH {
            Err(format!(
                "The file path component may only contain {} letters.",
                MAX_LENGTH_FILE_PATH
            ))
        } else if ILLEGAL_COMPONENTS.contains(&component) {
            Err(format!("The file path component may not be \"{}\".", component))
        } else if Self::component_contains_illegal_characters(component) {
            Err(format!(
                "The file path component may not contain any of the characters \"{}\".",
                ILLEGAL_COMPONENT_CHARACTERS.join(", ")
            ))
        } else {
            Ok(())
        }
    }

    /// Returns `true` if the path component contains illegal characters.
    fn component_contains_illegal_characters<T: AsRef<str>>(component: T) -> bool {
        ILLEGAL_COMPONENT_CHARACTERS
            .into_iter()
            .any(|illegal_character| component.as_ref().contains(illegal_character))
    }

    fn validate_path_components(&self) -> Result<(), String> {
        for component in &self.path_components {
            let result = Self::validate_file_path_component(component);
            if result.is_err() {
                return result;
            }
        }
        Ok(())
    }
}

impl UploadForm for GlobalDataFilePath {
    /**
     * Checks if the recieved upload data is valid and returns a corresponding
     * error message of not.
     */
    fn validate(&self) -> Result<(), String> {
        self.validate_file_path()
            .and(self.validate_path_components())
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_file_path() {
        let file_upload = GlobalDataFilePath {
            path_components: vec![
                "a".to_string(),
                "relative".to_string(),
                "path.file".to_string(),
            ],
        };
        let path = file_upload.file_path();
        let expected_path: PathBuf = "a/relative/path.file".into();
        assert_eq!(path, expected_path);
        assert!(path.is_relative());
    }

    #[test]
    fn test_validate_valid() {
        let file_upload = GlobalDataFilePath {
            path_components: vec![
                "a".to_string(),
                "relative".to_string(),
                "path.file".to_string(),
            ],
        };
        assert!(file_upload.validate().is_ok());
    }

    #[test]
    fn test_validate_invaild_empty() {
        let file_upload = GlobalDataFilePath {
            path_components: vec!["a".to_string(), "".to_string(), "path.file".to_string()],
        };
        assert!(file_upload.validate().is_err());
    }

    #[test]
    fn test_validate_invaild_length() {
        let file_upload = GlobalDataFilePath {
            path_components: vec![
                "0123456789".to_string(),
                "0123456789".to_string(),
                "0123456789".to_string(),
                "0123456789".to_string(),
                "0123456789".to_string(),
                "0123456789".to_string(),
                "0123456789".to_string(),
                "0123456789".to_string(),
                "0123456789".to_string(),
                "0123456789".to_string(),
                "0123456789".to_string(),
                "0123456789".to_string(),
                "0123456789".to_string(),
            ],
        };
        assert!(file_upload.validate().is_err());
    }

    #[test]
    fn test_validate_invaild_dot() {
        let file_upload = GlobalDataFilePath {
            path_components: vec!["a".to_string(), ".".to_string(), "path.file".to_string()],
        };
        assert!(file_upload.validate().is_err());
    }

    #[test]
    fn test_validate_invaild_dot_dot() {
        let file_upload = GlobalDataFilePath {
            path_components: vec!["a".to_string(), "..".to_string(), "path.file".to_string()],
        };
        assert!(file_upload.validate().is_err());
    }

    #[test]
    fn test_validate_invaild_tilde() {
        let file_upload = GlobalDataFilePath {
            path_components: vec!["a".to_string(), "~".to_string(), "path.file".to_string()],
        };
        assert!(file_upload.validate().is_err());
    }
}
