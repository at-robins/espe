use std::path::PathBuf;

use chrono::NaiveDateTime;
use serde::{Deserialize, Serialize};

use crate::model::db::global_data::{GlobalData, GlobalDataFile};

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

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GlobalDataFileDetails {
    /// The database ID.
    pub id: i32,
    /// The path components.
    pub path_components: Vec<String>,
    /// The time of creation.
    pub creation_time: NaiveDateTime,
}

impl From<GlobalDataFile> for GlobalDataFileDetails {
    fn from(value: GlobalDataFile) -> Self {
        let path: PathBuf = value.file_path.into();
        let components: Vec<String> = path
            .components()
            .map(|comp| {
                comp.as_os_str()
                    .to_str()
                    .expect(
                        "This path must be valid unicode as it was built from a unicode string.",
                    )
                    .to_string()
            })
            .collect();
        Self {
            id: value.id,
            path_components: components,
            creation_time: value.creation_time,
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_file_path() {
        let id = 42;
        let creation_time = chrono::Utc::now().naive_utc();
        let file_db = GlobalDataFile {
            id,
            global_data_id: 3,
            file_path: "path/to/file.txt".to_string(),
            creation_time,
        };
        let file_exchange: GlobalDataFileDetails = file_db.into();
        assert_eq!(&file_exchange.id, &id);
        assert_eq!(&file_exchange.creation_time, &creation_time);
        assert_eq!(
            file_exchange.path_components,
            vec!["path".to_string(), "to".to_string(), "file.txt".to_string()]
        )
    }
}
