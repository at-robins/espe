use chrono::NaiveDateTime;
use serde::{Deserialize, Serialize};

use crate::model::db::global_data::GlobalData;

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
