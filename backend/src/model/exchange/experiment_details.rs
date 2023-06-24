use chrono::NaiveDateTime;
use serde::{Deserialize, Serialize};

use crate::model::db::experiment::Experiment;

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ExperimentDetails {
    /// The database ID.
    pub id: i32,
    /// The display name.
    pub name: String,
    /// An optional comment on the experiment.
    pub comment: Option<String>,
    /// The ID of the selected pipeline if any.
    pub pipeline_id: Option<String>,
    /// An optional mail address to send info regarding the pipeline status to.
    pub mail: Option<String>,
    /// The time of creation.
    pub creation_time: NaiveDateTime,
}

impl From<Experiment> for ExperimentDetails {
    fn from(value: Experiment) -> Self {
        Self {
            id: value.id,
            name: value.experiment_name,
            comment: value.comment,
            creation_time: value.creation_time,
            pipeline_id: value.pipeline_id,
            mail: value.mail,
        }
    }
}
