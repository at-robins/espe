use crate::{
    model::exchange::experiment_upload::ExperimentUpload,
    schema::experiment::{self},
};
use chrono::{NaiveDateTime, Utc};
use diesel::{Identifiable, Insertable, Queryable};
use getset::{Getters, CopyGetters};

use super::pipeline::Pipeline;

#[derive(Identifiable, Queryable, Associations, Insertable, PartialEq, Debug)]
#[belongs_to(Pipeline)]
#[table_name = "experiment"]
/// A queryable experiment database entry.
pub struct Experiment {
    pub id: i32,
    pub experiment_name: String,
    pub mail: String,
    pub pipeline_id: i32,
    pub comment: String,
    pub creation_time: NaiveDateTime,
}

#[derive(Insertable, PartialEq, Debug, Getters, CopyGetters)]
#[table_name = "experiment"]
/// A new experiment database record.
pub struct NewExperiment {
    #[getset(get = "pub")]
    experiment_name: String,
    #[getset(get = "pub")]
    mail: String,
    #[getset(get_copy = "pub")]
    pipeline_id: i32,
    #[getset(get = "pub")]
    comment: String,
    #[getset(get = "pub")]
    creation_time: NaiveDateTime,
}

impl NewExperiment {
    /// Creates a new experiment record for insertion into the database.
    ///
    /// # Parameters
    ///
    /// * `name` - the experiment's name
    /// * `mail` - the E-mail address to notifiy on pipeline updates
    /// * `pipeline_id` - the referenced pipeline
    /// * `comment` - the comment describing the experiment
    pub fn new(name: String, mail: String, pipeline_id: i32, comment: String) -> Self {
        Self {
            experiment_name: name,
            mail,
            pipeline_id,
            comment,
            creation_time: Utc::now().naive_utc(),
        }
    }
}

impl From<ExperimentUpload> for NewExperiment {
    fn from(experiment_upload: ExperimentUpload) -> Self {
        NewExperiment::new(
            experiment_upload.name,
            experiment_upload.mail,
            experiment_upload.pipeline_id,
            experiment_upload.comment,
        )
    }
}
