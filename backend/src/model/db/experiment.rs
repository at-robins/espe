use crate::{
    model::exchange::experiment_upload::ExperimentUpload,
    schema::experiment::{self},
};
use chrono::{NaiveDateTime, Utc};
use diesel::{Identifiable, Insertable, Queryable};
use getset::{CopyGetters, Getters};

use super::pipeline::Pipeline;

#[derive(Identifiable, Queryable, Associations, Insertable, PartialEq, Debug)]
#[diesel(belongs_to(Pipeline))]
#[diesel(table_name = experiment)]
/// A queryable experiment database entry.
pub struct Experiment {
    pub id: i32,
    pub experiment_name: String,
    pub mail: Option<String>,
    pub pipeline_id: i32,
    pub comment: Option<String>,
    pub creation_time: NaiveDateTime,
}

#[derive(Insertable, PartialEq, Debug, Getters, CopyGetters)]
#[diesel(table_name = experiment)]
/// A new experiment database record.
pub struct NewExperiment {
    #[getset(get = "pub")]
    experiment_name: String,
    #[getset(get = "pub")]
    mail: Option<String>,
    #[getset(get_copy = "pub")]
    pipeline_id: i32,
    #[getset(get = "pub")]
    comment: Option<String>,
    #[getset(get = "pub")]
    creation_time: NaiveDateTime,
}

impl NewExperiment {
    /// Creates a new experiment record for insertion into the database.
    ///
    /// # Parameters
    ///
    /// * `name` - the experiment's name
    /// * `mail` - the optional E-mail address to notifiy on pipeline updates
    /// * `pipeline_id` - the referenced pipeline
    /// * `comment` - the optional comment describing the experiment
    pub fn new(
        name: String,
        mail: Option<String>,
        pipeline_id: i32,
        comment: Option<String>,
    ) -> Self {
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
