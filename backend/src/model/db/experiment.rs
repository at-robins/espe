use chrono::{DateTime, Utc};
use diesel::{Identifiable, Queryable};
use crate::schema::experiment;

use super::pipeline::Pipeline;

#[derive(Identifiable, Queryable, Associations, PartialEq, Debug)]
#[belongs_to(Pipeline)]
#[table_name = "experiment"]
/// A pipeline database entry.
pub struct Experiment {
    pub id: i64,
    pub experiment_name: String,
    pub mail: String,
    pub pipeline_id: i64,
    pub comment: String,
    pub creation_time: DateTime<Utc>,
}