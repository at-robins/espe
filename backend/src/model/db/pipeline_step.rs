use chrono::{DateTime, Utc};
use diesel::{Identifiable, Queryable};
use crate::schema::pipeline_step;

use super::pipeline::Pipeline;

#[derive(Identifiable, Queryable, Associations, PartialEq, Debug)]
#[belongs_to(Pipeline, foreign_key = "pipeline_id")]
#[table_name = "pipeline_step"]
/// A pipeline database entry.
pub struct PipelineStep {
    pub id: i64,
    pub pipeline_id: i64,
    pub execution_type: String,
    pub execution_configuration: String,
    pub ordering: i64,
    pub creation_time: DateTime<Utc>,
}