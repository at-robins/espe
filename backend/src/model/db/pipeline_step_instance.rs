use crate::schema::pipeline_step_instance;
use chrono::NaiveDateTime;
use diesel::{Identifiable, Queryable};

use super::experiment::Experiment;
use super::pipeline_step::PipelineStep;

#[derive(Identifiable, Queryable, Associations, PartialEq, Debug)]
#[belongs_to(PipelineStep, foreign_key = "pipeline_step_id")]
#[belongs_to(Experiment, foreign_key = "experiment_id")]
#[table_name = "pipeline_step_instance"]
/// A pipeline step instance database entry.
pub struct PipelineStepInstance {
    pub id: i32,
    pub pipeline_step_id: i32,
    pub experiment_id: i32,
    pub pipeline_step_status: String,
    pub creation_time: NaiveDateTime,
}
