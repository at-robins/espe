use crate::{
    application::error::{SeqError, SeqErrorType},
    schema::experiment_execution::{self},
};
use chrono::{NaiveDateTime, Utc};
use diesel::{
    BoolExpressionMethods, ExpressionMethods, Identifiable, Insertable, QueryDsl, Queryable,
    RunQueryDsl, SqliteConnection,
};
use getset::Getters;

#[derive(Identifiable, Queryable, Insertable, PartialEq, Debug)]
#[diesel(table_name = experiment_execution)]
/// A queryable experiment execution database entry.
pub struct ExperimentExecution {
    pub id: i32,
    pub experiment_id: i32,
    pub pipeline_id: String,
    pub pipeline_step_id: String,
    pub execution_status: String,
    pub start_time: Option<NaiveDateTime>,
    pub end_time: Option<NaiveDateTime>,
    pub creation_time: NaiveDateTime,
}

const EXECUTION_STATUS_ABORTED: &str = "Aborted";
const EXECUTION_STATUS_FAILED: &str = "Failed";
const EXECUTION_STATUS_FINISHED: &str = "Finished";
const EXECUTION_STATUS_RUNNING: &str = "Running";
const EXECUTION_STATUS_WAITING: &str = "Waiting";

#[derive(Debug, PartialEq)]
/// The execution status of a pipeline step.
pub enum ExecutionStatus {
    Aborted,
    Failed,
    Finished,
    Running,
    Waiting,
}

impl From<&ExecutionStatus> for String {
    fn from(value: &ExecutionStatus) -> String {
        match value {
            ExecutionStatus::Aborted => EXECUTION_STATUS_ABORTED.to_string(),
            ExecutionStatus::Failed => EXECUTION_STATUS_FAILED.to_string(),
            ExecutionStatus::Finished => EXECUTION_STATUS_FINISHED.to_string(),
            ExecutionStatus::Running => EXECUTION_STATUS_RUNNING.to_string(),
            ExecutionStatus::Waiting => EXECUTION_STATUS_WAITING.to_string(),
        }
    }
}

impl From<ExecutionStatus> for String {
    fn from(value: ExecutionStatus) -> String {
        String::from(&value)
    }
}

impl TryFrom<&str> for ExecutionStatus {
    type Error = SeqError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            EXECUTION_STATUS_ABORTED => Ok(ExecutionStatus::Aborted),
            EXECUTION_STATUS_FAILED => Ok(ExecutionStatus::Failed),
            EXECUTION_STATUS_FINISHED => Ok(ExecutionStatus::Finished),
            EXECUTION_STATUS_RUNNING => Ok(ExecutionStatus::Running),
            EXECUTION_STATUS_WAITING => Ok(ExecutionStatus::Waiting),
            _ => Err(SeqError::new(
                "Conversion error",
                SeqErrorType::InternalServerError,
                format!("{} is not a valid execution status.", value),
                "Invalid execution status.",
            )),
        }
    }
}

impl TryFrom<String> for ExecutionStatus {
    type Error = SeqError;

    fn try_from(value: String) -> Result<Self, Self::Error> {
        Self::try_from(value.as_str())
    }
}

impl TryFrom<&String> for ExecutionStatus {
    type Error = SeqError;

    fn try_from(value: &String) -> Result<Self, Self::Error> {
        Self::try_from(value.as_str())
    }
}

impl ExperimentExecution {
    /// Tries to parse the entities execution status.
    pub fn execution_status(&self) -> Result<ExecutionStatus, SeqError> {
        (&self.execution_status).try_into()
    }

    /// Returns `true` if the experiment has execution entries.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the experiment ID
    /// * `connection` - the database connection
    pub fn has_experiment_execution_entries(
        experiment_id: i32,
        connection: &mut SqliteConnection,
    ) -> Result<bool, diesel::result::Error> {
        diesel::select(diesel::dsl::exists(
            crate::schema::experiment_execution::table
                .filter(crate::schema::experiment_execution::experiment_id.eq(experiment_id)),
        ))
        .get_result(connection)
    }

    /// Returns all entities that belong to the specified experiment.
    ///
    /// # Parameters
    ///
    /// * `experiemnt_id` - the ID of the experiment
    /// * `connection` - the database connection
    pub fn get_by_experiment(
        experiment_id: i32,
        connection: &mut SqliteConnection,
    ) -> Result<Vec<ExperimentExecution>, diesel::result::Error> {
        crate::schema::experiment_execution::table
            .filter(crate::schema::experiment_execution::experiment_id.eq(experiment_id))
            .load::<ExperimentExecution>(connection)
    }

    /// Returns all entities with the specified [`ExecutionStatus`].
    ///
    /// # Parameters
    ///
    /// * `status` - the [`ExecutionStatus`]
    /// * `connection` - the database connection
    pub fn get_by_status(
        status: ExecutionStatus,
        connection: &mut SqliteConnection,
    ) -> Result<Vec<ExperimentExecution>, diesel::result::Error> {
        let status: String = status.into();
        crate::schema::experiment_execution::table
            .filter(crate::schema::experiment_execution::execution_status.eq(status))
            .load::<ExperimentExecution>(connection)
    }

    /// Deletes all entities that belong to the specified experiment.
    ///
    /// # Parameters
    ///
    /// * `experiemnt_id` - the ID of the experiment
    /// * `connection` - the database connection
    pub fn delete_by_experiment(
        experiment_id: i32,
        connection: &mut SqliteConnection,
    ) -> Result<(), diesel::result::Error> {
        diesel::delete(crate::schema::experiment_execution::table)
            .filter(crate::schema::experiment_execution::experiment_id.eq(experiment_id))
            .execute(connection)?;
        Ok(())
    }

    /// Returns all entities.
    ///
    /// # Parameters
    ///
    /// * `connection` - the database connection
    pub fn get_all(
        connection: &mut SqliteConnection,
    ) -> Result<Vec<ExperimentExecution>, diesel::result::Error> {
        crate::schema::experiment_execution::table.load::<ExperimentExecution>(connection)
    }

    /// Sets the [`ExecutionStatus`] of all steps belonging to the specified experiment that are not finished.
    ///
    /// # Parameters
    ///
    /// * `experiemnt_id` - the ID of the experiment
    /// * `connection` - the database connection
    pub fn update_scheduled_status_by_experiment(
        experiment_id: i32,
        status: ExecutionStatus,
        connection: &mut SqliteConnection,
    ) -> Result<usize, diesel::result::Error> {
        let status: String = status.into();
        let status_not_finished: Vec<String> = vec![
            ExecutionStatus::Running.into(),
            ExecutionStatus::Waiting.into(),
        ];
        diesel::update(
            crate::schema::experiment_execution::table.filter(
                crate::schema::experiment_execution::experiment_id
                    .eq(experiment_id)
                    .and(
                        crate::schema::experiment_execution::execution_status
                            .eq_any(status_not_finished),
                    ),
            ),
        )
        .set(crate::schema::experiment_execution::execution_status.eq(status))
        .execute(connection)
    }
}

#[derive(Insertable, PartialEq, Debug, Getters)]
#[diesel(table_name = experiment_execution)]
/// A new experiment execution database record.
pub struct NewExperimentExecution {
    #[getset(get = "pub")]
    experiment_id: i32,
    #[getset(get = "pub")]
    pipeline_id: String,
    #[getset(get = "pub")]
    pipeline_step_id: String,
    #[getset(get = "pub")]
    execution_status: String,
    #[getset(get = "pub")]
    start_time: Option<NaiveDateTime>,
    #[getset(get = "pub")]
    end_time: Option<NaiveDateTime>,
    #[getset(get = "pub")]
    creation_time: NaiveDateTime,
}

impl NewExperimentExecution {
    /// Creates a new experiment execution record for insertion into the database.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the according experiment
    /// * `pipeline_id` - the ID of the pipeline to be executed
    /// * `pipeline_step_id` - the ID of the pipeline step to be executed
    pub fn new<
        ExperimentIdType: Into<i32>,
        PipelineIdIdType: Into<String>,
        PipelineStepIdIdType: Into<String>,
    >(
        experiment_id: ExperimentIdType,
        pipeline_id: PipelineIdIdType,
        pipeline_step_id: PipelineStepIdIdType,
    ) -> Self {
        Self {
            experiment_id: experiment_id.into(),
            pipeline_id: pipeline_id.into(),
            pipeline_step_id: pipeline_step_id.into(),
            execution_status: ExecutionStatus::Waiting.into(),
            start_time: None,
            end_time: None,
            creation_time: Utc::now().naive_utc(),
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::{model::db::experiment::NewExperiment, test_utility::TestContext};

    use super::*;

    #[test]
    fn test_execution_status_conversion() {
        assert_eq!(
            ExecutionStatus::Aborted,
            ExecutionStatus::try_from(EXECUTION_STATUS_ABORTED).unwrap()
        );
        assert_eq!(
            ExecutionStatus::Failed,
            ExecutionStatus::try_from(EXECUTION_STATUS_FAILED).unwrap()
        );
        assert_eq!(
            ExecutionStatus::Finished,
            ExecutionStatus::try_from(EXECUTION_STATUS_FINISHED).unwrap()
        );
        assert_eq!(
            ExecutionStatus::Running,
            ExecutionStatus::try_from(EXECUTION_STATUS_RUNNING).unwrap()
        );
        assert_eq!(
            ExecutionStatus::Waiting,
            ExecutionStatus::try_from(EXECUTION_STATUS_WAITING).unwrap()
        );
        assert!(ExecutionStatus::try_from("invalid execution step").is_err());
        assert_eq!(String::from(ExecutionStatus::Aborted).as_str(), EXECUTION_STATUS_ABORTED);
        assert_eq!(String::from(ExecutionStatus::Failed).as_str(), EXECUTION_STATUS_FAILED);
        assert_eq!(String::from(ExecutionStatus::Finished).as_str(), EXECUTION_STATUS_FINISHED);
        assert_eq!(String::from(ExecutionStatus::Running).as_str(), EXECUTION_STATUS_RUNNING);
        assert_eq!(String::from(ExecutionStatus::Waiting).as_str(), EXECUTION_STATUS_WAITING);
    }

    #[test]
    fn test_get_by_experiment() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        // Create a dummy experiments.
        let experiment_0 = NewExperiment::new("0".to_string());
        let experiment_1 = NewExperiment::new("1".to_string());
        let experiment_id_0: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_0)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        let experiment_id_1: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_1)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        // Dummy record setup.
        assert!(ExperimentExecution::get_by_experiment(experiment_id_0, &mut connection)
            .unwrap()
            .is_empty());
        let number_of_records = 42;
        let new_records_expected: Vec<ExperimentExecution> = (0..number_of_records)
            .map(|id| ExperimentExecution {
                id,
                experiment_id: experiment_id_0,
                pipeline_id: id.to_string(),
                pipeline_step_id: id.to_string(),
                execution_status: ExecutionStatus::Waiting.into(),
                start_time: None,
                end_time: None,
                creation_time: chrono::Utc::now().naive_local(),
            })
            .collect();
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_expected)
            .execute(&mut connection)
            .unwrap();
        let new_records_other: Vec<ExperimentExecution> = (0..number_of_records)
            .map(|id| ExperimentExecution {
                id: number_of_records + id,
                experiment_id: experiment_id_1,
                pipeline_id: id.to_string(),
                pipeline_step_id: id.to_string(),
                execution_status: ExecutionStatus::Waiting.into(),
                start_time: None,
                end_time: None,
                creation_time: chrono::Utc::now().naive_local(),
            })
            .collect();
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_other)
            .execute(&mut connection)
            .unwrap();
        assert_eq!(
            new_records_expected,
            ExperimentExecution::get_by_experiment(experiment_id_0, &mut connection).unwrap()
        );
    }

    #[test]
    fn test_get_all() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        // Create a dummy experiments.
        let experiment_0 = NewExperiment::new("0".to_string());
        let experiment_1 = NewExperiment::new("1".to_string());
        let experiment_id_0: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_0)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        let experiment_id_1: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_1)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        // Dummy record setup.
        assert!(ExperimentExecution::get_all(&mut connection)
            .unwrap()
            .is_empty());
        let number_of_records = 42;
        let new_records_expected: Vec<ExperimentExecution> = (0..number_of_records)
            .map(|id| ExperimentExecution {
                id,
                experiment_id: if id % 2 == 0 {
                    experiment_id_0
                } else {
                    experiment_id_1
                },
                pipeline_id: id.to_string(),
                pipeline_step_id: id.to_string(),
                execution_status: ExecutionStatus::Waiting.into(),
                start_time: None,
                end_time: None,
                creation_time: chrono::Utc::now().naive_local(),
            })
            .collect();
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_expected)
            .execute(&mut connection)
            .unwrap();
        assert_eq!(new_records_expected, ExperimentExecution::get_all(&mut connection).unwrap());
    }

    #[test]
    fn test_delete_by_experiment() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        // Create a dummy experiments.
        let experiment_0 = NewExperiment::new("0".to_string());
        let experiment_1 = NewExperiment::new("1".to_string());
        let experiment_id_0: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_0)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        let experiment_id_1: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_1)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        // Dummy record setup.
        assert!(ExperimentExecution::get_all(&mut connection)
            .unwrap()
            .is_empty());
        let number_of_records = 42;
        let new_records_all: Vec<ExperimentExecution> = (0..number_of_records)
            .map(|id| ExperimentExecution {
                id,
                experiment_id: if id % 2 == 0 {
                    experiment_id_0
                } else {
                    experiment_id_1
                },
                pipeline_id: id.to_string(),
                pipeline_step_id: id.to_string(),
                execution_status: ExecutionStatus::Waiting.into(),
                start_time: None,
                end_time: None,
                creation_time: chrono::Utc::now().naive_local(),
            })
            .collect();
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();
        let records_expected: Vec<ExperimentExecution> = new_records_all
            .into_iter()
            .filter(|record| record.experiment_id == experiment_id_1)
            .collect();
        ExperimentExecution::delete_by_experiment(experiment_id_0, &mut connection).unwrap();
        assert_eq!(records_expected, ExperimentExecution::get_all(&mut connection).unwrap());
    }

    #[test]
    fn test_has_experiment_execution_entries() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        // Create a dummy experiments.
        let experiment_0 = NewExperiment::new("0".to_string());
        let experiment_1 = NewExperiment::new("1".to_string());
        // experiment.
        let experiment_id_0: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_0)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        let experiment_id_1: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_1)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        // Dummy record setup.
        assert!(!ExperimentExecution::has_experiment_execution_entries(
            experiment_id_0,
            &mut connection
        )
        .unwrap());
        assert!(!ExperimentExecution::has_experiment_execution_entries(
            experiment_id_1,
            &mut connection
        )
        .unwrap());
        let number_of_records = 42;
        let new_records_all: Vec<ExperimentExecution> = (0..number_of_records)
            .map(|id| ExperimentExecution {
                id,
                experiment_id: experiment_id_0,
                pipeline_id: id.to_string(),
                pipeline_step_id: id.to_string(),
                execution_status: ExecutionStatus::Waiting.into(),
                start_time: None,
                end_time: None,
                creation_time: chrono::Utc::now().naive_local(),
            })
            .collect();
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_all)
            .execute(&mut connection)
            .unwrap();
        assert!(ExperimentExecution::has_experiment_execution_entries(
            experiment_id_0,
            &mut connection
        )
        .unwrap());
        assert!(!ExperimentExecution::has_experiment_execution_entries(
            experiment_id_1,
            &mut connection
        )
        .unwrap());
    }

    #[test]
    fn test_get_by_status() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        // Create a dummy experiments.
        let experiment_0 = NewExperiment::new("0".to_string());
        let experiment_id_0: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_0)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        let new_records_expected: Vec<ExperimentExecution> = vec![
            ExecutionStatus::Aborted,
            ExecutionStatus::Failed,
            ExecutionStatus::Finished,
            ExecutionStatus::Running,
            ExecutionStatus::Finished,
            ExecutionStatus::Waiting,
        ]
        .iter()
        .enumerate()
        .map(|(id, status)| ExperimentExecution {
            id: id as i32,
            experiment_id: experiment_id_0,
            pipeline_id: id.to_string(),
            pipeline_step_id: id.to_string(),
            execution_status: status.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        })
        .collect();
        assert!(ExperimentExecution::get_by_status(ExecutionStatus::Finished, &mut connection)
            .unwrap()
            .is_empty());
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_expected)
            .execute(&mut connection)
            .unwrap();

        let finished_records =
            ExperimentExecution::get_by_status(ExecutionStatus::Finished, &mut connection).unwrap();
            assert_eq!(finished_records.len(), 2);
            assert_eq!(finished_records[0], new_records_expected[2]);
            assert_eq!(finished_records[1], new_records_expected[4]);
    }

    #[test]
    fn test_update_scheduled_status_by_experiment() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        // Create a dummy experiments.
        let experiment_0 = NewExperiment::new("0".to_string());
        let experiment_id_0: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_0)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        let new_records_expected: Vec<ExperimentExecution> = vec![
            ExecutionStatus::Aborted,
            ExecutionStatus::Failed,
            ExecutionStatus::Finished,
            ExecutionStatus::Running,
            ExecutionStatus::Waiting,
        ]
        .iter()
        .enumerate()
        .map(|(id, status)| ExperimentExecution {
            id: id as i32,
            experiment_id: experiment_id_0,
            pipeline_id: id.to_string(),
            pipeline_step_id: id.to_string(),
            execution_status: status.into(),
            start_time: None,
            end_time: None,
            creation_time: chrono::Utc::now().naive_local(),
        })
        .collect();
        diesel::insert_into(crate::schema::experiment_execution::table)
            .values(&new_records_expected)
            .execute(&mut connection)
            .unwrap();
        let all_records_before_update =
            ExperimentExecution::get_by_experiment(experiment_id_0, &mut connection).unwrap();
        assert_eq!(
            all_records_before_update[0].execution_status().unwrap(),
            ExecutionStatus::Aborted
        );
        assert_eq!(
            all_records_before_update[1].execution_status().unwrap(),
            ExecutionStatus::Failed
        );
        assert_eq!(
            all_records_before_update[2].execution_status().unwrap(),
            ExecutionStatus::Finished
        );
        assert_eq!(
            all_records_before_update[3].execution_status().unwrap(),
            ExecutionStatus::Running
        );
        assert_eq!(
            all_records_before_update[4].execution_status().unwrap(),
            ExecutionStatus::Waiting
        );
        ExperimentExecution::update_scheduled_status_by_experiment(
            experiment_id_0,
            ExecutionStatus::Failed,
            &mut connection,
        )
        .unwrap();
        let all_records_after_update =
            ExperimentExecution::get_by_experiment(experiment_id_0, &mut connection).unwrap();
        assert_eq!(
            all_records_after_update[0].execution_status().unwrap(),
            ExecutionStatus::Aborted
        );
        assert_eq!(
            all_records_after_update[1].execution_status().unwrap(),
            ExecutionStatus::Failed
        );
        assert_eq!(
            all_records_after_update[2].execution_status().unwrap(),
            ExecutionStatus::Finished
        );
        assert_eq!(
            all_records_after_update[3].execution_status().unwrap(),
            ExecutionStatus::Failed
        );
        assert_eq!(
            all_records_after_update[4].execution_status().unwrap(),
            ExecutionStatus::Failed
        );
    }
}
