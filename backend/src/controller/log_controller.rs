use crate::{
    application::{
        config::{Configuration, LogOutputType, LogProcessType},
        database::DatabaseManager,
        error::SeqError,
    },
    model::{
        db::experiment::Experiment,
        exchange::experiment_step_logs::{
            ExperimentStepLog, ExperimentStepLogRequest, ExperimentStepLogs,
        },
    },
};
use actix_web::web;

pub async fn get_experiment_step_logs(
    database_manager: web::Data<DatabaseManager>,
    app_config: web::Data<Configuration>,
    experiment_id: web::Path<i32>,
    info: web::Json<ExperimentStepLogRequest>,
) -> Result<web::Json<ExperimentStepLogs>, SeqError> {
    let experiment_id: i32 = experiment_id.into_inner();
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(experiment_id, &mut connection)?;

    let info: ExperimentStepLogRequest = info.into_inner();
    let log_reader = LogFileReader {
        config: app_config,
        experiment_id,
        pipeline_id: info.pipeline_id,
        step_id: info.step_id,
    };

    let build_logs = ExperimentStepLog {
        stdout: log_reader.get(LogProcessType::Build, LogOutputType::StdOut)?,
        stderr: log_reader.get(LogProcessType::Build, LogOutputType::StdErr)?,
        exit_code: log_reader.get(LogProcessType::Build, LogOutputType::ExitCode)?,
    };
    let run_logs = ExperimentStepLog {
        stdout: log_reader.get(LogProcessType::Run, LogOutputType::StdOut)?,
        stderr: log_reader.get(LogProcessType::Run, LogOutputType::StdErr)?,
        exit_code: log_reader.get(LogProcessType::Run, LogOutputType::ExitCode)?,
    };
    let logs = ExperimentStepLogs {
        build: build_logs,
        run: run_logs,
    };
    Ok(web::Json(logs))
}

/// A reader for log files.
struct LogFileReader {
    pub config: web::Data<Configuration>,
    pub experiment_id: i32,
    pub pipeline_id: String,
    pub step_id: String,
}

impl LogFileReader {
    /// Reads the respective log file to a [`String`] if it exists.
    ///
    /// # Parameters
    ///
    /// * `process_type` - the process type of the log file
    /// * `output_type` - the output type of the log file
    pub fn get(
        &self,
        process_type: LogProcessType,
        output_type: LogOutputType,
    ) -> Result<Option<String>, SeqError> {
        let path = self.config.experiment_log_path(
            self.experiment_id.to_string(),
            &self.pipeline_id,
            &self.step_id,
            process_type,
            output_type,
        );
        Ok(if !path.exists() {
            None
        } else {
            Some(std::fs::read_to_string(path)?)
        })
    }
}
