use std::{
    fs::File,
    io::Read,
    os::unix::fs::{FileExt, MetadataExt},
};

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

/// The maximum size of a single log file that is transmitted.
const MAX_TRANSMISSION_LOG_SIZE: usize = 512 * 1024;
/// The size of log read buffers if the transmission limit is exceeded.
const TRANSMISSION_BUFFER_SIZE: usize = MAX_TRANSMISSION_LOG_SIZE / 2;

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
            let mut file = File::open(&path)?;
            let file_size = file.metadata()?.size();
            if file_size > MAX_TRANSMISSION_LOG_SIZE as u64 {
                // Trims large log files so they do not block the frontend.
                log::warn!(
                    "The log file {} exceeds the limit of {} bytes and is truncated",
                    path.display(),
                    MAX_TRANSMISSION_LOG_SIZE
                );
                let mut start_buffer = [0u8; TRANSMISSION_BUFFER_SIZE];
                let mut end_buffer = [0u8; TRANSMISSION_BUFFER_SIZE];
                file.read_exact(&mut start_buffer)?;
                file.read_exact_at(&mut end_buffer, file_size - TRANSMISSION_BUFFER_SIZE as u64)?;
                let trimmed_content = format!(
                    "{}\n\n[ WARNING: The log content exceeded the size limit and has been trimmed. ]\n\n{}", 
                    String::from_utf8_lossy(&start_buffer),
                    String::from_utf8_lossy(&end_buffer)
                );
                Some(trimmed_content)
            } else {
                // Sends reasonably sized log files completely and handles potential errors.
                let file_content = std::fs::read(&path)?;
                match String::from_utf8(file_content) {
                    Ok(value) => Some(value),
                    Err(err) => {
                        // Log the error but still return the log file since most of the
                        // content is probably still readable.
                        // Invalid UTF-8 can for example be produced by container build logs
                        // being trimmed externally.
                        log::error!(
                            "The log file {} contains invalid UTF-8: {}\n{}",
                            path.display(),
                            err,
                            err.utf8_error()
                        );
                        let log_content_with_error: String = format!(
                            "[ ERROR: The log file containes invalid UTF-8. \
                            Please check the server logs for further details. ]\n\n{}",
                            String::from_utf8_lossy(err.as_bytes())
                        );
                        Some(log_content_with_error)
                    },
                }
            }
        })
    }
}
