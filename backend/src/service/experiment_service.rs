use std::fs::File;

use actix_web::web;

use crate::application::{
    config::{Configuration, LogOutputType, LogProcessType},
    error::{SeqError, SeqErrorType},
};

/// Deletes the specific pipeline step output of the specified experiment.
///
/// # Parameters
///
/// * `step_id` - the ID of the [`PipelineStepBlueprint`](crate::model::internal::pipeline_blueprint::PipelineStepBlueprint)
/// * `experiment_id` - the ID of the experiment
/// * `app_cofig` - the app [`Configuration`]
pub fn delete_step_output<S: AsRef<str>>(
    step_id: S,
    experiment_id: i32,
    app_config: web::Data<Configuration>,
) -> Result<(), SeqError> {
    let experiment_id = experiment_id.to_string();
    let output_path = app_config.experiment_step_path(&experiment_id, step_id);
    if output_path.exists() {
        std::fs::remove_dir_all(&output_path)?;
    }
    Ok(())
}

/// Deletes the specific pipeline step logs of the specified experiment.
///
/// # Parameters
///
/// * `pipeline_id` - the ID of the [`PipelineBlueprint`](crate::model::internal::pipeline_blueprint::PipelineBlueprint)
/// * `step_id` - the ID of the [`PipelineStepBlueprint`](crate::model::internal::pipeline_blueprint::PipelineStepBlueprint)
/// * `experiment_id` - the ID of the experiment
/// * `log_types` - the [`LogProcessType`]s to delete
/// * `app_cofig` - the app [`Configuration`]
pub fn delete_step_logs<P: AsRef<str>, S: AsRef<str>>(
    pipeline_id: P,
    step_id: S,
    experiment_id: i32,
    log_types: &[LogProcessType],
    app_config: web::Data<Configuration>,
) -> Result<(), SeqError> {
    let deletion_errors: Vec<String> = app_config
        .experiment_log_paths(experiment_id.to_string(), &pipeline_id, &step_id, log_types)
        .iter()
        .filter(|log_path| log_path.exists())
        .map(|log_path| {
            std::fs::remove_file(log_path).map_err(|error| {
                format!("Error while deleting log file {}: {}", log_path.display(), error)
            })
        })
        .filter(Result::is_err)
        .map(|deletion_error| deletion_error.unwrap_err())
        .collect();
    if deletion_errors.is_empty() {
        Ok(())
    } else {
        Err(SeqError::new(
            "Log deletion error",
            SeqErrorType::InternalServerError,
            deletion_errors
                .into_iter()
                .fold(String::new(), |acc, error_message| format!("{}{}\n\n", acc, error_message)),
            "Error during log deletion.",
        ))
    }
}

/// Creates context folder structure for running the
/// specific pipeline step of the specified experiment.
///
/// # Parameters
///
/// * `step_id` - the ID of the [`PipelineStepBlueprint`](crate::model::internal::pipeline_blueprint::PipelineStepBlueprint)
/// * `experiment_id` - the ID of the experiment
/// * `app_cofig` - the app [`Configuration`]
pub fn create_run_context<S: AsRef<str>>(
    step_id: S,
    experiment_id: i32,
    app_config: web::Data<Configuration>,
) -> Result<(), SeqError> {
    let experiment_id = experiment_id.to_string();
    // Creates input directory in case an empty pipeline is run.
    std::fs::create_dir_all(app_config.experiment_input_path(&experiment_id))?;
    // (Re-)Creates the output directory.
    let output_path = app_config.experiment_step_path(&experiment_id, &step_id);
    // Then create the output directory.
    std::fs::create_dir_all(&output_path)?;
    // Creates the log directory.
    let logs_path = app_config.experiment_logs_path(&experiment_id);
    std::fs::create_dir_all(&logs_path)?;
    Ok(())
}

/// Creates context folder structure for building the
/// specific pipeline step of the specified experiment.
///
/// # Parameters
///
/// * `experiment_id` - the ID of the experiment
/// * `app_cofig` - the app [`Configuration`]
pub fn create_build_context(
    experiment_id: i32,
    app_config: web::Data<Configuration>,
) -> Result<(), SeqError> {
    // Creates the log directory.
    let logs_path = app_config.experiment_logs_path(experiment_id.to_string());
    std::fs::create_dir_all(&logs_path)?;
    Ok(())
}

/// Deletes old context files and sets up the needed folder structure
/// for running the specific pipeline step of the specified experiment.
///
/// # Parameters
///
/// * `pipeline_id` - the ID of the [`PipelineBlueprint`](crate::model::internal::pipeline_blueprint::PipelineBlueprint)
/// * `step_id` - the ID of the [`PipelineStepBlueprint`](crate::model::internal::pipeline_blueprint::PipelineStepBlueprint)
/// * `experiment_id` - the ID of the experiment
/// * `app_cofig` - the app [`Configuration`]
pub fn prepare_context_for_run<P: AsRef<str>, S: AsRef<str>>(
    pipeline_id: P,
    step_id: S,
    experiment_id: i32,
    app_config: web::Data<Configuration>,
) -> Result<(), SeqError> {
    delete_step_logs(
        &pipeline_id,
        &step_id,
        experiment_id,
        &[LogProcessType::Run],
        web::Data::clone(&app_config),
    )?;
    delete_step_output(&step_id, experiment_id, web::Data::clone(&app_config))?;
    create_run_context(&step_id, experiment_id, app_config)?;
    Ok(())
}

/// Deletes old context files and sets up the needed folder structure
/// for building the specific pipeline step of the specified experiment.
///
/// # Parameters
///
/// * `pipeline_id` - the ID of the [`PipelineBlueprint`](crate::model::internal::pipeline_blueprint::PipelineBlueprint)
/// * `step_id` - the ID of the [`PipelineStepBlueprint`](crate::model::internal::pipeline_blueprint::PipelineStepBlueprint)
/// * `experiment_id` - the ID of the experiment
/// * `app_cofig` - the app [`Configuration`]
pub fn prepare_context_for_build<P: AsRef<str>, S: AsRef<str>>(
    pipeline_id: P,
    step_id: S,
    experiment_id: i32,
    app_config: web::Data<Configuration>,
) -> Result<(), SeqError> {
    delete_step_logs(
        &pipeline_id,
        &step_id,
        experiment_id,
        &[LogProcessType::Build],
        web::Data::clone(&app_config),
    )?;
    create_build_context(experiment_id, app_config)?;
    Ok(())
}

/// Opens the specified log file.
///
/// # Parameters
///
/// * `pipeline_id` - the ID of the [`PipelineBlueprint`](crate::model::internal::pipeline_blueprint::PipelineBlueprint)
/// * `step_id` - the ID of the [`PipelineStepBlueprint`](crate::model::internal::pipeline_blueprint::PipelineStepBlueprint)
/// * `experiment_id` - the ID of the experiment
/// * `log_process_type` - the [`LogProcessType`]
/// * `log_output_type` - the [`LogOutputType`]
/// * `app_cofig` - the app [`Configuration`]
pub fn open_step_log<P: AsRef<str>, S: AsRef<str>>(
    pipeline_id: P,
    step_id: S,
    experiment_id: i32,
    log_process_type: LogProcessType,
    log_output_type: LogOutputType,
    app_config: web::Data<Configuration>,
) -> Result<File, SeqError> {
    // Open stdout log file.
    let log_path = app_config.experiment_log_path(
        experiment_id.to_string(),
        pipeline_id,
        step_id,
        log_process_type,
        log_output_type,
    );
    let log_file = std::fs::OpenOptions::new()
        .create(true)
        .write(true)
        .append(false)
        .truncate(true)
        .open(log_path)?;

    Ok(log_file)
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::test_utility::TestContext;

    use super::*;

    /// Creates a dummy file at the specified path with dummy content.
    ///
    /// # Parameters
    ///
    /// * `file_path` - the [`Path`] to the file
    fn create_dummy_file<P: AsRef<Path>>(file_path: P) {
        if let Some(parent_path) = file_path.as_ref().parent() {
            std::fs::create_dir_all(parent_path).unwrap();
        }
        std::fs::write(file_path, "dummy content".as_bytes()).unwrap();
    }

    #[test]
    fn test_delete_step_output() {
        let context = TestContext::new();
        let app_config: Configuration = (&context).into();
        let experiment_id = 42;
        let step_id = "123456789";
        let output_path = app_config.experiment_step_path(experiment_id.to_string(), step_id);
        let test_output_path = output_path.join("test_dir/test.txt");
        create_dummy_file(&test_output_path);
        assert!(output_path.exists());
        assert!(test_output_path.exists());
        delete_step_output(step_id, experiment_id, web::Data::new(app_config)).unwrap();
        assert!(!output_path.exists());
        assert!(!test_output_path.exists());
    }

    #[test]
    fn test_delete_step_logs() {
        let context = TestContext::new();
        let app_config: Configuration = (&context).into();
        let experiment_id = 42;
        let step_id = "123456789";
        let pipeline_id = "test_pipeline";
        let build_log_paths = app_config.experiment_log_paths(
            experiment_id.to_string(),
            &pipeline_id,
            &step_id,
            &[LogProcessType::Build],
        );
        let run_log_paths = app_config.experiment_log_paths(
            experiment_id.to_string(),
            &pipeline_id,
            &step_id,
            &[LogProcessType::Run],
        );
        // Creates run and build log files.
        for log_path in build_log_paths.iter().chain(run_log_paths.iter()) {
            create_dummy_file(log_path);
        }
        // Deletes the run log files.
        delete_step_logs(
            pipeline_id,
            step_id,
            experiment_id,
            &[LogProcessType::Run],
            web::Data::new(app_config),
        )
        .unwrap();
        // Run logs should be deleted, while build logs should still exist.
        for build_log_path in build_log_paths {
            assert!(build_log_path.exists());
        }
        for run_log_path in run_log_paths {
            assert!(!run_log_path.exists());
        }
    }

    #[test]
    fn test_create_run_context() {
        let context = TestContext::new();
        let app_config: Configuration = (&context).into();
        let experiment_id = 42;
        let step_id = "123456789";

        let input_path = app_config.experiment_input_path(experiment_id.to_string());
        let output_path = app_config.experiment_step_path(experiment_id.to_string(), step_id);
        let logs_path = app_config.experiment_logs_path(experiment_id.to_string());

        assert!(!input_path.exists());
        assert!(!output_path.exists());
        assert!(!logs_path.exists());
        create_run_context(step_id, experiment_id, web::Data::new(app_config)).unwrap();
        assert!(input_path.exists());
        assert!(output_path.exists());
        assert!(logs_path.exists());
    }

    #[test]
    fn test_create_build_context() {
        let context = TestContext::new();
        let app_config: Configuration = (&context).into();
        let experiment_id = 42;
        let logs_path = app_config.experiment_logs_path(experiment_id.to_string());

        assert!(!logs_path.exists());
        create_build_context(experiment_id, web::Data::new(app_config)).unwrap();
        assert!(logs_path.exists());
    }

    #[test]
    fn test_prepare_context_for_run_empty() {
        let context = TestContext::new();
        let app_config: Configuration = (&context).into();
        let experiment_id = 42;
        let step_id = "123456789";
        let pipeline_id = "test_pipeline";

        let input_path = app_config.experiment_input_path(experiment_id.to_string());
        let output_path = app_config.experiment_step_path(experiment_id.to_string(), step_id);
        let logs_path = app_config.experiment_logs_path(experiment_id.to_string());

        assert!(!input_path.exists());
        assert!(!output_path.exists());
        assert!(!logs_path.exists());

        prepare_context_for_run(pipeline_id, step_id, experiment_id, web::Data::new(app_config))
            .unwrap();

        assert!(input_path.exists());
        assert!(output_path.exists());
        assert!(logs_path.exists());
    }

    #[test]
    fn test_prepare_context_for_run_filled() {
        let context = TestContext::new();
        let app_config: Configuration = (&context).into();
        let experiment_id = 42;
        let step_id = "123456789";
        let pipeline_id = "test_pipeline";

        let input_path = app_config.experiment_input_path(experiment_id.to_string());
        let output_path = app_config.experiment_step_path(experiment_id.to_string(), step_id);
        let logs_path = app_config.experiment_logs_path(experiment_id.to_string());
        let build_log_paths = app_config.experiment_log_paths(
            experiment_id.to_string(),
            &pipeline_id,
            &step_id,
            &[LogProcessType::Build],
        );
        let run_log_paths = app_config.experiment_log_paths(
            experiment_id.to_string(),
            &pipeline_id,
            &step_id,
            &[LogProcessType::Run],
        );

        // Dummy files.
        let step_file_path = output_path.join("dummy.txt");
        create_dummy_file(&step_file_path);
        // Creates run and build log files.
        for log_path in build_log_paths.iter().chain(run_log_paths.iter()) {
            create_dummy_file(log_path);
        }

        prepare_context_for_run(pipeline_id, step_id, experiment_id, web::Data::new(app_config))
            .unwrap();

        // All required directories should exist...
        assert!(input_path.exists());
        assert!(output_path.exists());
        assert!(logs_path.exists());
        // ... but old files should be deleted.
        assert!(!step_file_path.exists());
        // Run logs should be deleted, while build logs should still exist.
        for build_log_path in build_log_paths {
            assert!(build_log_path.exists());
        }
        for run_log_path in run_log_paths {
            assert!(!run_log_path.exists());
        }
    }

    #[test]
    fn test_prepare_context_for_build_empty() {
        let context = TestContext::new();
        let app_config: Configuration = (&context).into();
        let experiment_id = 42;
        let step_id = "123456789";
        let pipeline_id = "test_pipeline";

        let logs_path = app_config.experiment_logs_path(experiment_id.to_string());

        assert!(!logs_path.exists());

        prepare_context_for_build(pipeline_id, step_id, experiment_id, web::Data::new(app_config))
            .unwrap();

        assert!(logs_path.exists());
    }

    #[test]
    fn test_prepare_context_for_build_filled() {
        let context = TestContext::new();
        let app_config: Configuration = (&context).into();
        let experiment_id = 42;
        let step_id = "123456789";
        let pipeline_id = "test_pipeline";

        let input_path = app_config.experiment_input_path(experiment_id.to_string());
        let output_path = app_config.experiment_step_path(experiment_id.to_string(), step_id);
        let logs_path = app_config.experiment_logs_path(experiment_id.to_string());
        let build_log_paths = app_config.experiment_log_paths(
            experiment_id.to_string(),
            &pipeline_id,
            &step_id,
            &[LogProcessType::Build],
        );
        let run_log_paths = app_config.experiment_log_paths(
            experiment_id.to_string(),
            &pipeline_id,
            &step_id,
            &[LogProcessType::Run],
        );

        // Dummy files.
        let step_file_path = output_path.join("dummy.txt");
        create_dummy_file(&step_file_path);
        // Creates run and build log files.
        for log_path in build_log_paths.iter().chain(run_log_paths.iter()) {
            create_dummy_file(log_path);
        }

        prepare_context_for_build(pipeline_id, step_id, experiment_id, web::Data::new(app_config))
            .unwrap();

        // Only the log path should be created...
        assert!(!input_path.exists());
        assert!(logs_path.exists());
        // ... and old output files should still exist.
        assert!(output_path.exists());
        assert!(step_file_path.exists());
        // Build logs should be deleted, while run logs should still exist.
        for build_log_path in build_log_paths {
            assert!(!build_log_path.exists());
        }
        for run_log_path in run_log_paths {
            assert!(run_log_path.exists());
        }
    }
}
