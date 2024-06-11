use actix_web::web;

use crate::application::{
    config::Configuration,
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
/// * `app_cofig` - the app [`Configuration`]
pub fn delete_step_logs<P: AsRef<str>, S: AsRef<str>>(
    pipeline_id: P,
    step_id: S,
    experiment_id: i32,
    app_config: web::Data<Configuration>,
) -> Result<(), SeqError> {
    let deletion_errors: Vec<String> = app_config
        .experiment_log_paths_all(experiment_id.to_string(), &pipeline_id, &step_id)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_() {}
}
