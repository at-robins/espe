use std::{
    ffi::OsString,
    path::{Path, PathBuf},
    process::Command,
    sync::Arc,
};

use crate::{
    application::{
        config::{Configuration, PIPELINE_DEFINITION_FILE},
        error::{SeqError, SeqErrorType},
    },
    model::internal::pipeline_blueprint::{
        ContextualisedPipelineBlueprint, PipelineBlueprint, PipelineStepBlueprint,
    },
};

/// Returns all pipelines defined in the respective directory.
///
/// # Parameters
///
/// * `app_config` - the app [`Configuration`]
pub fn load_pipelines(
    app_config: Arc<Configuration>,
) -> Result<Vec<ContextualisedPipelineBlueprint>, SeqError> {
    let pipeline_path: PathBuf = app_config.pipeline_folder().into();
    let list = std::fs::read_dir(pipeline_path)?;
    let (dirs, errors): (Vec<std::fs::DirEntry>, Vec<std::io::Error>) =
        list.fold((Vec::new(), Vec::new()), |mut acc, entry| {
            if entry.is_ok() {
                acc.0.push(entry.unwrap());
            } else {
                acc.1.push(entry.unwrap_err());
            }
            acc
        });
    if errors.is_empty() {
        let pipeline_definition_files: Vec<PathBuf> = dirs
            .into_iter()
            .map(|dir| dir.path())
            .filter(|dir| dir.is_dir().into())
            .map(|mut dir| {
                dir.push(PIPELINE_DEFINITION_FILE);
                dir
            })
            .collect();
        let mut pipelines = Vec::with_capacity(pipeline_definition_files.len());
        for pipeline_def_file in pipeline_definition_files {
            let pipeline_def: PipelineBlueprint =
                serde_json::from_reader(std::fs::File::open(&pipeline_def_file)?)?;
            pipelines.push(ContextualisedPipelineBlueprint::new(
                pipeline_def,
                pipeline_def_file
                    .parent()
                    .expect("This unwrap of the parent path must work since we just appended the definition file."),
            ));
        }
        Ok(pipelines)
    } else {
        let combined_error = errors.into_iter().fold(String::new(), |mut acc, error| {
            acc.push_str(&error.to_string());
            acc.push('\n');
            acc
        });
        Err(SeqError::new(
            "std::io::Error",
            SeqErrorType::InternalServerError,
            combined_error,
            "An unforseen error occured during pipeline parsing. Please consult the logs.",
        ))
    }
}

// pub fn run_pipeline(pipeline: ContextualisedPipelineBlueprint) {
//     Command::new("docker").args(["build", &format!("-t {}", pipeline)]);
// }

/// Builds the specifc pipeline step container at the specified context.
///
/// # Parameters
///
/// * `step` - the [`PipelineStepBlueprint`] to build the container for
/// * `context` - the context directory contianing the pipeline
pub fn build_pipeline_step<P: AsRef<Path>>(
    step: &PipelineStepBlueprint,
    context: P,
) -> Result<std::process::Output, SeqError> {
    let mut pipeline_step_path = context.as_ref().to_path_buf();
    pipeline_step_path.push("container");
    pipeline_step_path.push(step.container());
    let build_arg: OsString = "build".into();
    let name_arg: OsString = format!("-t {}", step.id()).into();
    let output = Command::new("docker")
        .args([
            build_arg.as_os_str(),
            name_arg.as_os_str(),
            pipeline_step_path.as_os_str(),
        ])
        .output()?;
    Ok(output)
}

/// Runs the specifc pipeline step replacing all its previous output.
///
/// # Parameters
///
/// * `step` - the [`PipelineStepBlueprint`] to run
/// * `experiment_id` - the ID of the experiment
/// * `app_cofig` - the app [`Configuration`]
pub fn run_pipeline_step<P: AsRef<str>>(
    step: &PipelineStepBlueprint,
    experiment_id: P,
    app_config: Arc<Configuration>,
) -> Result<std::process::Output, SeqError> {
    // Create the output directory.
    let output_path = app_config.experiment_step_path(&experiment_id, step.id());
    // Clear the output folder if the step has been run before.
    if output_path.exists() {
        std::fs::remove_dir_all(&output_path)?;
    }
    // Then create the output directory.
    std::fs::create_dir_all(&output_path)?;
    // Set basic arguments.
    let mut arguments: Vec<OsString> = vec![
        "run".into(),
        format!("--name {}", step.id()).into(),
        "--rm".into(),
        pipeline_step_mount(output_path, "/output", false),
    ];
    // Set initial sample input mount.
    arguments.push(pipeline_step_mount(
        app_config.experiment_samples_path(&experiment_id),
        "/input/samples",
        true,
    ));
    // Set input mounts / dependencies.
    for dependency_id in step.dependencies() {
        arguments.push(pipeline_step_mount(
            app_config.experiment_step_path(&experiment_id, dependency_id),
            format!("/input/steps/{}", dependency_id),
            true,
        ));
    }
    // Set variables.
    todo!();

    // Set container to run.
    arguments.push(step.id().into());

    let output = Command::new("docker").args(arguments).output()?;
    Ok(output)
}

/// Creates an [`OsString`] for a container mount.
///
/// # Parameters
///
/// * `source` - the path to the local directory
/// * `target` - the path to the bound container directory
fn pipeline_step_mount<P: AsRef<Path>, Q: AsRef<Path>>(
    source: P,
    target: Q,
    readonly: bool,
) -> OsString {
    let mut mount_arg: OsString = "--mount type=bind,source=".into();
    mount_arg.push(source.as_ref());
    mount_arg.push(OsString::from(",target="));
    mount_arg.push(target.as_ref());
    if readonly {
        mount_arg.push(OsString::from(",readonly"));
    }
    mount_arg
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use serial_test::serial;

    use crate::{
        application::config::Configuration,
        test_utility::{TestContext, TEST_RESOURCES_PATH},
    };

    use super::*;

    #[test]
    fn test_pipeline_step_mount() {
        let input_arg = pipeline_step_mount(
            "/app/context/experiments/experiment_id/steps/step_id",
            "/input/steps/step_id",
            true,
        );
        let correct_input_arg: OsString =
            "--mount type=bind,source=/app/context/experiments/experiment_id/steps/step_id,target=/input/steps/step_id,readonly"
                .into();
        assert_eq!(input_arg, correct_input_arg);
        let output_arg = pipeline_step_mount(
            "/app/context/experiments/experiment_id/steps/step_id",
            "/output",
            false,
        );
        let correct_output_arg: OsString =
            "--mount type=bind,source=/app/context/experiments/experiment_id/steps/step_id,target=/output"
                .into();
        assert_eq!(output_arg, correct_output_arg);
    }

    // The tests need to be serial to prevent the same context UUID to be issued
    // to different tests at the same time.

    #[test]
    #[serial]
    fn test_load_pipelines_folder_not_exist() {
        let context = TestContext::new();
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let app_config: Arc<Configuration> = Arc::new((&context).into());
        let pipelines = load_pipelines(app_config);
        // The pipeline folder does not exist.
        assert!(pipelines.is_err());
    }

    #[test]
    #[serial]
    fn test_load_pipelines_empty() {
        let context = TestContext::new();
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let app_config: Arc<Configuration> = Arc::new((&context).into());
        std::fs::create_dir_all(context.pipeline_folder()).unwrap();
        let pipelines = load_pipelines(app_config).unwrap();
        assert!(pipelines.is_empty());
    }

    #[test]
    #[serial]
    fn test_load_pipelines() {
        let context = TestContext::new();
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let app_config: Arc<Configuration> = Arc::new(Configuration::new(
            context.database_url(),
            "info",
            "127.0.0.1",
            "8080",
            context.context_folder(),
            format!("{}/pipelines", TEST_RESOURCES_PATH),
        ));
        let pipelines = load_pipelines(app_config).unwrap();
        assert_eq!(pipelines.len(), 1);
        let test_pipeline = &pipelines[0];
        let pipeline_context: PathBuf =
            format!("{}/pipelines/testing_pipeline", TEST_RESOURCES_PATH,).into();
        assert_eq!(test_pipeline.context(), &pipeline_context);
        let pipeline: &PipelineBlueprint = test_pipeline.pipeline();
        assert_eq!(pipeline.id(), "testing_pipeline");
        assert_eq!(pipeline.name(), "Testing pipeline");
        assert_eq!(pipeline.description(), "This pipeline is for testing purposes.");
        assert_eq!(pipeline.steps().len(), 1);
        let step = &pipeline.steps()[0];
        assert_eq!(step.id(), "fastqc");
        assert_eq!(step.name(), "FastQC");
        assert_eq!(step.description(), "Performs a quality control.");
        assert_eq!(step.container(), "fastqc");
        assert!(step.dependencies().is_empty())
    }
}
