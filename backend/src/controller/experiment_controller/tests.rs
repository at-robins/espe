use super::*;

use crate::{
    application::config::Configuration,
    model::{
        db::{experiment::Experiment, pipeline_step_variable::NewPipelineStepVariable},
        internal::pipeline_blueprint::PipelineStepVariableCategory,
    },
    test_utility::{create_test_app, TestContext, TEST_RESOURCES_PATH},
};

use actix_web::{http::StatusCode, test};
use diesel::RunQueryDsl;

#[actix_web::test]
async fn test_create_experiment() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let exp_name = "Dummy experiment";
    let req = test::TestRequest::post()
        .uri("/api/experiments")
        .set_json(exp_name)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::CREATED);
    let id: i32 = test::read_body_json(resp).await;
    let test_data = Experiment::get(id, &mut connection).unwrap();
    assert_eq!(exp_name, &test_data.experiment_name);
    assert_eq!(&None, &test_data.comment);
    assert_eq!(id, test_data.id);
}

#[actix_web::test]
async fn test_create_experiment_title_too_long() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let exp_name: String = (0..513).fold(String::new(), |mut acc, _| {
        acc.push('a');
        acc
    });
    let req = test::TestRequest::post()
        .uri("/api/experiments")
        .set_json(exp_name)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
    assert!(Experiment::get_all(&mut connection).unwrap().is_empty());
}

#[actix_web::test]
async fn test_create_experiment_title_empty() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let req = test::TestRequest::post()
        .uri("/api/experiments")
        .set_json("")
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
    assert!(Experiment::get_all(&mut connection).unwrap().is_empty());
}

#[actix_web::test]
async fn test_delete_experiment() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
    let id = 42;
    let new_record = Experiment {
        id,
        experiment_name: "Dummy experiment".to_string(),
        comment: None,
        mail: None,
        pipeline_id: None,
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let experiment_path = app_config.experiment_path(id.to_string());
    let file_path = experiment_path.join("test_file.txt");
    std::fs::create_dir_all(&experiment_path).unwrap();
    std::fs::write(&file_path, "test content").unwrap();
    assert!(file_path.exists());
    assert!(!Experiment::get_all(&mut connection).unwrap().is_empty());
    let delete_req = test::TestRequest::delete()
        .uri(&format!("/api/experiments/{}", id))
        .to_request();
    let delete_resp = test::call_service(&app, delete_req).await;
    assert_eq!(delete_resp.status(), StatusCode::OK);
    assert!(Experiment::get_all(&mut connection).unwrap().is_empty());
    assert!(!experiment_path.exists());
}

#[actix_web::test]
async fn test_delete_experiment_non_existent() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    assert!(Experiment::get_all(&mut connection).unwrap().is_empty());
    let delete_req = test::TestRequest::delete()
        .uri("/api/experiments/42")
        .to_request();
    let delete_resp = test::call_service(&app, delete_req).await;
    assert_eq!(delete_resp.status(), StatusCode::NOT_FOUND);
    assert!(Experiment::get_all(&mut connection).unwrap().is_empty());
}

#[actix_web::test]
async fn test_get_experiment() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let new_record = Experiment {
        id,
        experiment_name: "Dummy record".to_string(),
        comment: Some("A comment".to_string()),
        mail: Some("a.b@c.de".to_string()),
        pipeline_id: Some("Dummy ID".to_string()),
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let req = test::TestRequest::get()
        .uri(&format!("/api/experiments/{}", id))
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    let fetched_data: ExperimentDetails = test::read_body_json(resp).await;
    assert_eq!(fetched_data.id, new_record.id);
    assert_eq!(fetched_data.name, new_record.experiment_name);
    assert_eq!(fetched_data.comment, new_record.comment);
    assert_eq!(fetched_data.creation_time, new_record.creation_time);
}

#[actix_web::test]
async fn test_get_experiment_non_existent() {
    let db_context = TestContext::new();
    let app = test::init_service(create_test_app(&db_context)).await;
    let req = test::TestRequest::get()
        .uri("/api/experiments/42")
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::NOT_FOUND);
}

#[actix_web::test]
async fn test_list_experiment() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let empty_req = test::TestRequest::get()
        .uri("/api/experiments")
        .to_request();
    let empty_resp = test::call_service(&app, empty_req).await;
    assert_eq!(empty_resp.status(), StatusCode::OK);
    let empty_fetched_data: Vec<ExperimentDetails> = test::read_body_json(empty_resp).await;
    assert!(empty_fetched_data.is_empty());
    let new_records: Vec<Experiment> = (0..42)
        .map(|id| Experiment {
            id,
            experiment_name: id.to_string(),
            comment: None,
            mail: None,
            pipeline_id: None,
            creation_time: chrono::Utc::now().naive_local(),
        })
        .collect();
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_records)
        .execute(&mut connection)
        .unwrap();
    let req = test::TestRequest::get()
        .uri("/api/experiments")
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    let fetched_data: Vec<ExperimentDetails> = test::read_body_json(resp).await;
    assert_eq!(new_records.len(), fetched_data.len());
    for i in 0..new_records.len() {
        assert_eq!(fetched_data[i].id, new_records[i].id);
        assert_eq!(fetched_data[i].name, new_records[i].experiment_name);
        assert_eq!(fetched_data[i].comment, new_records[i].comment);
        assert_eq!(fetched_data[i].mail, new_records[i].mail);
        assert_eq!(fetched_data[i].pipeline_id, new_records[i].pipeline_id);
        assert_eq!(fetched_data[i].creation_time, new_records[i].creation_time);
    }
}

#[actix_web::test]
async fn test_get_experiment_pipelines() {
    // Use a reference to the context, so the context is not dropped early
    // and messes up test context folder deletion.
    let mut db_context = TestContext::new();
    db_context.set_pipeline_folder(format!("{}/pipelines", TEST_RESOURCES_PATH));
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let new_experiment_record = Experiment {
        id,
        experiment_name: "Dummy record".to_string(),
        comment: Some("A comment".to_string()),
        mail: Some("a.b@c.de".to_string()),
        pipeline_id: Some("Dummy ID".to_string()),
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_experiment_record)
        .execute(&mut connection)
        .unwrap();
    let new_variable_records = vec![
        NewPipelineStepVariable::new(id, "testing_pipeline", "fastqc", "number", "123"),
        NewPipelineStepVariable::new(id, "testing_pipeline", "fastqc", "string", "abc"),
    ];
    diesel::insert_into(crate::schema::pipeline_step_variable::table)
        .values(&new_variable_records)
        .execute(&mut connection)
        .unwrap();

    let req = test::TestRequest::get()
        .uri(&format!("/api/experiments/{}/pipelines", id))
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    let fetched_data: Vec<ExperimentPipelineBlueprint> = test::read_body_json(resp).await;
    format!("{}/pipelines", TEST_RESOURCES_PATH);
    assert_eq!(fetched_data.len(), 1);
    let pipeline = &fetched_data[0];
    assert_eq!(pipeline.id(), "testing_pipeline");
    assert_eq!(pipeline.name(), "Testing pipeline");
    assert_eq!(pipeline.description(), "This pipeline is for testing purposes.");
    assert_eq!(pipeline.steps().len(), 1);
    let step = &pipeline.steps()[0];
    assert_eq!(step.id(), "fastqc");
    assert_eq!(step.name(), "FastQC");
    assert_eq!(step.description(), "Performs a quality control.");
    assert_eq!(step.container(), "fastqc");
    assert_eq!(step.dependencies(), &vec!["123", "456"]);
    assert_eq!(step.variables().len(), 5);
    assert_eq!(step.variables()[0].id(), "bool");
    assert_eq!(step.variables()[0].name(), "Boolean");
    assert_eq!(step.variables()[0].description(), "A boolean checkbox.");
    assert_eq!(step.variables()[0].category(), &PipelineStepVariableCategory::Boolean);
    assert_eq!(step.variables()[0].required(), &Some(true));
    assert_eq!(step.variables()[0].value(), &None);
    assert_eq!(step.variables()[1].id(), "global");
    assert_eq!(step.variables()[1].name(), "Global");
    assert_eq!(step.variables()[1].description(), "A global data reference.");
    assert_eq!(step.variables()[1].category(), &PipelineStepVariableCategory::Global);
    assert_eq!(step.variables()[1].required(), &Some(false));
    assert_eq!(step.variables()[1].value(), &None);
    assert_eq!(step.variables()[2].id(), "number");
    assert_eq!(step.variables()[2].name(), "Number");
    assert_eq!(step.variables()[2].description(), "A number field.");
    assert_eq!(step.variables()[2].category(), &PipelineStepVariableCategory::Number);
    assert_eq!(step.variables()[2].required(), &None);
    assert_eq!(step.variables()[2].value(), &Some("123".to_string()));
    assert_eq!(step.variables()[3].id(), "option");
    assert_eq!(step.variables()[3].name(), "Option");
    assert_eq!(step.variables()[3].description(), "An option dropdown.");
    if let PipelineStepVariableCategory::Option(options) = step.variables()[3].category() {
        assert_eq!(options.len(), 2);
        assert_eq!(options[0].name(), "Option 1");
        assert_eq!(options[0].value(), "option1");
        assert_eq!(options[1].name(), "Option 2");
        assert_eq!(options[1].value(), "option2");
    } else {
        panic!("Not an option variable!");
    }
    assert_eq!(step.variables()[3].required(), &None);
    assert_eq!(step.variables()[3].value(), &None);
    assert_eq!(step.variables()[4].id(), "string");
    assert_eq!(step.variables()[4].name(), "String");
    assert_eq!(step.variables()[4].description(), "A string text field.");
    assert_eq!(step.variables()[4].category(), &PipelineStepVariableCategory::String);
    assert_eq!(step.variables()[4].required(), &None);
    assert_eq!(step.variables()[4].value(), &Some("abc".to_string()));
}

#[actix_web::test]
async fn test_patch_experiment_comment() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let new_record = Experiment {
        id,
        experiment_name: "Dummy record".to_string(),
        comment: None,
        mail: None,
        pipeline_id: None,
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let new_comment = Some("A comment".to_string());
    let req = test::TestRequest::patch()
        .uri(&format!("/api/experiments/{}/comment", id))
        .set_json(&new_comment)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    assert_eq!(new_comment, Experiment::get(id, &mut connection).unwrap().comment);
    let clear_comment: Option<String> = None;
    let clear_req = test::TestRequest::patch()
        .uri(&format!("/api/experiments/{}/comment", id))
        .set_json(&clear_comment)
        .to_request();
    let clear_resp = test::call_service(&app, clear_req).await;
    assert_eq!(clear_resp.status(), StatusCode::OK);
    assert_eq!(clear_comment, Experiment::get(id, &mut connection).unwrap().comment);
}

#[actix_web::test]
async fn test_patch_global_comment_non_existent() {
    let db_context = TestContext::new();
    let app = test::init_service(create_test_app(&db_context)).await;
    let new_name = "A completely new name".to_string();
    let req = test::TestRequest::patch()
        .uri("/api/experiments/42/comment")
        .set_json(&new_name)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::NOT_FOUND);
}

#[actix_web::test]
async fn test_patch_experiment_name() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let new_record = Experiment {
        id,
        experiment_name: "Dummy record".to_string(),
        comment: None,
        mail: None,
        pipeline_id: None,
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let new_name = "A completely new name".to_string();
    let req = test::TestRequest::patch()
        .uri(&format!("/api/experiments/{}/name", id))
        .set_json(&new_name)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    assert_eq!(
        new_name,
        Experiment::get(id, &mut connection)
            .unwrap()
            .experiment_name
    );
}

#[actix_web::test]
async fn test_patch_experiment_name_non_existent() {
    let db_context = TestContext::new();
    let app = test::init_service(create_test_app(&db_context)).await;
    let new_name = "A completely new name".to_string();
    let req = test::TestRequest::patch()
        .uri("/api/experiments/42/name")
        .set_json(&new_name)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::NOT_FOUND);
}

#[actix_web::test]
async fn test_patch_experiment_name_empty() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let old_name = "Dummy record".to_string();
    let new_record = Experiment {
        id,
        experiment_name: old_name.clone(),
        comment: None,
        mail: None,
        pipeline_id: None,
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let req = test::TestRequest::patch()
        .uri(&format!("/api/experiments/{}/name", id))
        .set_json("")
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
    assert_eq!(
        old_name,
        Experiment::get(id, &mut connection)
            .unwrap()
            .experiment_name
    );
}

#[actix_web::test]
async fn test_patch_experiment_name_too_long() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let old_name = "Dummy record".to_string();
    let new_record = Experiment {
        id,
        experiment_name: old_name.clone(),
        comment: None,
        mail: None,
        pipeline_id: None,
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let new_name: String = (0..513).fold(String::new(), |mut acc, _| {
        acc.push('a');
        acc
    });
    let req = test::TestRequest::patch()
        .uri(&format!("/api/experiments/{}/name", id))
        .set_json(&new_name)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
    assert_eq!(
        old_name,
        Experiment::get(id, &mut connection)
            .unwrap()
            .experiment_name
    );
}

#[actix_web::test]
async fn test_patch_experiment_mail() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let new_record = Experiment {
        id,
        experiment_name: "Dummy record".to_string(),
        comment: None,
        mail: None,
        pipeline_id: None,
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let new_mail = "A.completely@new.mail".to_string();
    let req = test::TestRequest::patch()
        .uri(&format!("/api/experiments/{}/mail", id))
        .set_json(&new_mail)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    assert_eq!(new_mail, Experiment::get(id, &mut connection).unwrap().mail.unwrap());
}

#[actix_web::test]
async fn test_patch_experiment_mail_non_existent() {
    let db_context = TestContext::new();
    let app = test::init_service(create_test_app(&db_context)).await;
    let new_mail = "a.completely@new.mail".to_string();
    let req = test::TestRequest::patch()
        .uri("/api/experiments/42/mail")
        .set_json(&new_mail)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::NOT_FOUND);
}

#[actix_web::test]
async fn test_patch_experiment_mail_empty() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let old_mail = "an@old.mail".to_string();
    let new_record = Experiment {
        id,
        experiment_name: "Dummy record".to_string(),
        comment: None,
        mail: Some(old_mail.clone()),
        pipeline_id: None,
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let req = test::TestRequest::patch()
        .uri(&format!("/api/experiments/{}/mail", id))
        .set_json("")
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
    assert_eq!(old_mail, Experiment::get(id, &mut connection).unwrap().mail.unwrap());
}

#[actix_web::test]
async fn test_patch_experiment_mail_too_long() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let old_mail = "an@old.mail".to_string();
    let new_record = Experiment {
        id,
        experiment_name: "Dummy record".to_string(),
        comment: None,
        mail: Some(old_mail.clone()),
        pipeline_id: None,
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let new_mail_component: String = (0..=254).fold(String::new(), |mut acc, _| {
        acc.push('a');
        acc
    });
    let new_mail = format!("{}@{}.com", new_mail_component, new_mail_component);
    let req = test::TestRequest::patch()
        .uri(&format!("/api/experiments/{}/mail", id))
        .set_json(&new_mail)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
    assert_eq!(old_mail, Experiment::get(id, &mut connection).unwrap().mail.unwrap());
}

#[actix_web::test]
async fn test_patch_experiment_mail_invalid() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let old_mail = "an@old.mail".to_string();
    let new_record = Experiment {
        id,
        experiment_name: "Dummy record".to_string(),
        comment: None,
        mail: Some(old_mail.clone()),
        pipeline_id: None,
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let req = test::TestRequest::patch()
        .uri(&format!("/api/experiments/{}/mail", id))
        .set_json("a@b@c.de".to_string())
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
    assert_eq!(old_mail, Experiment::get(id, &mut connection).unwrap().mail.unwrap());
}
