use crate::test_utility::{create_test_app, TestContext};

use super::*;
use actix_web::{
    http::{header::ContentType, StatusCode},
    test,
};
use diesel::RunQueryDsl;
use mime;

#[actix_web::test]
async fn test_get_experiment_files() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
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
    let experiment_path = app_config.experiment_path(id.to_string());
    let folders = vec![
        experiment_path.join("1/11/111"),
        experiment_path.join("1/11/112"),
        experiment_path.join("1/11/113"),
        experiment_path.join("2/21"),
        experiment_path.join("2/22"),
        experiment_path.join("3"),
    ];
    for folder in folders {
        std::fs::create_dir_all(folder).unwrap();
    }
    let test_file_path = experiment_path.join("1/11/112/test_file.txt");
    std::fs::write(test_file_path, "test_content").unwrap();

    let req = test::TestRequest::get()
        .uri(&format!("/api/files/experiments/{}", id))
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    let fetched_data: Vec<FileDetails> = test::read_body_json(resp).await;
    let expected_data: Vec<FileDetails> = vec![
        FileDetails {
            path_components: vec!["1".to_string()],
            is_file: false,
        },
        FileDetails {
            path_components: vec!["1".to_string(), "11".to_string()],
            is_file: false,
        },
        FileDetails {
            path_components: vec!["1".to_string(), "11".to_string(), "111".to_string()],
            is_file: false,
        },
        FileDetails {
            path_components: vec!["1".to_string(), "11".to_string(), "112".to_string()],
            is_file: false,
        },
        FileDetails {
            path_components: vec![
                "1".to_string(),
                "11".to_string(),
                "112".to_string(),
                "test_file.txt".to_string(),
            ],
            is_file: true,
        },
        FileDetails {
            path_components: vec!["1".to_string(), "11".to_string(), "113".to_string()],
            is_file: false,
        },
        FileDetails {
            path_components: vec!["2".to_string()],
            is_file: false,
        },
        FileDetails {
            path_components: vec!["2".to_string(), "21".to_string()],
            is_file: false,
        },
        FileDetails {
            path_components: vec!["2".to_string(), "22".to_string()],
            is_file: false,
        },
        FileDetails {
            path_components: vec!["3".to_string()],
            is_file: false,
        },
    ];
    assert_eq!(fetched_data, expected_data);
}

#[actix_web::test]
async fn test_get_experiment_files_non_existent() {
    let context = TestContext::new();
    let app = test::init_service(create_test_app(&context)).await;
    let req = test::TestRequest::get()
        .uri("/api/files/experiments/42")
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::NOT_FOUND);
}

#[actix_web::test]
async fn test_delete_experiment_files_by_path() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
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
    let experiment_path = app_config.experiment_path(id.to_string());
    let folders = vec![
        experiment_path.join("1/11/111"),
        experiment_path.join("1/11/112"),
        experiment_path.join("1/11/113"),
        experiment_path.join("2/21"),
        experiment_path.join("2/22"),
        experiment_path.join("3"),
    ];
    for folder in folders {
        std::fs::create_dir_all(folder).unwrap();
    }
    std::fs::write(experiment_path.join("1/11/112/test_file_1.txt"), "test_content").unwrap();
    std::fs::write(experiment_path.join("3/test_file_2.txt"), "test_content").unwrap();
    // Assert that all files and folder exist.
    assert!(experiment_path.join("1/11/111").exists());
    assert!(experiment_path.join("1/11/112").exists());
    assert!(experiment_path.join("1/11/112/test_file_1.txt").exists());
    assert!(experiment_path.join("1/11/113").exists());
    assert!(experiment_path.join("2/21").exists());
    assert!(experiment_path.join("2/22").exists());
    assert!(experiment_path.join("3").exists());
    assert!(experiment_path.join("3/test_file_2.txt").exists());
    // Delete a folder without content.
    let terminal_folder_path = FilePath {
        path_components: vec!["2".to_string(), "21".to_string()],
    };
    let req = test::TestRequest::delete()
        .uri(&format!("/api/files/experiments/{}", id))
        .set_json(terminal_folder_path)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    assert!(experiment_path.join("1/11/111").exists());
    assert!(experiment_path.join("1/11/112").exists());
    assert!(experiment_path.join("1/11/112/test_file_1.txt").exists());
    assert!(experiment_path.join("1/11/113").exists());
    assert!(!experiment_path.join("2/21").exists());
    assert!(experiment_path.join("2/22").exists());
    assert!(experiment_path.join("3").exists());
    assert!(experiment_path.join("3/test_file_2.txt").exists());
    // Delete a single file.
    let terminal_file_path = FilePath {
        path_components: vec!["3".to_string(), "test_file_2.txt".to_string()],
    };
    let req = test::TestRequest::delete()
        .uri(&format!("/api/files/experiments/{}", id))
        .set_json(terminal_file_path)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    assert!(experiment_path.join("1/11/111").exists());
    assert!(experiment_path.join("1/11/112").exists());
    assert!(experiment_path.join("1/11/112/test_file_1.txt").exists());
    assert!(experiment_path.join("1/11/113").exists());
    assert!(!experiment_path.join("2/21").exists());
    assert!(experiment_path.join("2/22").exists());
    assert!(experiment_path.join("3").exists());
    assert!(!experiment_path.join("3/test_file_2.txt").exists());
    // Delete a folder containing other files and folders.
    let super_folder_path = FilePath {
        path_components: vec!["1".to_string()],
    };
    let req = test::TestRequest::delete()
        .uri(&format!("/api/files/experiments/{}", id))
        .set_json(super_folder_path)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    assert!(!experiment_path.join("1").exists());
    assert!(!experiment_path.join("1/11").exists());
    assert!(!experiment_path.join("1/11/111").exists());
    assert!(!experiment_path.join("1/11/112").exists());
    assert!(!experiment_path.join("1/11/112/test_file_1.txt").exists());
    assert!(!experiment_path.join("1/11/113").exists());
    assert!(!experiment_path.join("2/21").exists());
    assert!(experiment_path.join("2/22").exists());
    assert!(experiment_path.join("3").exists());
    assert!(!experiment_path.join("3/test_file_2.txt").exists());
}

#[actix_web::test]
async fn test_delete_experiment_files_by_path_all() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
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
    let experiment_path = app_config.experiment_path(id.to_string());
    let folders = vec![
        experiment_path.join("1/11/111"),
        experiment_path.join("1/11/112"),
        experiment_path.join("1/11/113"),
        experiment_path.join("2/21"),
        experiment_path.join("2/22"),
        experiment_path.join("3"),
    ];
    for folder in folders {
        std::fs::create_dir_all(folder).unwrap();
    }
    std::fs::write(experiment_path.join("1/11/112/test_file_1.txt"), "test_content").unwrap();
    std::fs::write(experiment_path.join("3/test_file_2.txt"), "test_content").unwrap();
    // Assert that all files and folder exist.
    assert!(experiment_path.join("1/11/111").exists());
    assert!(experiment_path.join("1/11/112").exists());
    assert!(experiment_path.join("1/11/112/test_file_1.txt").exists());
    assert!(experiment_path.join("1/11/113").exists());
    assert!(experiment_path.join("2/21").exists());
    assert!(experiment_path.join("2/22").exists());
    assert!(experiment_path.join("3").exists());
    assert!(experiment_path.join("3/test_file_2.txt").exists());
    let root_path = FilePath {
        path_components: vec![],
    };
    let req = test::TestRequest::delete()
        .uri(&format!("/api/files/experiments/{}", id))
        .set_json(root_path)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    assert!(!experiment_path.exists());
}

#[actix_web::test]
async fn test_delete_experiment_files_non_existent() {
    let context = TestContext::new();
    let app = test::init_service(create_test_app(&context)).await;
    let file_path = FilePath {
        path_components: vec!["1".to_string(), "test_file.txt".to_string()],
    };
    let req = test::TestRequest::delete()
        .uri("/api/files/experiments/42")
        .set_json(file_path)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::NOT_FOUND);
}

#[actix_web::test]
async fn test_post_experiment_add_file_sub_folder() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
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
    let experiment_path = app_config.experiment_path(id.to_string());
    let payload_file = "../testing_resources/requests/file_upload/multipart_file_subfolders";
    let payload = std::fs::read(payload_file).unwrap();
    let content_type: mime::Mime =
        "multipart/form-data; boundary=---------------------------5851692324164894962235391524"
            .parse()
            .unwrap();
    assert!(!experiment_path.join("1").exists());
    assert!(!experiment_path.join("1/2/test_file.txt").exists());
    let req = test::TestRequest::post()
        .uri(&format!("/api/files/experiments/{}", id))
        .insert_header(ContentType(content_type))
        .set_payload(payload)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::CREATED);
    assert!(experiment_path.join("1/2/test_file.txt").exists());
    assert!(experiment_path.join("1/2/test_file.txt").is_file());
    assert_eq!(
        std::fs::read_to_string(experiment_path.join("1/2/test_file.txt")).unwrap(),
        "test content".to_string()
    );
}

#[actix_web::test]
async fn test_post_experiment_add_file_root_folder() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
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
    let experiment_path = app_config.experiment_path(id.to_string());
    let payload_file = "../testing_resources/requests/file_upload/multipart_file_root";
    let payload = std::fs::read(payload_file).unwrap();
    let content_type: mime::Mime =
        "multipart/form-data; boundary=---------------------------5851692324164894962235391524"
            .parse()
            .unwrap();
    assert!(!experiment_path.join("test_file.txt").exists());
    let req = test::TestRequest::post()
        .uri(&format!("/api/files/experiments/{}", id))
        .insert_header(ContentType(content_type))
        .set_payload(payload)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::CREATED);
    assert!(experiment_path.join("test_file.txt").exists());
    assert!(experiment_path.join("test_file.txt").is_file());
    assert_eq!(
        std::fs::read_to_string(experiment_path.join("test_file.txt")).unwrap(),
        "test content".to_string()
    );
}

#[actix_web::test]
async fn test_post_experiment_add_file_non_existent() {
    let context = TestContext::new();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
    let id = 42;
    let experiment_path = app_config.experiment_path(id.to_string());
    let payload_file = "../testing_resources/requests/file_upload/multipart_file_root";
    let payload = std::fs::read(payload_file).unwrap();
    let content_type: mime::Mime =
        "multipart/form-data; boundary=---------------------------5851692324164894962235391524"
            .parse()
            .unwrap();
    assert!(!experiment_path.join("test_file.txt").exists());
    let req = test::TestRequest::post()
        .uri(&format!("/api/files/experiments/{}", id))
        .insert_header(ContentType(content_type))
        .set_payload(payload)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::NOT_FOUND);
    assert!(!experiment_path.join("test_file.txt").exists());
}

#[actix_web::test]
async fn test_post_experiment_add_file_already_existent() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
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
    let experiment_path = app_config.experiment_path(id.to_string());
    let file_path = experiment_path.join("test_file.txt");
    let payload_file = "../testing_resources/requests/file_upload/multipart_file_root";
    let payload = std::fs::read(payload_file).unwrap();
    let content_type: mime::Mime =
        "multipart/form-data; boundary=---------------------------5851692324164894962235391524"
            .parse()
            .unwrap();
    std::fs::create_dir_all(&experiment_path).unwrap();
    let initial_content = "initial content";
    std::fs::write(&file_path, initial_content).unwrap();
    assert!(file_path.exists());
    assert!(file_path.is_file());
    let req = test::TestRequest::post()
        .uri(&format!("/api/files/experiments/{}", id))
        .insert_header(ContentType(content_type))
        .set_payload(payload)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::CONFLICT);
    assert!(file_path.exists());
    assert!(file_path.is_file());
    assert_eq!(std::fs::read_to_string(file_path).unwrap(), initial_content.to_string());
}

#[actix_web::test]
async fn test_post_experiment_add_file_empty_path() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
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
    let payload_file = "../testing_resources/requests/file_upload/multipart_file_empty_path";
    let payload = std::fs::read(payload_file).unwrap();
    let content_type: mime::Mime =
        "multipart/form-data; boundary=---------------------------5851692324164894962235391524"
            .parse()
            .unwrap();
    let req = test::TestRequest::post()
        .uri(&format!("/api/files/experiments/{}", id))
        .insert_header(ContentType(content_type))
        .set_payload(payload)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
}

#[actix_web::test]
async fn test_post_experiment_add_folder_sub_folder() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
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
    let experiment_path = app_config.experiment_path(id.to_string());
    let folder_path = FilePath {
        path_components: vec!["1".to_string(), "2".to_string(), "3".to_string()],
    };
    assert!(!experiment_path.join("1").exists());
    assert!(!experiment_path.join("1/2").exists());
    assert!(!experiment_path.join("1/2/3").exists());
    let req = test::TestRequest::post()
        .uri(&format!("/api/folders/experiments/{}", id))
        .set_json(folder_path)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::CREATED);
    assert!(experiment_path.join("1/2/3").exists());
    assert!(experiment_path.join("1/2/3").is_dir());
}

#[actix_web::test]
async fn test_post_experiment_add_folder_root() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
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
    let experiment_path = app_config.experiment_path(id.to_string());
    let folder_path = FilePath {
        path_components: vec!["1".to_string()],
    };
    assert!(!experiment_path.join("1").exists());
    let req = test::TestRequest::post()
        .uri(&format!("/api/folders/experiments/{}", id))
        .set_json(folder_path)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::CREATED);
    assert!(experiment_path.join("1").exists());
    assert!(experiment_path.join("1").is_dir());
}

#[actix_web::test]
async fn test_post_experiment_add_folder_non_existent() {
    let context = TestContext::new();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
    let id = 42;
    let experiment_path = app_config.experiment_path(id.to_string());
    let folder_path = FilePath {
        path_components: vec!["1".to_string()],
    };
    assert!(!experiment_path.join("1").exists());
    let req = test::TestRequest::post()
        .uri(&format!("/api/folders/experiments/{}", id))
        .set_json(folder_path)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::NOT_FOUND);
    assert!(!experiment_path.join("1").exists());
}

#[actix_web::test]
async fn test_post_experiment_add_folder_already_existent() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
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
    let experiment_path = app_config.experiment_path(id.to_string());
    let folder_path = experiment_path.join("1");
    let folder_path_data = FilePath {
        path_components: vec!["1".to_string()],
    };
    std::fs::create_dir_all(&folder_path).unwrap();
    assert!(folder_path.exists());
    assert!(folder_path.is_dir());
    let req = test::TestRequest::post()
        .uri(&format!("/api/folders/experiments/{}", id))
        .set_json(folder_path_data)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::CONFLICT);
    assert!(folder_path.exists());
    assert!(folder_path.is_dir());
}

#[actix_web::test]
async fn test_post_experiment_add_folder_empty_path() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
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
    let folder_path = FilePath {
        path_components: vec![],
    };
    let req = test::TestRequest::post()
        .uri(&format!("/api/folders/experiments/{}", id))
        .set_json(folder_path)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
}
