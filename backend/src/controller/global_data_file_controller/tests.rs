use crate::{
    model::db::pipeline::NewPipeline,
    test_utility::{create_test_app, TestContext},
};

use super::*;
use actix_web::{
    http::{header::ContentType, StatusCode},
    test,
};
use diesel::RunQueryDsl;
use mime;
use serial_test::serial;

#[actix_web::test]
async fn test_get_global_data_files() {
    let context = TestContext::new();
    let mut connection = context.get_connection();
    let app = test::init_service(create_test_app(&context)).await;
    let app_config: Configuration = (&context).into();
    let id = 42;
    let new_record = GlobalData {
        id,
        global_data_name: "Dummy record".to_string(),
        comment: None,
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::global_data::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let global_data_path = app_config.global_data_path(id.to_string());
    let folders = vec![
        global_data_path.join("1/11/111"),
        global_data_path.join("1/11/112"),
        global_data_path.join("1/11/113"),
        global_data_path.join("2/21"),
        global_data_path.join("2/22"),
        global_data_path.join("3"),
    ];
    for folder in folders {
        std::fs::create_dir_all(folder).unwrap();
    }
    let test_file_path = global_data_path.join("1/11/112/test_file.txt");
    std::fs::write(test_file_path, "test_content").unwrap();

    let req = test::TestRequest::get()
        .uri(&format!("/api/globals/{}/files", id))
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    let fetched_data: Vec<GlobalDataFileDetails> = test::read_body_json(resp).await;
    let expected_data: Vec<GlobalDataFileDetails> = vec![
        GlobalDataFileDetails {
            path_components: vec!["1".to_string()],
            is_file: false,
        },
        GlobalDataFileDetails {
            path_components: vec!["1".to_string(), "11".to_string()],
            is_file: false,
        },
        GlobalDataFileDetails {
            path_components: vec!["1".to_string(), "11".to_string(), "111".to_string()],
            is_file: false,
        },
        GlobalDataFileDetails {
            path_components: vec!["1".to_string(), "11".to_string(), "112".to_string()],
            is_file: false,
        },
        GlobalDataFileDetails {
            path_components: vec![
                "1".to_string(),
                "11".to_string(),
                "112".to_string(),
                "test_file.txt".to_string(),
            ],
            is_file: true,
        },
        GlobalDataFileDetails {
            path_components: vec!["1".to_string(), "11".to_string(), "113".to_string()],
            is_file: false,
        },
        GlobalDataFileDetails {
            path_components: vec!["2".to_string()],
            is_file: false,
        },
        GlobalDataFileDetails {
            path_components: vec!["2".to_string(), "21".to_string()],
            is_file: false,
        },
        GlobalDataFileDetails {
            path_components: vec!["2".to_string(), "22".to_string()],
            is_file: false,
        },
        GlobalDataFileDetails {
            path_components: vec!["3".to_string()],
            is_file: false,
        },
    ];
    assert_eq!(fetched_data, expected_data);
}
