use crate::test_utility::{create_test_app, TestContext};

use super::*;
use actix_web::{http::StatusCode, test};

#[actix_web::test]
async fn test_create_global_data() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let repo_name = "Dummy repo";
    let req = test::TestRequest::post()
        .uri("/api/globals")
        .set_json(repo_name)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::CREATED);
    let id: i32 = test::read_body_json(resp).await;
    let test_data = GlobalData::get(id, &mut connection).unwrap();
    assert_eq!(repo_name, &test_data.global_data_name);
    assert_eq!(&None, &test_data.comment);
    assert_eq!(id, test_data.id);
}

#[actix_web::test]
async fn test_create_global_data_title_too_long() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let repo_name: String = (0..513).fold(String::new(), |mut acc, _| {
        acc.push('a');
        acc
    });
    let req = test::TestRequest::post()
        .uri("/api/globals")
        .set_json(repo_name)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
    assert!(GlobalData::get_all(&mut connection).unwrap().is_empty());
}

#[actix_web::test]
async fn test_create_global_data_title_empty() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let req = test::TestRequest::post()
        .uri("/api/globals")
        .set_json("")
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
    assert!(GlobalData::get_all(&mut connection).unwrap().is_empty());
}

#[actix_web::test]
async fn test_delete_global_data() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let create_req = test::TestRequest::post()
        .uri("/api/globals")
        .set_json("Dummy repo")
        .to_request();
    let create_resp = test::call_service(&app, create_req).await;
    assert_eq!(create_resp.status(), StatusCode::CREATED);
    let id: i32 = test::read_body_json(create_resp).await;
    assert!(!GlobalData::get_all(&mut connection).unwrap().is_empty());
    let delete_req = test::TestRequest::delete()
        .uri(&format!("/api/globals/{}", id))
        .to_request();
    let delete_resp = test::call_service(&app, delete_req).await;
    assert_eq!(delete_resp.status(), StatusCode::OK);
    assert!(GlobalData::get_all(&mut connection).unwrap().is_empty());
}

#[actix_web::test]
async fn test_delete_global_data_non_existant() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    assert!(GlobalData::get_all(&mut connection).unwrap().is_empty());
    let delete_req = test::TestRequest::delete()
        .uri("/api/globals/42")
        .to_request();
    let delete_resp = test::call_service(&app, delete_req).await;
    assert_eq!(delete_resp.status(), StatusCode::NOT_FOUND);
    assert!(GlobalData::get_all(&mut connection).unwrap().is_empty());
}

#[actix_web::test]
async fn test_get_global_data() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let new_record = GlobalData {
        id,
        global_data_name: "Dummy record".to_string(),
        comment: Some("A comment".to_string()),
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::global_data::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let req = test::TestRequest::get()
        .uri(&format!("/api/globals/{}", id))
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    let fetched_data: GlobalDataDetails = test::read_body_json(resp).await;
    assert_eq!(fetched_data.id, new_record.id);
    assert_eq!(fetched_data.name, new_record.global_data_name);
    assert_eq!(fetched_data.comment, new_record.comment);
    assert_eq!(fetched_data.creation_time, new_record.creation_time);
}

#[actix_web::test]
async fn test_get_global_data_not_existant() {
    let db_context = TestContext::new();
    let app = test::init_service(create_test_app(&db_context)).await;
    let req = test::TestRequest::get().uri("/api/globals/42").to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::NOT_FOUND);
}

#[actix_web::test]
async fn test_list_global_data() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let empty_req = test::TestRequest::get().uri("/api/globals").to_request();
    let empty_resp = test::call_service(&app, empty_req).await;
    assert_eq!(empty_resp.status(), StatusCode::OK);
    let empty_fetched_data: Vec<GlobalDataDetails> = test::read_body_json(empty_resp).await;
    assert!(empty_fetched_data.is_empty());
    let new_records: Vec<GlobalData> = (0..42)
        .map(|id| GlobalData {
            id,
            global_data_name: id.to_string(),
            comment: None,
            creation_time: chrono::Utc::now().naive_local(),
        })
        .collect();
    diesel::insert_into(crate::schema::global_data::table)
        .values(&new_records)
        .execute(&mut connection)
        .unwrap();
    let req = test::TestRequest::get().uri("/api/globals").to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    let fetched_data: Vec<GlobalDataDetails> = test::read_body_json(resp).await;
    assert_eq!(new_records.len(), fetched_data.len());
    for i in 0..new_records.len() {
        assert_eq!(fetched_data[i].id, new_records[i].id);
        assert_eq!(fetched_data[i].name, new_records[i].global_data_name);
        assert_eq!(fetched_data[i].comment, new_records[i].comment);
        assert_eq!(fetched_data[i].creation_time, new_records[i].creation_time);
    }
}

#[actix_web::test]
async fn test_patch_global_data_comment() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
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
    let new_comment = Some("A comment".to_string());
    let req = test::TestRequest::patch()
        .uri(&format!("/api/globals/{}/comment", id))
        .set_json(&new_comment)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    assert_eq!(new_comment, GlobalData::get(id, &mut connection).unwrap().comment);
    let clear_comment: Option<String> = None;
    let clear_req = test::TestRequest::patch()
        .uri(&format!("/api/globals/{}/comment", id))
        .set_json(&clear_comment)
        .to_request();
    let clear_resp = test::call_service(&app, clear_req).await;
    assert_eq!(clear_resp.status(), StatusCode::OK);
    assert_eq!(clear_comment, GlobalData::get(id, &mut connection).unwrap().comment);
}

#[actix_web::test]
async fn test_patch_global_data_name() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
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
    let new_name = "A completely new name".to_string();
    let req = test::TestRequest::patch()
        .uri(&format!("/api/globals/{}/name", id))
        .set_json(&new_name)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::OK);
    assert_eq!(
        new_name,
        GlobalData::get(id, &mut connection)
            .unwrap()
            .global_data_name
    );
}

#[actix_web::test]
async fn test_patch_global_data_name_empty() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let old_name = "Dummy record".to_string();
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
    let req = test::TestRequest::patch()
        .uri(&format!("/api/globals/{}/name", id))
        .set_json("")
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
    assert_eq!(
        old_name,
        GlobalData::get(id, &mut connection)
            .unwrap()
            .global_data_name
    );
}

#[actix_web::test]
async fn test_patch_global_data_name_too_long() {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let app = test::init_service(create_test_app(&db_context)).await;
    let id = 42;
    let old_name = "Dummy record".to_string();
    let new_record = GlobalData {
        id,
        global_data_name: old_name.clone(),
        comment: None,
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::global_data::table)
        .values(&new_record)
        .execute(&mut connection)
        .unwrap();
    let new_name: String = (0..513).fold(String::new(), |mut acc, _| {
        acc.push('a');
        acc
    });
    let req = test::TestRequest::patch()
        .uri(&format!("/api/globals/{}/name", id))
        .set_json(&new_name)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);
    assert_eq!(
        old_name,
        GlobalData::get(id, &mut connection)
            .unwrap()
            .global_data_name
    );
}

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
