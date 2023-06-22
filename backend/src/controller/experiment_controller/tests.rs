use crate::{test_utility::{create_test_app, TestContext}, model::db::experiment::Experiment};

use actix_web::{http::StatusCode, test};

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