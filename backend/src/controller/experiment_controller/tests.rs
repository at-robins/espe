use crate::{
    model::db::pipeline::NewPipeline,
    test_utility::{create_test_app, TestContext},
};

use super::*;
use actix_web::{
    http::{header::ContentType, StatusCode},
    test,
};
use mime;
use serial_test::serial;

#[actix_web::test]
#[serial]
async fn test_upload_sample_post() {
    expect_dummy_request_response("sample_submission_multipart", StatusCode::CREATED).await;
}

#[actix_web::test]
#[serial]
async fn test_upload_sample_post_invalid_mail() {
    expect_dummy_request_response(
        "sample_submission_multipart_invalid_mail",
        StatusCode::BAD_REQUEST,
    )
    .await;
}

#[actix_web::test]
#[serial]
async fn test_upload_sample_post_empty_mail() {
    expect_dummy_request_response(
        "sample_submission_multipart_empty_mail",
        StatusCode::BAD_REQUEST,
    )
    .await;
}

#[actix_web::test]
#[serial]
async fn test_upload_sample_post_invalid_pipeline() {
    expect_dummy_request_response(
        "sample_submission_multipart_invalid_pipeline",
        StatusCode::BAD_REQUEST,
    )
    .await;
}

#[actix_web::test]
#[serial]
async fn test_upload_sample_post_empty_name() {
    expect_dummy_request_response(
        "sample_submission_multipart_empty_name",
        StatusCode::BAD_REQUEST,
    )
    .await;
}

#[actix_web::test]
#[serial]
async fn test_upload_sample_post_no_comment_and_mail() {
    expect_dummy_request_response(
        "sample_submission_multipart_no_comment_mail",
        StatusCode::CREATED,
    )
    .await;
}

async fn expect_dummy_request_response(file: &str, expected_code: StatusCode) {
    let db_context = TestContext::new();
    let mut connection = db_context.get_connection();
    let dummy_pipeline = NewPipeline::new("test pipeline", "test comment");
    diesel::insert_into(crate::schema::pipeline::table)
        .values(dummy_pipeline)
        .execute(&mut connection)
        .unwrap();
    let app = test::init_service(create_test_app(&db_context)).await;
    let payload_file = format!("../testing_resources/requests/sample_submission/{}", file);
    let payload = std::fs::read(payload_file).unwrap();
    let content_type: mime::Mime =
        "multipart/form-data; boundary=---------------------------5851692324164894962235391524"
            .parse()
            .unwrap();
    let req = test::TestRequest::post()
        .uri("/api/experiment")
        .insert_header(ContentType(content_type))
        .set_payload(payload)
        .to_request();
    let resp = test::call_service(&app, req).await;
    assert_eq!(resp.status(), expected_code);
}
