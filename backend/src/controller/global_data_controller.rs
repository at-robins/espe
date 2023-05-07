use std::sync::Arc;

use crate::{
    application::{
        config::Configuration,
        error::{SeqError, SeqErrorType},
    },
    diesel::{ExpressionMethods, QueryDsl, RunQueryDsl},
    model::{
        db::global_data::{GlobalData, NewGlobalData},
        exchange::global_data_details::GlobalDataDetails,
    },
};
use actix_web::{web, HttpRequest, HttpResponse};

pub async fn create_global_data(
    request: HttpRequest,
    name: actix_web::web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let name: String = name.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    let new_record = NewGlobalData::new(name, None);
    let inserted_id: i32 = diesel::insert_into(crate::schema::global_data::table)
        .values(&new_record)
        .returning(crate::schema::global_data::id)
        .get_result(&mut connection)?;
    // Return the ID of the created global data.
    Ok(HttpResponse::Created().json(inserted_id))
}

pub async fn delete_global_data(
    request: HttpRequest,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    let exists: bool = diesel::select(diesel::dsl::exists(
        crate::schema::global_data::table.filter(crate::schema::global_data::id.eq(id)),
    ))
    .get_result(&mut connection)?;
    if exists {
        connection.immediate_transaction(|connection| {
            diesel::delete(crate::schema::global_data::table)
                .filter(crate::schema::global_data::id.eq(id))
                .execute(connection)
        })?;
        Ok(HttpResponse::Ok().finish())
    } else {
        Err(
            SeqError::new(
                "Invalid DELETE request",
                SeqErrorType::NotFoundError,
                format!("Global data with ID {} does not exist and can thereby not be deleted by request {:?}", id, request), 
                "The entity does not exist."
            )
        )
    }
}

pub async fn get_global_data(
    request: HttpRequest,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    let exists: bool = diesel::select(diesel::dsl::exists(
        crate::schema::global_data::table.filter(crate::schema::global_data::id.eq(id)),
    ))
    .get_result(&mut connection)?;
    if exists {
        let global_repo_details: GlobalDataDetails = crate::schema::global_data::table
            .find(id)
            .first::<GlobalData>(&mut connection)?
            .into();
        Ok(HttpResponse::Ok().json(global_repo_details))
    } else {
        Err(
            SeqError::new(
                "Invalid GET request",
                SeqErrorType::NotFoundError,
                format!("Global data with ID {} does not exist and can thereby not be fetched by request {:?}", id, request), 
                "The entity does not exist."
            )
        )
    }
}

pub async fn patch_global_data_name(
    request: HttpRequest,
    id: web::Path<i32>,
    new_name: web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_name = new_name.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    let exists: bool = diesel::select(diesel::dsl::exists(
        crate::schema::global_data::table.filter(crate::schema::global_data::id.eq(id)),
    ))
    .get_result(&mut connection)?;
    if exists {
        connection.immediate_transaction(|connection| {
            diesel::update(crate::schema::global_data::table.find(id))
                .set(crate::schema::global_data::global_data_name.eq(new_name))
                .execute(connection)
        })?;
        Ok(HttpResponse::Ok().finish())
    } else {
        Err(
            SeqError::new(
                "Invalid PATCH request",
                SeqErrorType::NotFoundError,
                format!("Global data with ID {} does not exist and can thereby not be patched (field: global_data_name) by request {:?}", id, request), 
                "The entity does not exist."
            )
        )
    }
}

pub async fn patch_global_data_comment(
    request: HttpRequest,
    id: web::Path<i32>,
    new_comment: web::Json<Option<String>>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_comment = new_comment.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    let exists: bool = diesel::select(diesel::dsl::exists(
        crate::schema::global_data::table.filter(crate::schema::global_data::id.eq(id)),
    ))
    .get_result(&mut connection)?;
    if exists {
        connection.immediate_transaction(|connection| {
            diesel::update(crate::schema::global_data::table.find(id))
                .set(crate::schema::global_data::comment.eq(new_comment))
                .execute(connection)
        })?;
        Ok(HttpResponse::Ok().finish())
    } else {
        Err(
            SeqError::new(
                "Invalid PATCH request",
                SeqErrorType::NotFoundError,
                format!("Global data with ID {} does not exist and can thereby not be patched (field: comment) by request {:?}", id, request), 
                "The entity does not exist."
            )
        )
    }
}

pub async fn list_global_data(
    request: HttpRequest,
) -> Result<web::Json<Vec<GlobalDataDetails>>, SeqError> {
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    let global_repos: Vec<GlobalDataDetails> = crate::schema::global_data::table
        .load::<GlobalData>(&mut connection)?
        .into_iter()
        .map(|val| val.into())
        .collect();
    Ok(web::Json(global_repos))
}
