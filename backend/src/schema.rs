table! {
    experiment (id) {
        id -> Integer,
        experiment_name -> Text,
        mail -> Nullable<Text>,
        pipeline_id -> Nullable<Text>,
        comment -> Nullable<Text>,
        creation_time -> Timestamp,
    }
}

table! {
    global_data (id) {
        id -> Integer,
        global_data_name -> Text,
        comment -> Nullable<Text>,
        creation_time -> Timestamp,
    }
}

table! {
    pipeline_step_variable (id) {
        id -> Integer,
        experiment_id -> Integer,
        pipeline_id -> Text,
        pipeline_step_id -> Text,
        variable_id -> Text,
        variable_value -> Nullable<Text>,
        creation_time -> Timestamp,
    }
}

joinable!(pipeline_step_variable -> experiment (experiment_id));

allow_tables_to_appear_in_same_query!(
    experiment,
    global_data,
    pipeline_step_variable,
);
