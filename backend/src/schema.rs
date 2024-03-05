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
    experiment_execution (id) {
        id -> Integer,
        experiment_id -> Integer,
        pipeline_id -> Text,
        pipeline_step_id -> Text,
        execution_status -> Text,
        start_time -> Nullable<Timestamp>,
        end_time -> Nullable<Timestamp>,
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
    pipeline_build_register (id) {
        id -> Integer,
        pipeline_id -> Text,
        pipeline_step_id -> Text,
        pipeline_version -> Text,
        creation_time -> Timestamp,
    }
}

table! {
    pipeline_global_variable (id) {
        id -> Integer,
        experiment_id -> Integer,
        pipeline_id -> Text,
        variable_id -> Text,
        variable_value -> Nullable<Text>,
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

joinable!(experiment_execution -> experiment (experiment_id));
joinable!(pipeline_global_variable -> experiment (experiment_id));
joinable!(pipeline_step_variable -> experiment (experiment_id));

allow_tables_to_appear_in_same_query!(
    experiment,
    experiment_execution,
    global_data,
    pipeline_build_register,
    pipeline_global_variable,
    pipeline_step_variable,
);
