table! {
    experiment (id) {
        id -> Integer,
        experiment_name -> Text,
        mail -> Text,
        pipeline_id -> Integer,
        comment -> Text,
        creation_time -> Timestamp,
    }
}

table! {
    pipeline (id) {
        id -> Integer,
        pipeline_name -> Text,
        comment -> Text,
        creation_time -> Timestamp,
    }
}

table! {
    pipeline_step (id) {
        id -> Integer,
        pipeline_id -> Integer,
        execution_type -> Text,
        execution_configuration -> Text,
        ordering -> Integer,
        creation_time -> Timestamp,
    }
}

table! {
    pipeline_step_instance (id) {
        id -> Integer,
        pipeline_step_id -> Integer,
        experiment_id -> Integer,
        pipeline_step_status -> Text,
        creation_time -> Timestamp,
    }
}

joinable!(experiment -> pipeline (pipeline_id));
joinable!(pipeline_step -> pipeline (pipeline_id));
joinable!(pipeline_step_instance -> experiment (experiment_id));
joinable!(pipeline_step_instance -> pipeline_step (pipeline_step_id));

allow_tables_to_appear_in_same_query!(
    experiment,
    pipeline,
    pipeline_step,
    pipeline_step_instance,
);
