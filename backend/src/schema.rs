table! {
    experiment (id) {
        id -> Integer,
        experiment_name -> Text,
        mail -> Text,
        pipeline_id -> Nullable<Integer>,
        comment -> Text,
        creation_time -> Text,
    }
}

table! {
    pipeline (id) {
        id -> Integer,
        pipeline_name -> Text,
        comment -> Text,
        creation_time -> Text,
    }
}

table! {
    pipeline_step (id) {
        id -> Integer,
        pipeline_id -> Nullable<Integer>,
        execution_type -> Text,
        execution_configuration -> Text,
        ordering -> Integer,
        creation_time -> Text,
    }
}

joinable!(experiment -> pipeline (pipeline_id));
joinable!(pipeline_step -> pipeline (pipeline_id));

allow_tables_to_appear_in_same_query!(
    experiment,
    pipeline,
    pipeline_step,
);
