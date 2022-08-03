CREATE TABLE experiment (
    id INTEGER PRIMARY KEY NOT NULL,
    experiment_name TEXT NOT NULL,
    mail TEXT NOT NULL,
    pipeline_id INTEGER NOT NULL,
    comment TEXT NOT NULL,
    creation_time DATETIME NOT NULL,
    FOREIGN KEY (pipeline_id) REFERENCES pipeline (id) ON UPDATE CASCADE ON DELETE
    SET NULL
);
CREATE TABLE pipeline (
    id INTEGER PRIMARY KEY NOT NULL,
    pipeline_name TEXT NOT NULL,
    comment TEXT NOT NULL,
    creation_time DATETIME NOT NULL
);
CREATE TABLE pipeline_step (
    id INTEGER PRIMARY KEY NOT NULL,
    pipeline_id INTEGER NOT NULL,
    execution_type TEXT NOT NULL,
    execution_configuration TEXT NOT NULL,
    ordering INTEGER NOT NULL,
    creation_time DATETIME NOT NULL,
        FOREIGN KEY (pipeline_id) REFERENCES pipeline (id) ON UPDATE CASCADE ON DELETE CASCADE
);
CREATE TABLE pipeline_step_instance (
    id INTEGER PRIMARY KEY NOT NULL,
    pipeline_step_id INTEGER NOT NULL,
    experiment_id INTEGER NOT NULL,
    pipeline_step_status TEXT NOT NULL,
    creation_time DATETIME NOT NULL,
        FOREIGN KEY (pipeline_step_id) REFERENCES pipeline_step (id) ON UPDATE CASCADE ON DELETE CASCADE,
        FOREIGN KEY (experiment_id) REFERENCES experiment (id) ON UPDATE CASCADE ON DELETE CASCADE
)