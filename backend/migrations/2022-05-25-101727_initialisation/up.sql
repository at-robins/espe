CREATE TABLE experiment (
    id INTEGER PRIMARY KEY NOT NULL,
    experiment_name TEXT NOT NULL,
    mail TEXT,
    pipeline_id TEXT,
    comment TEXT,
    creation_time DATETIME NOT NULL
);
CREATE TABLE experiment_execution (
    id INTEGER PRIMARY KEY NOT NULL,
    experiment_id INTEGER NOT NULL,
    pipeline_id TEXT NOT NULL,
    pipeline_step_id TEXT NOT NULL,
    execution_status TEXT NOT NULL,
    start_time DATETIME,
    end_time DATETIME,
    creation_time DATETIME NOT NULL,
    FOREIGN KEY (experiment_id) REFERENCES experiment (id) ON UPDATE CASCADE ON DELETE CASCADE
);
CREATE TABLE global_data (
    id INTEGER PRIMARY KEY NOT NULL,
    global_data_name TEXT NOT NULL,
    comment TEXT,
    creation_time DATETIME NOT NULL
);
CREATE TABLE pipeline_step_variable (
    id INTEGER PRIMARY KEY NOT NULL,
    experiment_id INTEGER NOT NULL,
    pipeline_id TEXT NOT NULL,
    pipeline_step_id TEXT NOT NULL,
    variable_id TEXT NOT NULL,
    variable_value TEXT,
    creation_time DATETIME NOT NULL,
    FOREIGN KEY (experiment_id) REFERENCES experiment (id) ON UPDATE CASCADE ON DELETE CASCADE
);