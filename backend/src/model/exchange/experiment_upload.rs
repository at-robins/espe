use regex::Regex;
use serde::{Deserialize, Serialize};

const MAX_LENGTH_NAME: usize = 128;
const MAX_LENGTH_MAIL: usize = 128;
const MAX_LENGTH_COMMENT: usize = 256;

lazy_static! {
    static ref MAIL_VALIDATION_PATTERN: Regex = Regex::new("(?:[a-z0-9!#$%&'*+/=?^_`{|}~-]+(?:\\.[a-z0-9!#$%&'*+/=?^_`{|}~-]+)*|\"(?:[\\x01-\\x08\\x0b\\x0c\\x0e-\\x1f\\x21\\x23-\\x5b\\x5d-\\x7f]|\\\\[\\x01-\\x09\\x0b\\x0c\\x0e-\\x7f])*\")@(?:(?:[a-z0-9](?:[a-z0-9-]*[a-z0-9])?\\.)+[a-z0-9](?:[a-z0-9-]*[a-z0-9])?|\\[(?:(?:(2(5[0-5]|[0-4][0-9])|1[0-9][0-9]|[1-9]?[0-9]))\\.){3}(?:(2(5[0-5]|[0-4][0-9])|1[0-9][0-9]|[1-9]?[0-9])|[a-z0-9-]*[a-z0-9]:(?:[\\x01-\\x08\\x0b\\x0c\\x0e-\\x1f\\x21-\\x5a\\x53-\\x7f]|\\\\[\\x01-\\x09\\x0b\\x0c\\x0e-\\x7f])+)\\])").unwrap();
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ExperimentUpload {
    pub name: String,
    pub mail: Option<String>,
    pub comment: Option<String>,
    pub pipeline_id: i32,
}

impl ExperimentUpload {
    /**
     * Checks if the recieved upload data is valid and returns a corresponding
     * error message of not.
     */
    pub fn validate(&self) -> Result<(), String> {
        self.validate_name()
            .and(self.validate_mail())
            .and(self.validate_comment())
    }

    fn validate_name(&self) -> Result<(), String> {
        if self.name.is_empty() {
            Err("A sample name must be set.".to_string())
        } else if self.name.len() > MAX_LENGTH_NAME {
            Err(format!("The sample name may only contain {} letters.", MAX_LENGTH_NAME))
        } else {
            Ok(())
        }
    }

    fn validate_mail(&self) -> Result<(), String> {
        if let Some(email) = &self.mail {
            if !MAIL_VALIDATION_PATTERN.is_match(email) {
                return Err("The mail is invalid.".to_string());
            } else if email.len() > MAX_LENGTH_MAIL {
                return Err(format!("The mail may only contain {} letters.", MAX_LENGTH_MAIL));
            }
        }
        Ok(())
    }

    fn validate_comment(&self) -> Result<(), String> {
        if let Some(comm) = &self.comment {
            if comm.len() > MAX_LENGTH_COMMENT {
                return Err(format!(
                    "The comment may only contain {} letters.",
                    MAX_LENGTH_COMMENT
                ));
            }
        }
        Ok(())
    }
}
