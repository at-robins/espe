use crate::application::error::{SeqError, SeqErrorType};

/// The maximum length a name of an entity that is persisted in the database.
pub const MAXIMUM_ENTITY_NAME_LENGTH: usize = 512;

/// The maximum length of a comment that is persisted in the database.
const MAXIMUM_COMMENT_LENGTH: usize = 20480;

/// Validates a the specified mail address.
///
/// # Parameters
///
/// * `mail` - the mail address to validate
pub fn validate_mail<T: AsRef<str>>(mail: T) -> Result<(), SeqError> {
    let mail: &str = mail.as_ref();
    if !email_address::EmailAddress::is_valid(mail) {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            format!("The mail address {} is invalid.", mail),
            "The mail address is invalid.",
        ));
    }
    Ok(())
}

/// Validates a the specified comment.
///
/// # Parameters
///
/// * `comment` - the comment to validate
pub fn validate_comment<T: AsRef<str>>(comment: T) -> Result<(), SeqError> {
    let comment: &str = comment.as_ref();
    if comment.len() > MAXIMUM_COMMENT_LENGTH {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            format!(
                "The comment {} exceeds the limit of {} characters.",
                comment, MAXIMUM_COMMENT_LENGTH
            ),
            "The comment is invalid.",
        ));
    }
    Ok(())
}

/// Validates the name of a databse entity.
///
/// # Parameters
///
/// * `name` - the entity name to validate
pub fn validate_entity_name<T: AsRef<str>>(name: T) -> Result<(), SeqError> {
    let name: &str = name.as_ref();
    if name.is_empty() {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            "The name may not be empty.",
            "The name is invalid.",
        ));
    }
    if name.len() > MAXIMUM_ENTITY_NAME_LENGTH {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            format!(
                "The name {} exceeds the limit of {} characters.",
                name, MAXIMUM_ENTITY_NAME_LENGTH
            ),
            "The name is invalid.",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_validate_entity_name() {
        assert!(validate_entity_name("Valid entity name").is_ok());
        assert!(validate_entity_name("").is_err());
        let name_too_long: String =
            (0..=MAXIMUM_ENTITY_NAME_LENGTH).fold(String::new(), |mut acc, _| {
                acc.push('a');
                acc
            });
        assert!(validate_entity_name(name_too_long).is_err());
    }

    #[test]
    fn test_validate_mail() {
        assert!(validate_mail("a@bc.de").is_ok());
        assert!(validate_mail("").is_err());
        let mail_component_too_long: String = (0..=254).fold(String::new(), |mut acc, _| {
            acc.push('a');
            acc
        });
        let mail_too_long = format!("{}@{}.com", mail_component_too_long, mail_component_too_long);
        assert!(validate_mail(mail_too_long).is_err());
        assert!(validate_mail("an@invalid@mail.address").is_err());
        assert!(validate_mail("no.separator.mail.address").is_err());
    }

    #[test]
    fn test_validate_comment() {
        assert!(validate_comment("Valid comment").is_ok());
        assert!(validate_comment("<div>Another valid comment</div>").is_ok());
        assert!(validate_comment("").is_ok());
        let name_too_long: String =
            (0..=MAXIMUM_COMMENT_LENGTH).fold(String::new(), |mut acc, _| {
                acc.push('a');
                acc
            });
        assert!(validate_comment(name_too_long).is_err());
    }
}
