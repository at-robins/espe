import type { ErrorResponse } from "./types";

/**
 * Compares two objects for shallow equality.
 */
export function equality_shallow_object(
  a: Record<string, unknown>,
  b: Record<string, unknown>
): boolean {
  const keys = Object.keys(a);
  if (keys.length !== Object.keys(b).length) {
    return false;
  }
  return keys.some((key) => a[key] !== b[key]) ? false : true;
}

/**
 * Converts an error into a string.
 *
 * @param error the error to convert into a string
 */
export function error_to_string(error: ErrorResponse): string {
  return error.uuid + " | " + error.name + ": " + error.message;
}

/**
 * Checks if the passed value is an error response.
 *
 * @param data the value to check
 */
export function is_error_response(data: any): boolean {
  return (
    data &&
    typeof data === "object" &&
    typeof data["code"] === "number" &&
    typeof data["uuid"] === "string" &&
    typeof data["name"] === "string" &&
    typeof data["message"] === "string"
  );
}
