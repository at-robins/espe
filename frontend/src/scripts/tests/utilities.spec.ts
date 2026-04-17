import { expect, test } from "vitest";
import { is_error_response } from "../utilities";

test("Test is_error_response", () => {
  // Not an object.
  expect(is_error_response(null)).toBe(false);
  expect(is_error_response(undefined)).toBe(false);
  expect(is_error_response([1, 2, 3])).toBe(false);
  expect(is_error_response(NaN)).toBe(false);
  expect(is_error_response("test")).toBe(false);

  // Missing fields.
  expect(is_error_response({})).toBe(false);
  expect(
    is_error_response({
      uuid: "dummy ID",
      name: "Dummy name",
      message: "Dummy message",
    })
  ).toBe(false);
  expect(
    is_error_response({
      code: 1,
      name: "Dummy name",
      message: "Dummy message",
    })
  ).toBe(false);
  expect(
    is_error_response({
      code: 1,
      uuid: "dummy ID",
      message: "Dummy message",
    })
  ).toBe(false);
  expect(
    is_error_response({
      code: 1,
      uuid: "dummy ID",
      name: "Dummy name",
    })
  ).toBe(false);

  // Wrong field types.
  expect(
    is_error_response({
      code: "123",
      uuid: "dummy ID",
      name: "Dummy name",
      message: "Dummy message",
    })
  ).toBe(false);
  expect(
    is_error_response({
      code: 1,
      uuid: 1,
      name: "Dummy name",
      message: "Dummy message",
    })
  ).toBe(false);
  expect(
    is_error_response({
      code: 1,
      uuid: "dummy ID",
      name: 1,
      message: "Dummy message",
    })
  ).toBe(false);
  expect(
    is_error_response({
      code: 1,
      uuid: "dummy ID",
      name: "Dummy name",
      message: 1,
    })
  ).toBe(false);

  // Correct fields.
  expect(
    is_error_response({
      code: 1,
      uuid: "dummy ID",
      name: "Dummy name",
      message: "Dummy message",
    })
  ).toBe(true);
  expect(
    is_error_response({
      code: 1,
      uuid: "dummy ID",
      additional_field: 123,
      name: "Dummy name",
      message: "Dummy message",
    })
  ).toBe(true);
});
