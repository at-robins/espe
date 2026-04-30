import { expect, test } from "vitest";
import { areRequiredVariablesSet, hasRequiredVariables, isVariableSet, type PipelineStepBlueprintVariable } from "../pipeline-blueprint";

test("Test isVariableSet", () => {

  let nullVariable: PipelineStepBlueprintVariable = {
    id: "null",
    name: "null",
    description: "null",
    category: {
      tag: "String",
      content: undefined
    },
    required: false,
    value: null
  }

  let undefinedVariable: PipelineStepBlueprintVariable = {
    id: "undefined",
    name: "undefined",
    description: "undefined",
    category: {
      tag: "String",
      content: undefined
    },
    required: false,
    value: undefined
  }

  let setVariable: PipelineStepBlueprintVariable = {
    id: "set",
    name: "set",
    description: "set",
    category: {
      tag: "String",
      content: undefined
    },
    required: false,
    value: "set"
  }

  expect(isVariableSet(nullVariable)).toBe(false);
  expect(isVariableSet(undefinedVariable)).toBe(false);
  expect(isVariableSet(setVariable)).toBe(true);
});

test("Test areRequiredVariablesSet", () => {

  let nullRequired: PipelineStepBlueprintVariable = {
    id: "null",
    name: "null",
    description: "null",
    category: {
      tag: "String",
      content: undefined
    },
    required: true,
    value: null
  }

  let nullOpional: PipelineStepBlueprintVariable = {
    id: "null",
    name: "null",
    description: "null",
    category: {
      tag: "String",
      content: undefined
    },
    required: false,
    value: null
  }

  let setOptional: PipelineStepBlueprintVariable = {
    id: "set",
    name: "set",
    description: "set",
    category: {
      tag: "String",
      content: undefined
    },
    required: false,
    value: "set"
  }

  let setRequired: PipelineStepBlueprintVariable = {
    id: "set",
    name: "set",
    description: "set",
    category: {
      tag: "String",
      content: undefined
    },
    required: true,
    value: "set"
  }

  expect(areRequiredVariablesSet([nullRequired, setOptional])).toBe(false);
  expect(areRequiredVariablesSet([nullRequired, setRequired])).toBe(false);
  expect(areRequiredVariablesSet([nullOpional, setOptional])).toBe(true);
  expect(areRequiredVariablesSet([nullOpional, setRequired])).toBe(true);
});

test("Test hasRequiredVariables", () => {

  let nullRequired: PipelineStepBlueprintVariable = {
    id: "null",
    name: "null",
    description: "null",
    category: {
      tag: "String",
      content: undefined
    },
    required: true,
    value: null
  }

  let nullOpional: PipelineStepBlueprintVariable = {
    id: "null",
    name: "null",
    description: "null",
    category: {
      tag: "String",
      content: undefined
    },
    required: false,
    value: null
  }

  let setOptional: PipelineStepBlueprintVariable = {
    id: "set",
    name: "set",
    description: "set",
    category: {
      tag: "String",
      content: undefined
    },
    required: false,
    value: "set"
  }

  let setRequired: PipelineStepBlueprintVariable = {
    id: "set",
    name: "set",
    description: "set",
    category: {
      tag: "String",
      content: undefined
    },
    required: true,
    value: "set"
  }

  expect(hasRequiredVariables([nullRequired, setOptional])).toBe(true);
  expect(hasRequiredVariables([nullRequired, setRequired])).toBe(true);
  expect(hasRequiredVariables([nullOpional, setOptional])).toBe(false);
  expect(hasRequiredVariables([nullOpional, setRequired])).toBe(true);
});