import { expect, test } from "vitest";
import {
  areAllRequiredVariablesSet,
  areRequiredVariablesSet,
  hasRequiredVariables,
  isVariableSet,
  PipelineStepStatus,
  type PipelineBlueprint,
  type PipelineStepBlueprintVariable,
} from "../pipeline-blueprint";

test("Test isVariableSet", () => {
  let nullVariable: PipelineStepBlueprintVariable = {
    id: "null",
    name: "null",
    description: "null",
    category: {
      tag: "String",
      content: undefined,
    },
    required: false,
    value: null,
  };

  let undefinedVariable: PipelineStepBlueprintVariable = {
    id: "undefined",
    name: "undefined",
    description: "undefined",
    category: {
      tag: "String",
      content: undefined,
    },
    required: false,
    value: undefined,
  };

  let setVariable: PipelineStepBlueprintVariable = {
    id: "set",
    name: "set",
    description: "set",
    category: {
      tag: "String",
      content: undefined,
    },
    required: false,
    value: "set",
  };

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
      content: undefined,
    },
    required: true,
    value: null,
  };

  let nullOpional: PipelineStepBlueprintVariable = {
    id: "null",
    name: "null",
    description: "null",
    category: {
      tag: "String",
      content: undefined,
    },
    required: false,
    value: null,
  };

  let setOptional: PipelineStepBlueprintVariable = {
    id: "set",
    name: "set",
    description: "set",
    category: {
      tag: "String",
      content: undefined,
    },
    required: false,
    value: "set",
  };

  let setRequired: PipelineStepBlueprintVariable = {
    id: "set",
    name: "set",
    description: "set",
    category: {
      tag: "String",
      content: undefined,
    },
    required: true,
    value: "set",
  };

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
      content: undefined,
    },
    required: true,
    value: null,
  };

  let nullOpional: PipelineStepBlueprintVariable = {
    id: "null",
    name: "null",
    description: "null",
    category: {
      tag: "String",
      content: undefined,
    },
    required: false,
    value: null,
  };

  let setOptional: PipelineStepBlueprintVariable = {
    id: "set",
    name: "set",
    description: "set",
    category: {
      tag: "String",
      content: undefined,
    },
    required: false,
    value: "set",
  };

  let setRequired: PipelineStepBlueprintVariable = {
    id: "set",
    name: "set",
    description: "set",
    category: {
      tag: "String",
      content: undefined,
    },
    required: true,
    value: "set",
  };

  expect(hasRequiredVariables([nullRequired, setOptional])).toBe(true);
  expect(hasRequiredVariables([nullRequired, setRequired])).toBe(true);
  expect(hasRequiredVariables([nullOpional, setOptional])).toBe(false);
  expect(hasRequiredVariables([nullOpional, setRequired])).toBe(true);
});

test("Test areAllRequiredVariablesSet: no pipeline", () => {
  expect(areAllRequiredVariablesSet(null)).toBe(false);
});

test("Test areAllRequiredVariablesSet: no required variables", () => {
  let testPipeline: PipelineBlueprint = {
    id: "pipeline1",
    sanitised_id: "pipeline1",
    name: "pipeline1",
    version: "pipeline1",
    description: "pipeline1",
    global_variables: [
      {
        id: "global1",
        name: "global1",
        description: "",
        category: {
          tag: "String",
          content: undefined,
        },
        required: false,
        value: "testval",
      },
      {
        id: "global2",
        name: "global2",
        description: "",
        category: {
          tag: "String",
          content: undefined,
        },
        required: false,
        value: null,
      },
    ],
    steps: [
      {
        id: "step1",
        name: "step1",
        sanitised_id: "step1",
        description: "step1",
        container: "step1",
        dependencies: [],
        variables: [
          {
            id: "step1var1",
            name: "step1var1",
            description: "",
            category: {
              tag: "String",
              content: undefined,
            },
            required: false,
            value: null,
          },
        ],
        status: PipelineStepStatus.Waiting,
      },
      {
        id: "step2",
        name: "step2",
        sanitised_id: "step2",
        description: "step2",
        container: "step2",
        dependencies: [],
        variables: [
          {
            id: "step2var1",
            name: "step2var1",
            description: "",
            category: {
              tag: "String",
              content: undefined,
            },
            required: false,
            value: null,
          },
          {
            id: "step2var2",
            name: "step2var2",
            description: "",
            category: {
              tag: "String",
              content: undefined,
            },
            required: false,
            value: null,
          },
        ],
        status: PipelineStepStatus.Waiting,
      },
    ],
  };

  expect(areAllRequiredVariablesSet(testPipeline)).toBe(true);
});

test("Test areAllRequiredVariablesSet: required global variables", () => {
  let testPipeline: PipelineBlueprint = {
    id: "pipeline1",
    sanitised_id: "pipeline1",
    name: "pipeline1",
    version: "pipeline1",
    description: "pipeline1",
    global_variables: [
      {
        id: "global1",
        name: "global1",
        description: "",
        category: {
          tag: "String",
          content: undefined,
        },
        required: false,
        value: "testval",
      },
      {
        id: "global2",
        name: "global2",
        description: "",
        category: {
          tag: "String",
          content: undefined,
        },
        required: true,
        value: null,
      },
    ],
    steps: [
      {
        id: "step1",
        name: "step1",
        sanitised_id: "step1",
        description: "step1",
        container: "step1",
        dependencies: [],
        variables: [
          {
            id: "step1var1",
            name: "step1var1",
            description: "",
            category: {
              tag: "String",
              content: undefined,
            },
            required: false,
            value: null,
          },
        ],
        status: PipelineStepStatus.Waiting,
      },
      {
        id: "step2",
        name: "step2",
        sanitised_id: "step2",
        description: "step2",
        container: "step2",
        dependencies: [],
        variables: [
          {
            id: "step2var1",
            name: "step2var1",
            description: "",
            category: {
              tag: "String",
              content: undefined,
            },
            required: false,
            value: null,
          },
          {
            id: "step2var2",
            name: "step2var2",
            description: "",
            category: {
              tag: "String",
              content: undefined,
            },
            required: false,
            value: null,
          },
        ],
        status: PipelineStepStatus.Waiting,
      },
    ],
  };

  expect(areAllRequiredVariablesSet(testPipeline)).toBe(false);
  testPipeline.global_variables[1].value = "testval";
  expect(areAllRequiredVariablesSet(testPipeline)).toBe(true);
});

test("Test areAllRequiredVariablesSet: required step variables", () => {
  let testPipeline: PipelineBlueprint = {
    id: "pipeline1",
    sanitised_id: "pipeline1",
    name: "pipeline1",
    version: "pipeline1",
    description: "pipeline1",
    global_variables: [
      {
        id: "global1",
        name: "global1",
        description: "",
        category: {
          tag: "String",
          content: undefined,
        },
        required: false,
        value: "testval",
      },
      {
        id: "global2",
        name: "global2",
        description: "",
        category: {
          tag: "String",
          content: undefined,
        },
        required: false,
        value: null,
      },
    ],
    steps: [
      {
        id: "step1",
        name: "step1",
        sanitised_id: "step1",
        description: "step1",
        container: "step1",
        dependencies: [],
        variables: [
          {
            id: "step1var1",
            name: "step1var1",
            description: "",
            category: {
              tag: "String",
              content: undefined,
            },
            required: false,
            value: null,
          },
        ],
        status: PipelineStepStatus.Waiting,
      },
      {
        id: "step2",
        name: "step2",
        sanitised_id: "step2",
        description: "step2",
        container: "step2",
        dependencies: [],
        variables: [
          {
            id: "step2var1",
            name: "step2var1",
            description: "",
            category: {
              tag: "String",
              content: undefined,
            },
            required: false,
            value: null,
          },
          {
            id: "step2var2",
            name: "step2var2",
            description: "",
            category: {
              tag: "String",
              content: undefined,
            },
            required: true,
            value: null,
          },
        ],
        status: PipelineStepStatus.Waiting,
      },
    ],
  };

  expect(areAllRequiredVariablesSet(testPipeline)).toBe(false);
  testPipeline.steps[1].variables[1].value = "testval";
  expect(areAllRequiredVariablesSet(testPipeline)).toBe(true);
});
