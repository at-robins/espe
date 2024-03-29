export type PipelineBlueprint = {
  id: string;
  name: string;
  version: string;
  description: string;
  global_variables: PipelineStepBlueprintVariable[];
  steps: PipelineStepBlueprint[];
};

export type PipelineStepBlueprint = {
  id: string;
  name: string;
  description: string;
  container: string;
  dependencies: string[];
  variables: PipelineStepBlueprintVariable[];
  status: PipelineStepStatus | null | undefined;
};

export type PipelineStepBlueprintVariable = {
  id: string;
  name: string;
  description: string;
  category: PipelineStepBlueprintVariableCategory;
  required: boolean | null | undefined;
  value: string | null | undefined;
};

export type PipelineStepBlueprintVariableCategory = {
  tag: string;
  content: unknown | undefined;
};

export type PipelineStepBlueprintVariableOption = {
  name: string;
  value: string;
};

export type PipelineStepVariableUpload = {
  pipelineId: string;
  pipelineStepId: string;
  variableId: string;
  variableValue: string | null | undefined;
};

export type PipelineGlobalVariableUpload = {
  pipelineId: string;
  variableId: string;
  variableValue: string | null | undefined;
};

export enum PipelineStepStatus {
  Aborted = "Aborted",
  Failed = "Failed",
  Finished = "Finished",
  Running = "Running",
  Waiting = "Waiting",
}

/**
 * Returns ```true``` if the variable is a boolean checkbox.
 *
 * @param variable the variable to check the type of
 */
export function isBoolean(variable: PipelineStepBlueprintVariable): boolean {
  return variable.category.tag === "Boolean";
}

/**
 * Returns ```true``` if the variable is a global data repository reference.
 *
 * @param variable the variable to check the type of
 */
export function isGlobal(variable: PipelineStepBlueprintVariable): boolean {
  return variable.category.tag === "Global";
}

/**
 * Returns ```true``` if the variable is a number field.
 *
 * @param variable the variable to check the type of
 */
export function isNumber(variable: PipelineStepBlueprintVariable): boolean {
  return variable.category.tag === "Number";
}

/**
 * Returns ```true``` if the variable is a option select.
 *
 * @param variable the variable to check the type of
 */
export function isOption(variable: PipelineStepBlueprintVariable): boolean {
  return variable.category.tag === "Option";
}

/**
 * Returns ```true``` if the variable is a string field.
 *
 * @param variable the variable to check the type of
 */
export function isString(variable: PipelineStepBlueprintVariable): boolean {
  return variable.category.tag === "String";
}

/**
 * Returns the content of a variable category as array of options.
 *
 * @param category the category to convert
 */
export function contentAsOptions(
  category: PipelineStepBlueprintVariableCategory
): PipelineStepBlueprintVariableOption[] {
  return category.content as PipelineStepBlueprintVariableOption[];
}

/**
 * Returns `true` if any of the pipeline step variables is required.
 *
 * @param step the step to check for requirements
 */
export function hasRequiredStepVariable(step: PipelineStepBlueprint): boolean {
  return step.variables.some((stepVar) => stepVar.required);
}

/**
 * Returns `true` if any of the global pipeline variables is required.
 *
 * @param pipeline the pipeline to check for requirements
 */
export function hasRequiredGlobalVariable(
  pipeline: PipelineBlueprint
): boolean {
  return pipeline.global_variables.some((globalVar) => globalVar.required);
}
