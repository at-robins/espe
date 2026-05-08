export type PipelineBlueprint = {
  id: string;
  sanitised_id: string;
  name: string;
  version: string;
  description: string;
  global_variables: PipelineStepBlueprintVariable[];
  steps: PipelineStepBlueprint[];
};

export type PipelineStepBlueprint = {
  id: string;
  name: string;
  sanitised_id: string;
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
 * Returns `true` if any variables in the array are required to be set before pipeline execution.
 *
 * @param pipeline_variables the variables to check
 */
export function hasRequiredVariables(
  pipeline_variables: PipelineStepBlueprintVariable[]
): boolean {
  return pipeline_variables.some((stepVar) => stepVar.required);
}

/**
 * Returns `true` if the variable value has been set.
 *
 * @param pipeline_variable the variable to check
 */
export function isVariableSet(
  pipeline_variable: PipelineStepBlueprintVariable
): boolean {
  return (
    pipeline_variable.value !== null && pipeline_variable.value !== undefined
  );
}

/**
 * Returns `true` if all variables in the array that are required have been set.
 *
 * @param pipeline_variables the variables to check
 */
export function areRequiredVariablesSet(
  pipeline_variables: PipelineStepBlueprintVariable[]
): boolean {
  return !pipeline_variables.some(
    (stepVar) => stepVar.required && !isVariableSet(stepVar)
  );
}

/**
 * Returns `true` if all variables in the array that are required have been set.
 *
 * @param pipeline_variables the variables to check
 */
export function areAllRequiredVariablesSet(
  pipeline: PipelineBlueprint | null
): boolean {
  return (
    pipeline != null &&
    areRequiredVariablesSet(pipeline.global_variables) &&
    pipeline.steps.every((step) => areRequiredVariablesSet(step.variables))
  );
}
