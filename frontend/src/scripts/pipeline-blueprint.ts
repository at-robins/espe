export class PipelineBlueprint {
  id!: string;
  name!: string;
  description!: string;
  steps!: PipelineStepBlueprint[];
}

export class PipelineStepBlueprint {
  id!: string;
  name!: string;
  description!: string;
  container!: string;
  dependencies!: string[];
  variables!: PipelineStepBlueprintVariable[];
}

export class PipelineStepBlueprintVariable {
  id!: string;
  name!: string;
  description!: string;
  category!: PipelineStepBlueprintVariableCategory;
  required!: boolean | undefined;

  /**
   * Returns ```true``` if the variable is a boolean checkbox.
   */
  isBoolean(): boolean {
    return this.category.tag === "Boolean";
  }

  /**
   * Returns ```true``` if the variable is a global data repository reference.
   */
  isGlobal(): boolean {
    return this.category.tag === "Global";
  }

  /**
   * Returns ```true``` if the variable is a number field.
   */
  isNumber(): boolean {
    return this.category.tag === "Number";
  }

  /**
   * Returns ```true``` if the variable is a option select.
   */
  isOption(): boolean {
    return this.category.tag === "Option";
  }

  /**
   * Returns ```true``` if the variable is a string field.
   */
  isString(): boolean {
    return this.category.tag === "String";
  }
}

export class PipelineStepBlueprintVariableCategory {
  tag!: string;
  content!: unknown | undefined;
}
