<template>
  <div>
    <div class="row flex-center">
      <div class="col">
        <div v-if="isBoolean(pipelineVariable)">
          <q-checkbox
            v-model="booleanModel"
            :label="pipelineVariable.name"
            toggle-indeterminate
            :color="booleanModel === null ? 'grey-4' : 'primary'"
            @update:model-value="updateModelValueBooleanNumberString"
          />
        </div>
        <div v-else-if="isGlobal(pipelineVariable)">
          <q-select
            option-label="name"
            option-value="id"
            clearable
            outlined
            v-model="globalModel"
            :options="globalOptions"
            :label="pipelineVariable.name"
            @update:model-value="updateModelValueGlobal"
          />
        </div>
        <div v-else-if="isNumber(pipelineVariable)">
          <q-input
            outlined
            v-model="numberModel"
            :label="pipelineVariable.name"
            type="number"
            @update:model-value="updateModelValueBooleanNumberString"
          >
            <template v-slot:append>
              <q-icon
                :name="matClose"
                @click="numberModel = null"
                class="cursor-pointer"
              />
            </template>
          </q-input>
        </div>
        <div v-else-if="isOption(pipelineVariable)">
          <q-select
            option-label="name"
            option-value="value"
            clearable
            outlined
            v-model="optionModel"
            :options="contentAsOptions(pipelineVariable.category)"
            :label="pipelineVariable.name"
            @update:model-value="updateModelValueOption"
          />
        </div>
        <div v-else-if="isString(pipelineVariable)">
          <q-input
            outlined
            v-model="stringModel"
            :label="pipelineVariable.name"
            @update:model-value="updateModelValueBooleanNumberString"
          >
            <template v-slot:append>
              <q-icon
                :name="matClose"
                @click="stringModel = null"
                class="cursor-pointer"
              />
            </template>
          </q-input>
        </div>
      </div>
      <div v-if="pipelineVariable.required" class="q-ml-md">
        <q-icon color="warning" :name="matPriorityHigh" size="md">
          <q-tooltip>Specifying this variable is required.</q-tooltip>
        </q-icon>
      </div>
    </div>
    <div v-if="pipelineVariable.description" class="q-mt-md q-ml-xs row">
      <div v-html="pipelineVariable.description" />
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, type PropType, type Ref, watch } from "vue";
import {
  type PipelineStepBlueprintVariable,
  contentAsOptions,
  type PipelineStepBlueprintVariableOption,
  isBoolean,
  isGlobal,
  isNumber,
  isOption,
  isString,
} from "@/scripts/pipeline-blueprint";
import { matClose, matPriorityHigh } from "@quasar/extras/material-icons";
import type { GlobalDataDetails } from "@/scripts/types";

const props = defineProps({
  pipelineVariable: {
    type: Object as PropType<PipelineStepBlueprintVariable>,
    required: true,
  },
  globalOptions: {
    type: Array as PropType<GlobalDataDetails[] | undefined>,
    default: undefined,
    required: false,
  },
});

const booleanModel: Ref<boolean | null> = ref(null);
const globalModel: Ref<GlobalDataDetails | null> = ref(null);
const numberModel: Ref<number | null> = ref(null);
const optionModel: Ref<PipelineStepBlueprintVariableOption | null> = ref(null);
const stringModel: Ref<string | null> = ref(null);

watch(props.pipelineVariable, (newValue) => {
  if (isBoolean(newValue)) {
    if (newValue === undefined || newValue.value === null) {
      booleanModel.value = null;
    } else {
      booleanModel.value = newValue.value?.toLowerCase() === "true";
    }
  } else if (isGlobal(newValue)) {
    if (newValue === undefined || newValue.value === null) {
      globalModel.value = null;
    } else {
      const query = props.globalOptions?.find(
        (option) => option.id === Number(newValue.value)
      );
      globalModel.value = query === undefined ? null : query;
    }
  } else if (isNumber(newValue)) {
    if (newValue === undefined || newValue.value === null) {
      numberModel.value = null;
    } else {
      numberModel.value = Number(newValue.value);
    }
  } else if (isOption(newValue)) {
    if (newValue === undefined || newValue.value === null) {
      optionModel.value = null;
    } else {
      const query = contentAsOptions(newValue.category).find(
        (option) => option.value === newValue.value
      );
      optionModel.value = query === undefined ? null : query;
    }
  } else if (isString(newValue)) {
    if (newValue === undefined || newValue.value === null) {
      stringModel.value = null;
    } else {
      stringModel.value = newValue.value === undefined ? null : newValue.value;
    }
  }
});

const emit = defineEmits<{
  (event: "update:modelValue", value: string | null): void;
}>();

function updateModelValueBooleanNumberString(
  newValue: boolean | number | string | null
) {
  if (newValue === null) {
    emit("update:modelValue", null);
  } else {
    emit("update:modelValue", newValue.toString());
  }
}

function updateModelValueGlobal(newValue: GlobalDataDetails | null) {
  if (newValue === null) {
    emit("update:modelValue", null);
  } else {
    emit("update:modelValue", newValue.id.toString());
  }
}

function updateModelValueOption(
  newValue: PipelineStepBlueprintVariableOption | null
) {
  if (newValue === null) {
    emit("update:modelValue", null);
  } else {
    emit("update:modelValue", newValue.value.toString());
  }
}
</script>
<style scoped lang="scss"></style>
