<template>
  <div class="row no-wrap">
    <div class="q-ma-md">
      <q-btn :icon="symOutlinedVariables" :loading="isLoadingGlobals" round>
        <q-tooltip>Specify pipeline specific variables.</q-tooltip>
      </q-btn>
    </div>
    <q-separator vertical class="q-ml-md q-mr-md" />
    <div class="col">
      <div v-if="pipeline.global_variables.length > 0" class="row">
        <div class="col">
          <q-expansion-item expand-separator>
            <template v-slot:header>
              <q-item-section avatar>
                <q-icon
                  :name="
                    hasRequiredGlobalVariable(pipeline)
                      ? matPriorityHigh
                      : undefined
                  "
                  :color="
                    hasRequiredGlobalVariable(pipeline) ? 'warning' : 'primary'
                  "
                />
                <q-tooltip v-if="hasRequiredGlobalVariable(pipeline)">
                  This pipeline contains variables that must be specified for to
                  enusre proper execution.
                </q-tooltip>
              </q-item-section>
              <q-item-section>
                <b> Global variables </b>
              </q-item-section>
            </template>
            <div
              v-for="(
                pipelineVariable, variableIndex
              ) in pipeline.global_variables"
              :key="pipeline.id + pipelineVariable.id"
            >
              <div class="q-ma-md row">
                <entity-pipeline-variable
                  :pipeline-variable="pipelineVariable"
                  :global-options="loadedGlobalRepos"
                  @update:model-value="
                    uploadVariableGlobal(
                      pipeline.id,
                      pipelineVariable.id,
                      $event
                    )
                  "
                />
              </div>
              <q-separator
                v-if="variableIndex < pipeline.global_variables.length - 1"
                inset
              />
            </div>
          </q-expansion-item>
        </div>
      </div>
      <div
        v-for="pipelineStep in pipeline.steps"
        :key="pipelineStep.id"
        class="row"
      >
        <div v-if="pipelineStep.variables.length > 0" class="col">
          <q-expansion-item expand-separator>
            <template v-slot:header>
              <q-item-section avatar>
                <q-icon
                  :name="
                    hasRequiredStepVariable(pipelineStep)
                      ? matPriorityHigh
                      : undefined
                  "
                  :color="
                    hasRequiredStepVariable(pipelineStep)
                      ? 'warning'
                      : 'primary'
                  "
                />
                <q-tooltip v-if="hasRequiredStepVariable(pipelineStep)">
                  This step contains variables that must be specified for the
                  pipeline to work.
                </q-tooltip>
              </q-item-section>
              <q-item-section>
                <b>
                  {{ pipelineStep.name }}
                </b>
              </q-item-section>
            </template>
            <div
              v-for="(
                pipelineVariable, variableIndex
              ) in pipelineStep.variables"
              :key="pipelineStep.id + pipelineVariable.id"
            >
              <div class="q-ma-md row">
                <entity-pipeline-variable
                  :pipeline-variable="pipelineVariable"
                  :global-options="loadedGlobalRepos"
                  @update:model-value="
                    uploadVariableStep(
                      pipeline.id,
                      pipelineStep.id,
                      pipelineVariable.id,
                      $event
                    )
                  "
                />
              </div>
              <q-separator
                v-if="variableIndex < pipelineStep.variables.length - 1"
                inset
              />
            </div>
          </q-expansion-item>
        </div>
      </div>
    </div>
    <q-dialog v-model="showLoadingError" v-if="loadingError">
      <error-popup :error-response="loadingError" />
    </q-dialog>
  </div>
</template>

<script setup lang="ts">
import { ref, type PropType, type Ref, onMounted } from "vue";
import { symOutlinedVariables } from "@quasar/extras/material-symbols-outlined";
import {
  hasRequiredStepVariable,
  hasRequiredGlobalVariable,
  type PipelineBlueprint,
  type PipelineGlobalVariableUpload,
  type PipelineStepVariableUpload,
} from "@/scripts/pipeline-blueprint";
import EntityPipelineVariable from "@/components/shared/EntityPipelineVariable.vue";
import type { ErrorResponse, GlobalDataDetails } from "@/scripts/types";
import axios from "axios";
import { matPriorityHigh } from "@quasar/extras/material-icons";

const props = defineProps({
  pipeline: {
    type: Object as PropType<PipelineBlueprint>,
    required: true,
  },
  entityId: { type: Number, required: true },
  endpointType: { type: String, required: true },
});

const loadedGlobalRepos: Ref<Array<GlobalDataDetails>> = ref([]);
const isLoadingGlobals = ref(false);
const loadingError: Ref<ErrorResponse | null> = ref(null);
const showLoadingError = ref(false);

onMounted(() => {
  loadGlobalDataDetails();
});

/**
 * Initial loading of global data repo details from the server.
 */
function loadGlobalDataDetails() {
  isLoadingGlobals.value = true;
  loadingError.value = null;
  axios
    .get("/api/globals")
    .then((response) => {
      loadedGlobalRepos.value = response.data;
    })
    .catch((error) => {
      loadedGlobalRepos.value = [];
      loadingError.value = error.response.data;
      showLoadingError.value = true;
    })
    .finally(() => {
      isLoadingGlobals.value = false;
    });
}

function uploadVariableStep(
  pipelineId: string,
  pipelineStepId: string,
  variableId: string,
  variableValue: string | null
) {
  const variableUpload: PipelineStepVariableUpload = {
    pipelineId,
    pipelineStepId,
    variableId,
    variableValue,
  };
  const formData = JSON.stringify(variableUpload);
  const config = {
    headers: {
      "Content-Type": "application/json",
    },
  };
  axios
    .post(
      "/api/" + props.endpointType + "/" + props.entityId + "/variable/step",
      formData,
      config
    )
    .catch((error) => {
      loadingError.value = error.response.data;
      showLoadingError.value = true;
    });
}

function uploadVariableGlobal(
  pipelineId: string,
  variableId: string,
  variableValue: string | null
) {
  const variableUpload: PipelineGlobalVariableUpload = {
    pipelineId,
    variableId,
    variableValue,
  };
  const formData = JSON.stringify(variableUpload);
  const config = {
    headers: {
      "Content-Type": "application/json",
    },
  };
  axios
    .post(
      "/api/" + props.endpointType + "/" + props.entityId + "/variable/global",
      formData,
      config
    )
    .catch((error) => {
      loadingError.value = error.response.data;
      showLoadingError.value = true;
    });
}
</script>
<style scoped lang="scss"></style>
