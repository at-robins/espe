<template>
  <div class="row no-wrap">
    <div class="q-ma-md">
      <q-btn :icon="symOutlinedVariables" :loading="isLoadingGlobals" round>
        <q-tooltip>Specify pipeline specific variables.</q-tooltip>
      </q-btn>
    </div>
    <q-separator vertical class="q-ml-md q-mr-md" />
    <div class="col">
      <div
        v-for="pipelineStep in pipeline.steps"
        :key="pipelineStep.id"
        class="row"
      >
        <div class="col">
          <q-expansion-item expand-separator>
            <template v-slot:header>
              <q-item-section avatar>
                <q-icon
                  :name="
                    hasRequiredVariable(pipelineStep)
                      ? matPriorityHigh
                      : undefined
                  "
                  :color="
                    hasRequiredVariable(pipelineStep) ? 'warning' : 'primary'
                  "
                />
                <q-tooltip>
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
  hasRequiredVariable,
  type PipelineBlueprint,
} from "@/scripts/pipeline-blueprint";
import EntityPipelineVariable from "@/components/shared/EntityPipelineVariable.vue";
import type { ErrorResponse, GlobalDataDetails } from "@/scripts/types";
import axios from "axios";
import { matPriorityHigh } from "@quasar/extras/material-icons";

defineProps({
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
</script>
<style scoped lang="scss"></style>
