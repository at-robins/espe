<template>
  <div>
    <q-card>
      <q-card-section>
        <div v-if="selectedStep === null" class="text-h6">
          Select a step to display further information.
        </div>
        <div v-else class="text-h6">{{ selectedStep.name }}</div>
      </q-card-section>
      <q-card-section>
        <div v-if="selectedStep" class="q-pl-md">
          <div v-html="selectedStep.description" />
        </div>
      </q-card-section>
      <q-card-section v-if="selectedStep && pipeline">
        <q-expansion-item
          expand-separator
          :icon="symOutlinedTerminal"
          label="Display pipeline step logs"
          class="shadow-1 overflow-hidden"
          header-class="bg-secondary text-white"
          style="border-radius: 3px"
        >
          <q-card>
            <q-card-section>
              <experiment-step-logs
                :experiment-id="id"
                :pipeline-id="pipeline.id"
                :step-id="selectedStep.id"
              />
            </q-card-section>
          </q-card>
        </q-expansion-item>
      </q-card-section>
      <div v-if="selectedStep && pipeline" class="q-gutter-md q-pa-md col">
        <div class="row">
          <q-btn
            label="Download output"
            outline
            :icon="matDownload"
            color="primary"
            :disable="selectedStep.status !== PipelineStepStatus.Finished"
            :href="
              '/api/experiments/' +
              id +
              '/archive/' +
              pipeline.sanitised_id +
              selectedStep.sanitised_id
            "
          >
            <template v-slot:loading>
              <span class="block">
                <q-spinner class="on-left" />
                Generating archive
              </span>
            </template>
            <q-tooltip>
              <div>Downloads the output of the execution step.</div>
            </q-tooltip>
          </q-btn>

          <q-btn
            v-if="canBeStarted(selectedStep)"
            :icon="matRestartAlt"
            label="Restart step"
            class="q-ml-md"
            :color="restartingError ? 'negative' : 'positive'"
            :loading="isRestarting"
            @click="restartStep(selectedStep)"
          >
            <q-tooltip>
              <div v-if="restartingError" class="text-black">
                <error-popup :error-response="restartingError" />
              </div>
              <div v-else>Restarts the experiment execution step.</div>
            </q-tooltip>
          </q-btn>
        </div>
      </div>
    </q-card>
  </div>
</template>

<script setup lang="ts">
import { type ErrorResponse, type ExperimentDetails } from "@/scripts/types";
import axios from "axios";
import { ref, onMounted, type Ref, computed, type PropType, watch } from "vue";
import ErrorPopup from "@/components/ErrorPopup.vue";
import {
  PipelineStepStatus,
  type PipelineBlueprint,
  type PipelineStepBlueprint,
} from "@/scripts/pipeline-blueprint";
import { symOutlinedTerminal } from "@quasar/extras/material-symbols-outlined";
import { matDownload, matRestartAlt } from "@quasar/extras/material-icons";
import ExperimentStepLogs from "./ExperimentStepLogs.vue";
import Poller from "../shared/Poller.vue";

const isRestarting = ref(false);
const restartingError: Ref<ErrorResponse | null> = ref(null);

const props = defineProps({
  id: { type: String, required: true },
  pipeline: {
    type: Object as PropType<PipelineBlueprint | null>,
    required: false,
    default: null,
  },
  selectedStep: {
    type: Object as PropType<PipelineStepBlueprint | null>,
    required: false,
    default: null,
  },
});

// Resets the component when the step is switched.
watch(
  () => props.selectedStep,
  () => {
    isRestarting.value = false;
    restartingError.value = null;
  },
  { immediate: true }
);

/**
 * Returns ```true``` if the specified step can be (re-)started.
 *
 * @param step the step to check
 */
function canBeStarted(step: PipelineStepBlueprint | null): boolean {
  if (!props.pipeline) {
    return false;
  }
  if (!step) {
    return false;
  }
  const satisfied_dependencies = props.pipeline.steps
    .filter(
      (s) =>
        s.status === PipelineStepStatus.Finished ||
        s.status === PipelineStepStatus.Running ||
        s.status === PipelineStepStatus.Waiting
    )
    .map((s) => s.id);
  const isDependecySatisfied = step.dependencies.every((dependency) =>
    satisfied_dependencies.includes(dependency)
  );
  return (
    step.status !== PipelineStepStatus.Running &&
    step.status !== PipelineStepStatus.Waiting &&
    isDependecySatisfied
  );
}

/**
 * Tries to restart the specified step.
 *
 * @param step the step to restart
 */
function restartStep(step: PipelineStepBlueprint | null) {
  if (step && !isRestarting.value) {
    isRestarting.value = true;
    restartingError.value = null;
    const config = {
      headers: {
        "content-type": "application/json",
      },
    };
    axios
      .post(
        "/api/experiments/" + props.id + "/rerun",
        JSON.stringify(step.id),
        config
      )
      .then(() => (step.status = PipelineStepStatus.Waiting))
      .catch((error) => {
        restartingError.value = error.response.data;
      })
      .finally(() => {
        isRestarting.value = false;
      });
  }
}
</script>
<style scoped lang="scss"></style>
