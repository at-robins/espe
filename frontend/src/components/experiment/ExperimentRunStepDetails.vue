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
            v-if="canBeStarted"
            :icon="matRestartAlt"
            label="Restart step"
            class="q-ml-md"
            :color="
              restartingError.has(selectedStep.id) ? 'negative' : 'positive'
            "
            :loading="isRestarting.has(selectedStep.id)"
            :disable="!areAllVariablesSet"
            @click="restartStep(selectedStep)"
          >
            <q-tooltip>
              <div
                v-if="restartingError.has(selectedStep.id)"
                class="text-black"
              >
                <error-popup
                  :error-response="restartingError.get(selectedStep.id)!"
                />
              </div>
              <div v-else-if="!areAllVariablesSet">
                All required pipeline variables must be set before execution.
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
import { type ErrorResponse } from "@/scripts/types";
import axios from "axios";
import { ref, type PropType, computed } from "vue";
import ErrorPopup from "@/components/ErrorPopup.vue";
import {
  areAllRequiredVariablesSet,
  PipelineStepStatus,
  type PipelineBlueprint,
  type PipelineStepBlueprint,
} from "@/scripts/pipeline-blueprint";
import { symOutlinedTerminal } from "@quasar/extras/material-symbols-outlined";
import { matDownload, matRestartAlt } from "@quasar/extras/material-icons";
import ExperimentStepLogs from "./ExperimentStepLogs.vue";

const isRestarting = ref(new Set<String>([]));
const restartingError = ref(new Map<String, ErrorResponse>());

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

const areAllVariablesSet = computed(() => {
  return !!props.pipeline && areAllRequiredVariablesSet(props.pipeline);
});

/**
 * Returns ```true``` if the specified step can be (re-)started.
 *
 */
const canBeStarted = computed(() => {
  if (!props.pipeline || !props.selectedStep) {
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
  const isDependecySatisfied = props.selectedStep.dependencies.every(
    (dependency) => satisfied_dependencies.includes(dependency)
  );
  return (
    props.selectedStep.status !== PipelineStepStatus.Running &&
    props.selectedStep.status !== PipelineStepStatus.Waiting &&
    isDependecySatisfied
  );
});

/**
 * Tries to restart the specified step.
 *
 * @param step the step to restart
 */
function restartStep(step: PipelineStepBlueprint | null) {
  if (step && !isRestarting.value.has(step.id)) {
    isRestarting.value.add(step.id);
    restartingError.value.delete(step.id);
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
        restartingError.value.set(step.id, error.response.data);
      })
      .finally(() => {
        isRestarting.value.delete(step.id);
      });
  }
}
</script>
<style scoped lang="scss"></style>
