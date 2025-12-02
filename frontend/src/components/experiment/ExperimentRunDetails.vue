<template>
  <div class="q-pa-md q-gutter-md">
    <q-card>
      <q-card-section>
        <div class="text-h6">Experiment: {{ experiment_name }}</div>
      </q-card-section>
      <div class="q-pa-md gutter-md no-wrap row" v-if="!loadingError">
        <div
          v-if="isLoadingPipelineDetails"
          class="flex-center col q-ma-lg"
          style="display: flex"
        >
          <q-spinner size="xl" color="primary" />
        </div>
        <div class="row no-wrap q-pa-xs" style="overflow: auto">
          <div
            v-for="(stepGroup, indexGroup) in sortedSteps"
            :key="
              stepGroup.map((s) => s.id).reduce((p, c) => p + c, 'stepGroup')
            "
            class="q-pb-md row no-wrap"
          >
            <div class="row flex-center no-wrap">
              <div class="col">
                <div
                  v-for="(step, indexStep) in stepGroup"
                  :key="step.id"
                  :class="{ 'q-pb-md': indexStep < stepGroup.length - 1 }"
                >
                  <q-btn
                    rounded
                    flat
                    no-caps
                    no-wrap
                    @click="selectStep(step)"
                    :class="getChipColour(step)"
                  >
                    <div v-if="!step.status">
                      <q-icon
                        :name="symOutlinedNotStarted"
                        color="primary"
                        left
                      />
                    </div>
                    <div v-else-if="step.status == PipelineStepStatus.Aborted">
                      <q-icon
                        :name="symOutlinedStopCircle"
                        color="warning"
                        left
                      />
                    </div>
                    <div v-else-if="step.status == PipelineStepStatus.Failed">
                      <q-icon :name="symOutlinedError" color="negative" left />
                    </div>
                    <div v-else-if="step.status == PipelineStepStatus.Finished">
                      <q-icon
                        :name="symOutlinedCheckCircle"
                        color="positive"
                        left
                      />
                    </div>
                    <div v-else-if="step.status == PipelineStepStatus.Running">
                      <q-spinner-orbit color="primary" class="on-left" />
                    </div>
                    <div v-else-if="step.status == PipelineStepStatus.Waiting">
                      <q-icon
                        :name="symOutlinedSchedule"
                        color="primary"
                        left
                      />
                    </div>
                    <div class="text-center">{{ step.name }}</div>
                  </q-btn>
                </div>
              </div>
              <div v-if="indexGroup < sortedSteps.length - 1" class="q-ma-md">
                <q-icon name="trending_flat" size="lg" />
              </div>
            </div>
          </div>
        </div>
      </div>
      <div v-else>
        <error-popup :error-response="loadingError" />
      </div>
    </q-card>
    <experiment-run-step-details :id=id :pipeline="pipeline" :selected-step="selectedStep"/>
    <poller
      v-if="enableRunDetailsPoller"
      :url="run_details_url"
      @success="setPipelineDetails"
    ></poller>
  </div>
</template>

<script setup lang="ts">
import { type ErrorResponse, type ExperimentDetails } from "@/scripts/types";
import axios from "axios";
import { ref, onMounted, type Ref, computed } from "vue";
import ErrorPopup from "@/components/ErrorPopup.vue";
import {
  PipelineStepStatus,
  type PipelineBlueprint,
  type PipelineStepBlueprint,
} from "@/scripts/pipeline-blueprint";
import {
  symOutlinedCheckCircle,
  symOutlinedError,
  symOutlinedNotStarted,
  symOutlinedSchedule,
  symOutlinedStopCircle,
} from "@quasar/extras/material-symbols-outlined";
import ExperimentRunStepDetails from "./ExperimentRunStepDetails.vue";
import Poller from "../shared/Poller.vue";

const experiment: Ref<ExperimentDetails | null> = ref(null);
const pipeline: Ref<PipelineBlueprint | null> = ref(null);
const sortedSteps: Ref<PipelineStepBlueprint[][]> = ref([]);
const isLoadingPipelineDetails = ref(false);
const enableRunDetailsPoller = ref(false);
const loadingError: Ref<ErrorResponse | null> = ref(null);
const selectedStep: Ref<PipelineStepBlueprint | null> = ref(null);

const props = defineProps({
  id: { type: String, required: true },
});

const experiment_name = computed(() => {
  return experiment.value ? experiment.value.name : props.id;
});
const run_details_url = computed(() => {
  return "/api/experiments/" + props.id + "/run";
});

onMounted(() => {
  loadPipelineDetails();
});

/**
 * Initial loading of details from the server.
 */
function loadPipelineDetails() {
  isLoadingPipelineDetails.value = true;
  loadingError.value = null;
  axios
    .get("/api/experiments/" + props.id)
    .then((response) => {
      experiment.value = response.data;
      enableRunDetailsPoller.value = true;
    })
    .catch((error) => {
      pipeline.value = null;
      loadingError.value = error.response.data;
      enableRunDetailsPoller.value = false;
    })
    .finally(() => {
      isLoadingPipelineDetails.value = false;
    });
}

function setPipelineDetails(response: PipelineBlueprint | null) {
  pipeline.value = response;
  // Groupes pipeline steps based on dependencies.
  // This sorting algorithm has a bad time complexity, but since it is
  // executed asynchronously it does not really matter.
  if (pipeline.value) {
    const stepsByDependency: PipelineStepBlueprint[][] = [];
    const satisfiedDependencies: string[] = [];
    let remainingSteps = [...pipeline.value.steps];
    while (remainingSteps.length > 0) {
      const numberOfRemainingSteps = remainingSteps.length;
      // Obtaines steps with satisfied dependencies.
      const steps_with_satisfied_dependencies = remainingSteps.filter((step) =>
        step.dependencies.every((dependency) =>
          satisfiedDependencies.includes(dependency)
        )
      );
      // Removes the obtained steps from the remaining steps.
      remainingSteps = remainingSteps.filter(
        (step) =>
          !step.dependencies.every((dependency) =>
            satisfiedDependencies.includes(dependency)
          )
      );
      // Updates the dependencies which have already been
      /// satisfied with the newly obtained values.
      for (const step of steps_with_satisfied_dependencies) {
        satisfiedDependencies.push(step.id);
      }
      stepsByDependency.push(steps_with_satisfied_dependencies);
      // If for any reason there remain invalid pipeline steps that
      // have dependiencies, which cannot be satisfied, append all
      // of them to the end.
      if (numberOfRemainingSteps == remainingSteps.length) {
        stepsByDependency.push(remainingSteps);
        break;
      }
    }
    // Update at the end to avoid inconsitent UI state.
    sortedSteps.value = stepsByDependency;
  } else {
    sortedSteps.value = [];
  }
  if (selectedStep.value) {
    const update_selected = response?.steps.find(
      (detail) => detail.id === selectedStep.value?.id
    );
    selectedStep.value = !update_selected ? null : update_selected;
  }
}

function getChipColour(step: PipelineStepBlueprint) {
  if (!selectedStep.value) {
    return "chip-unselected";
  } else if (step === selectedStep.value) {
    return "chip-selected";
  } else if (selectedStep.value.dependencies.includes(step.id)) {
    return "chip-dependency";
  } else {
    return "chip-unselected";
  }
}

/**
 * Selects the specified pipeline step to display related information.
 */
function selectStep(step: PipelineStepBlueprint) {
  if (selectedStep.value === step) {
    selectedStep.value = null;
  } else {
    selectedStep.value = step;
  }
}
</script>
<style scoped lang="scss">
.chip-unselected {
  box-shadow:
    0 1px 5px rgba(0, 0, 0, 0.2),
    0 2px 2px rgba(0, 0, 0, 0.141),
    0 3px 1px -2px rgba(0, 0, 0, 0.122);
  transition: box-shadow 0.3s ease-in-out;
}
.chip-selected {
  box-shadow:
    0 1px 5px rgba(0, 100, 255, 0.4),
    0 2px 2px rgba(0, 100, 255, 0.282),
    0 3px 1px -2px rgba(0, 100, 255, 0.244);
  transition: box-shadow 0.3s ease-in-out;
}
.chip-dependency {
  box-shadow:
    0 1px 5px rgba(255, 190, 0, 0.4),
    0 2px 2px rgba(255, 190, 0, 0.282),
    0 3px 1px -2px rgba(255, 190, 0, 0.244);
  transition: box-shadow 0.3s ease-in-out;
}
</style>
