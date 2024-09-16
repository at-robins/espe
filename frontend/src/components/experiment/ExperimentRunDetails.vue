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
      <div v-if="selectedStep" class="q-gutter-md q-pa-md col">
        <div class="row">
          <q-btn
            label="Download output"
            outline
            :icon="matDownload"
            :color="downloadError ? 'negative' : 'primary'"
            :loading="isArchiving"
            :disable="selectedStep.status !== PipelineStepStatus.Finished"
            @click="downloadStepResults(selectedStep)"
          >
            <template v-slot:loading>
              <span class="block">
                <q-spinner class="on-left" />
                Generating archive
              </span>
            </template>
            <q-tooltip>
              <div v-if="downloadError">
                <error-popup :error-response="downloadError" />
              </div>
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
              <div v-if="restartingError">
                <error-popup :error-response="restartingError" />
              </div>
              <div>Restarts the experiment execution step.</div>
            </q-tooltip>
          </q-btn>
        </div>
      </div>
    </q-card>
    <q-dialog v-model="showPollingError" v-if="pollingError">
      <error-popup :error-response="pollingError" />
    </q-dialog>
  </div>
</template>

<script setup lang="ts">
import { type ErrorResponse, type ExperimentDetails } from "@/scripts/types";
import axios from "axios";
import { ref, onMounted, type Ref, computed } from "vue";
import ErrorPopup from "@/components/ErrorPopup.vue";
import { onBeforeRouteLeave, useRouter } from "vue-router";
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
  symOutlinedTerminal,
} from "@quasar/extras/material-symbols-outlined";
import { matDownload, matRestartAlt } from "@quasar/extras/material-icons";
import ExperimentStepLogs from "./ExperimentStepLogs.vue";

// The intervall in which pipeline updates are requested from the server.
const POLLING_INTERVALL_MILLISECONDS = 10000;

const experiment: Ref<ExperimentDetails | null> = ref(null);
const pipeline: Ref<PipelineBlueprint | null> = ref(null);
const sortedSteps: Ref<PipelineStepBlueprint[][]> = ref([]);
const isLoadingPipelineDetails = ref(false);
const isRestarting = ref(false);
const loadingError: Ref<ErrorResponse | null> = ref(null);
const isPollingPipelineDetails = ref(false);
const pollingError: Ref<ErrorResponse | null> = ref(null);
const restartingError: Ref<ErrorResponse | null> = ref(null);
const selectedStep: Ref<PipelineStepBlueprint | null> = ref(null);
const showPollingError = ref(false);
const router = useRouter();
const this_route = router.currentRoute.value.fullPath;
const pollingTimer: Ref<number | null> = ref(null);
const isArchiving = ref(false);
const downloadError: Ref<ErrorResponse | null> = ref(null);

const props = defineProps({
  id: { type: String, required: true },
});

const experiment_name = computed(() => {
  return experiment.value ? experiment.value.name : props.id;
});

onMounted(() => {
  loadPipelineDetails();
});

onBeforeRouteLeave(() => {
  if (pollingTimer.value !== null) {
    clearTimeout(pollingTimer.value);
  }
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
      return axios.get("/api/experiments/" + props.id + "/run");
    })
    .then((response) => {
      setPipelineDetails(response.data);
      pollingTimer.value = window.setTimeout(
        pollDetailsChanges,
        POLLING_INTERVALL_MILLISECONDS
      );
    })
    .catch((error) => {
      pipeline.value = null;
      loadingError.value = error.response.data;
    })
    .finally(() => {
      isLoadingPipelineDetails.value = false;
    });
}

/**
 * Conitinuesly polls changes from the server.
 */
function pollDetailsChanges() {
  if (
    !isPollingPipelineDetails.value &&
    !loadingError.value &&
    !pollingError.value &&
    // Stop polling if the route changes.
    router.currentRoute.value.fullPath === this_route
  ) {
    pollingError.value = null;
    axios
      .get("/api/experiments/" + props.id + "/run")
      .then((response) => {
        setPipelineDetails(response.data);
        pollingTimer.value = window.setTimeout(
          pollDetailsChanges,
          POLLING_INTERVALL_MILLISECONDS
        );
      })
      .catch((error) => {
        showPollingError.value = true;
        pollingError.value = error.response.data;
      })
      .finally(() => {
        isPollingPipelineDetails.value = false;
      });
  }
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

/**
 * Returns ```true``` if the specified step can be (re-)started.
 *
 * @param step the step to check
 */
function canBeStarted(step: PipelineStepBlueprint | null): boolean {
  if (!pipeline.value) {
    return false;
  }
  if (!step) {
    return false;
  }
  const satisfied_dependencies = pipeline.value.steps
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

/**
 * Tries to download the specified step.
 *
 * @param step the step to download
 */
function downloadStepResults(step: PipelineStepBlueprint | null) {
  if (step && step.status == PipelineStepStatus.Finished) {
    isArchiving.value = true;
    downloadError.value = null;
    const config = {
      headers: {
        "content-type": "application/json",
      },
    };
    axios
      .post(
        "/api/experiments/" + props.id + "/archive",
        JSON.stringify(step.id),
        config
      )
      .then((response) => {
        window.location.href =
          "/api/experiments/" + props.id + "/download/" + response.data;
      })
      .catch((error) => {
        downloadError.value = error.response.data;
      })
      .finally(() => {
        isArchiving.value = false;
      });
  }
}
</script>
<style scoped lang="scss">
.chip-unselected {
  box-shadow: 0 1px 5px rgba(0, 0, 0, 0.2), 0 2px 2px rgba(0, 0, 0, 0.141),
    0 3px 1px -2px rgba(0, 0, 0, 0.122);
  transition: box-shadow 0.3s ease-in-out;
}
.chip-selected {
  box-shadow: 0 1px 5px rgba(0, 100, 255, 0.4),
    0 2px 2px rgba(0, 100, 255, 0.282), 0 3px 1px -2px rgba(0, 100, 255, 0.244);
  transition: box-shadow 0.3s ease-in-out;
}
.chip-dependency {
  box-shadow: 0 1px 5px rgba(255, 190, 0, 0.4),
    0 2px 2px rgba(255, 190, 0, 0.282), 0 3px 1px -2px rgba(255, 190, 0, 0.244);
  transition: box-shadow 0.3s ease-in-out;
}
</style>
