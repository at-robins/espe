<template>
  <div class="no-wrap">
    <q-tabs v-model="tab" narrow-indicator dense inline-label align="justify">
      <q-tab
        class="text-orange"
        name="run"
        :icon="symOutlinedRunCircle"
        label="Step run process"
      />
      <q-tab
        class="text-purple"
        name="build"
        :icon="symOutlinedBuildCircle"
        label="Container build process"
      />
    </q-tabs>

    <q-tab-panels v-model="tab" animated>
      <q-tab-panel name="run">
        <split-log-display v-if="logs" :log="logs.run" />
      </q-tab-panel>

      <q-tab-panel name="build">
        <split-log-display v-if="logs" :log="logs.build" />
      </q-tab-panel>
    </q-tab-panels>
  </div>
</template>

<script setup lang="ts">
import { type ErrorResponse, type ExperimentStepLogs } from "@/scripts/types";
import SplitLogDisplay from "@/components/shared/SplitLogDisplay.vue";
import { onBeforeRouteLeave, useRouter } from "vue-router";
import { ref, watch, type Ref, onBeforeUnmount } from "vue";
import axios from "axios";
import {
  symOutlinedBuildCircle,
  symOutlinedRunCircle,
} from "@quasar/extras/material-symbols-outlined";

// The intervall in which log updates are requested from the server.
const POLLING_INTERVALL_MILLISECONDS = 10000;

const isPolling = ref(false);
const pollingError: Ref<ErrorResponse | null> = ref(null);
const showPollingError = ref(false);
const router = useRouter();
const this_route = router.currentRoute.value.fullPath;
const pollingTimer: Ref<number | null> = ref(null);
const logs: Ref<ExperimentStepLogs | null> = ref(null);
const tab = ref("run");

const props = defineProps({
  experimentId: { type: String, required: true },
  pipelineId: { type: String, required: true },
  stepId: { type: String, required: true },
});

watch(
  () => props.stepId,
  () => {
    stopPolling();
    logs.value = null;
    pollLogChanges();
  },
  { immediate: true }
);

// Clears the timer if the route is changed.
onBeforeRouteLeave(() => {
  stopPolling();
});

onBeforeUnmount(() => {
  stopPolling();
});

function stopPolling() {
  if (pollingTimer.value !== null) {
    clearTimeout(pollingTimer.value);
    pollingTimer.value = null;
  }
}

/**
 * Conitinuesly polls changes from the server.
 */
function pollLogChanges() {
  if (
    !isPolling.value &&
    !pollingError.value &&
    // Stop polling if the route changes.
    router.currentRoute.value.fullPath === this_route
  ) {
    isPolling.value = true;
    pollingError.value = null;
    const config = {
      headers: {
        "content-type": "application/json",
      },
    };
    axios
      .post(
        "/api/experiments/" + props.experimentId + "/logs",
        JSON.stringify({
          pipelineId: props.pipelineId,
          stepId: props.stepId,
        }),
        config
      )
      .then((response) => {
        logs.value = response.data;
        pollingTimer.value = window.setTimeout(
          pollLogChanges,
          POLLING_INTERVALL_MILLISECONDS
        );
      })
      .catch((error) => {
        showPollingError.value = true;
        pollingError.value = error.response.data;
      })
      .finally(() => {
        isPolling.value = false;
      });
  }
}
</script>
<style scoped lang="scss"></style>
