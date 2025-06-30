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
    <div class="row">
      <q-spinner
        v-if="logs === null"
        class="col flex-center"
        color="primary"
        size="xl"
      />
    </div>

    <q-tab-panels v-model="tab" animated>
      <q-tab-panel name="run">
        <split-log-display v-if="logs" :log="logs.run" />
      </q-tab-panel>

      <q-tab-panel name="build">
        <split-log-display v-if="logs" :log="logs.build" />
      </q-tab-panel>
    </q-tab-panels>
    <poller
      :url="logsUrl"
      :postData="postData"
      @success="logs = $event"
    ></poller>
  </div>
</template>

<script setup lang="ts">
import { type ExperimentStepLogs } from "@/scripts/types";
import SplitLogDisplay from "@/components/shared/SplitLogDisplay.vue";
import { ref, watch, type Ref, computed } from "vue";
import {
  symOutlinedBuildCircle,
  symOutlinedRunCircle,
} from "@quasar/extras/material-symbols-outlined";
import Poller from "../shared/Poller.vue";

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
    logs.value = null;
  },
  { immediate: true }
);

const logsUrl = computed(() => {
  return "/api/experiments/" + props.experimentId + "/logs";
});
const postData = computed(() => {
  return {
    config: {
      headers: {
        "content-type": "application/json",
      },
    },
    data: JSON.stringify({
      pipelineId: props.pipelineId,
      stepId: props.stepId,
    }),
  };
});
</script>
<style scoped lang="scss"></style>
