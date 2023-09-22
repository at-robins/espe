<template>
  <div class="row no-wrap">
    <div class="q-ma-md">
      <q-btn round
        ><div v-if="status == ExperimentExecutionStatus.None">
          <q-icon :name="matNotStarted" color="primary" />
        </div>
        <div v-else-if="status == ExperimentExecutionStatus.Aborted">
          <q-icon :name="matStop" color="warning" />
        </div>
        <div v-else-if="status == ExperimentExecutionStatus.Failed">
          <q-icon :name="matError" color="negative" />
        </div>
        <div v-else-if="status == ExperimentExecutionStatus.Finished">
          <q-icon :name="matDone" color="positive" />
        </div>
        <div v-else-if="status == ExperimentExecutionStatus.Running">
          <q-spinner-orbit color="primary" />
        </div>
        <div v-else-if="status == ExperimentExecutionStatus.Waiting">
          <q-spinner-hourglass color="primary" />
        </div>
        <q-tooltip>The current execution status of the experiment.</q-tooltip>
      </q-btn>
    </div>
    <q-separator vertical class="q-ml-md q-mr-md" />

    <div class="col" style="display: flex; align-items: center">
      <div v-if="status == ExperimentExecutionStatus.None">
        The experiment has not yet been submitted for execution.
      </div>
      <div v-else-if="status == ExperimentExecutionStatus.Aborted">
        The pipeline was aborted.
      </div>
      <div v-else-if="status == ExperimentExecutionStatus.Failed">
        The pipeline failed.
      </div>
      <div v-else-if="status == ExperimentExecutionStatus.Finished">
        The pipeline finished successfully.
      </div>
      <div v-else-if="status == ExperimentExecutionStatus.Running">
        The pipeline is currently running.
      </div>
      <div v-else-if="status == ExperimentExecutionStatus.Waiting">
        The next pipeline step is ready and waiting for execution.
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import { type PropType } from "vue";
import {
  matDone,
  matError,
  matNotStarted,
  matStop,
} from "@quasar/extras/material-icons";
import { ExperimentExecutionStatus } from "@/scripts/types";

defineProps({
  status: {
    type: Object as PropType<ExperimentExecutionStatus>,
    required: true,
  },
});
</script>
<style scoped lang="scss"></style>
