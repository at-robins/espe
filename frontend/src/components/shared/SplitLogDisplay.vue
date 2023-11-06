<template>
  <div style="width: 100%">
    <q-splitter v-model="splitterModel" :limits="[0, 100]">
      <template v-slot:before>
        <log-display header="Standard output" :body="log.stdout" />
      </template>

      <template v-slot:after>
        <log-display header="Standard error" :body="log.stderr" />
      </template>
    </q-splitter>
  </div>
</template>

<script setup lang="ts">
import { ref, type PropType, type Ref } from "vue";
import LogDisplay from "./LogDisplay.vue";
import type { ExperimentStepLog } from "@/scripts/types";

defineProps({
  log: { type: Object as PropType<ExperimentStepLog>, required: true },
});

const splitterModel: Ref<number> = ref(50);
</script>
<style scoped lang="scss"></style>
