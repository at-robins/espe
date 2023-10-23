<script setup lang="ts">
import ExperimentRunDetails from "@/components/experiment/ExperimentRunDetails.vue";
defineProps({
  id: {
    type: String,
    required: true,
  },
});
</script>

<template>
  <main>
    <experiment-run-details :id="id" />
  </main>
</template>
