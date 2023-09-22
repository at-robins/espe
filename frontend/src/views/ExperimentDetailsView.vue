<script setup lang="ts">
import ExperimentDetails from "@/components/experiment/ExperimentDetails.vue";
defineProps({
  id: {
    type: String,
    required: true,
  },
});
</script>

<template>
  <main>
    <experiment-details :id="id" />
  </main>
</template>
