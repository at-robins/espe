<script setup lang="ts">
import PipelineDetails from "../components/PipelineDetails.vue";
defineProps({
  id: {
    type: String,
    required: true,
  },
});
</script>

<template>
  <main>
    <pipeline-details :id="id" />
  </main>
</template>
