<script setup lang="ts">
import GlobalDataDetails from "../components/global_data/GlobalDataDetails.vue";
defineProps({
  id: {
    type: String,
    required: true,
  },
});
</script>

<template>
  <main>
    <global-data-details :id="id" />
  </main>
</template>
