<template>
  <div>
    <div class="text-h6 q-ma-md">{{ header }}</div>
    <q-separator />
    <div v-if="body === null || body === undefined" class="q-pa-md">
      No logs are currently present.
    </div>
    <div
      v-else
      v-text="body"
      class="q-pa-md"
      style="
        white-space: pre;
        overflow: auto;
        max-height: 80vh;
        font-family: 'Courier New', Courier, monospace;
      "
    />
  </div>
</template>

<script setup lang="ts">
import type { PropType } from "vue";

defineProps({
  header: { type: String, required: true },
  body: { type: String as PropType<string | null | undefined>, required: true },
});
</script>
<style scoped lang="scss"></style>
