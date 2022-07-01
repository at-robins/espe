<template>
  <q-card>
    <q-card-section class="row items-center q-pb-none">
      <div class="text-h6">
        <q-icon name="error" size="md" color="negative" class="q-pr-sm" />{{
          errorResponse.name
        }}
      </div>
      <q-space />
      <q-btn icon="close" flat round dense v-close-popup />
    </q-card-section>

    <q-card-section>
      <p>{{ errorResponse.message }}</p>
      <p>Error-ID: {{ errorResponse.uuid }}</p>
    </q-card-section>
  </q-card>
</template>

<script setup lang="ts">
import type { ErrorResponse } from "@/scripts/types";
import type { PropType } from "vue";

defineProps({
  errorResponse: { type: Object as PropType<ErrorResponse>, required: true },
});
</script>
