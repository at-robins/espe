<script setup lang="ts">
import ExperimentCreation from "@/components/experiment/ExperimentCreation.vue";
import ExperimentList from "@/components/experiment/ExperimentList.vue";
</script>

<template>
  <main>
    <div class="q-ma-md">
      <q-card>
        <q-card-section>
          <div class="text-h6">Experiments</div>
        </q-card-section>
        <q-card-section class="q-pt-none">
          The Experiments section lists all user created experiments. Here new
          experiments can be created, old experiments can be accessed, modified
          or deleted. Each experiment contains intial sample information, a
          respecitve pipeline to run on the sample information and variables
          related to the specified pipeline.
        </q-card-section>
        <q-card-section>
          <experiment-creation class="q-ma-md" />
        </q-card-section>
        <q-card-section>
          <experiment-list />
        </q-card-section>
      </q-card>
    </div>
  </main>
</template>
