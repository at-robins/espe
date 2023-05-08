<script setup lang="ts">
import GlobalDataCreation from "@/components/GlobalDataCreation.vue";
import GlobalDataList from "@/components/GlobalDataList.vue";
</script>

<template>
  <main>
    <div class="q-ma-md">
      <q-card flat bordered class="my-card">
        <q-card-section>
          <div class="text-h6">Globla data repositories</div>
        </q-card-section>
        <q-card-section class="q-pt-none">
          Global data repositories contain data that needs to be shared between
          multiple pipeline instances such as genomes that are used for read
          alignment.
        </q-card-section>
      </q-card>
    </div>
    <global-data-creation class="q-ma-md" />
    <global-data-list />
  </main>
</template>
