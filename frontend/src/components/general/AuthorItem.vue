<template>
  <q-item>
    <q-item-section side>
      <q-avatar
        :icon="
          props.imageUrl && !hasImageLoadingError
            ? undefined
            : symOutlinedAccountCircle
        "
        font-size="1em"
        size="xl"
      >
        <q-img
          v-if="props.imageUrl && !hasImageLoadingError"
          :src="props.imageUrl"
          spinner-color="primary"
          ratio="1"
          @error="hasImageLoadingError = true"
        />
      </q-avatar>
    </q-item-section>
    <q-item-section>
      <q-item-label>{{ props.name }}</q-item-label>
      <div v-for="(info, infoIndex) in props.infos" :key="infoIndex">
        <q-item-label v-if="info" caption>
          <span class="row">
            <div>{{ info.name + ":&nbsp;" }}</div>
            <a v-if="info.url" :href="info.url">
              {{ info.value }}
            </a>
            <div v-else>
              {{ info.value }}
            </div>
          </span>
        </q-item-label>
      </div>
    </q-item-section>
  </q-item>
</template>

<script setup lang="ts">
import { symOutlinedAccountCircle } from "@quasar/extras/material-symbols-outlined";
import { type PropType, ref } from "vue";

const props = defineProps({
  name: { type: String, required: true },
  imageUrl: {
    type: String as PropType<string | null | undefined>,
    required: true,
  },
  infos: {
    type: Array as PropType<Array<InfoType>>,
    default: () => [],
    required: false,
  },
});

type InfoType = {
  name: string;
  value: string;
  url: string | null | undefined;
};

const hasImageLoadingError = ref(false);
</script>
<style scoped lang="scss"></style>
