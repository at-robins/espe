<template>
  <div></div>
</template>

<script setup lang="ts">
import axios from "axios";
import { ref, onMounted, type Ref, onUnmounted } from "vue";
import { matClose } from "@quasar/extras/material-icons";
import { useQuasar } from "quasar";
import { error_to_string, is_error_response } from "@/scripts/utilities";

const isStopping = ref(false);
const pollingTimer: Ref<number | null> = ref(null);
const lastPollWasError = ref(false);

const props = defineProps({
  url: { type: String, required: true },
  interval: {
    type: Number,
    default: 10000,
    required: false,
  },
});

const emit = defineEmits<{
  (event: "success", response: any): void;
}>();

const $q = useQuasar();

onUnmounted(() => {
  stopPolling();
});

onMounted(() => {
  pollChanges();
});

/**
 * Stops polling changes from the server.
 */
function stopPolling() {
  isStopping.value = true;
  if (pollingTimer.value !== null) {
    clearTimeout(pollingTimer.value);
    pollingTimer.value = null;
  }
}

/**
 * Conitinuesly polls changes from the server.
 */
function pollChanges() {
  // Stop polling if the component was shut down.
  if (!isStopping.value) {
    axios
      .get(props.url)
      .then((response) => {
        emit("success", response.data);
        // Sends a brief message that everything works again after a previous error.
        if (lastPollWasError.value) {
          $q.notify({
            message: "Retrieving data from " + props.url + " works again.",
            color: "positive",
            position: "top",
            timeout: props.interval * 0.9,
            actions: [
              {
                icon: matClose,
                color: "white",
                round: true,
                handler: () => {},
              },
            ],
          });
        }
        lastPollWasError.value = false;
      })
      .catch((error) => {
        lastPollWasError.value = true;
        let error_message =
          "Retrieving data from " +
          props.url +
          " failed. If this error persists, please check the logs for further information.";
        let error_caption = is_error_response(error.response.data)
          ? error_to_string(error.response.data)
          : String(error.response.data);
        $q.notify({
          message: error_message,
          caption: error_caption,
          color: "negative",
          position: "top",
          timeout: props.interval * 0.9,
          actions: [
            {
              icon: matClose,
              color: "white",
              round: true,
              handler: () => {},
            },
          ],
        });
      })
      .finally(() => {
        pollingTimer.value = window.setTimeout(pollChanges, props.interval);
      });
  }
}
</script>
