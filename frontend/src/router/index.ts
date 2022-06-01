import PipelineDetailsView from "@/views/PipelineDetailsView.vue";
import { createRouter, createWebHistory } from "vue-router";
import HomeView from "../views/HomeView.vue";

const router = createRouter({
  history: createWebHistory(import.meta.env.BASE_URL),
  routes: [
    {
      path: "/",
      redirect: "/ui/",
    },
    {
      path: "/ui/",
      name: "home",
      component: HomeView,
    },
    {
      path: "/ui/pipeline/:id",
      name: "pipeline_details",
      component: PipelineDetailsView,
      props: true,
    },
  ],
});

export default router;
