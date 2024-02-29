import { createRouter, createWebHistory } from "vue-router";
import AuthorView from "../views/AuthorView.vue";
import HomeView from "../views/HomeView.vue";
import GlobalDataView from "@/views/GlobalDataView.vue";
import GlobalDataDetailsView from "@/views/GlobalDataDetailsView.vue";
import ExperimentView from "@/views/ExperimentView.vue";
import ExperimentDetailsView from "@/views/ExperimentDetailsView.vue";
import ExperimentRunDetailsView from "@/views/ExperimentRunDetailsView.vue";

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
      path: "/ui/author",
      name: "author",
      component: AuthorView,
    },
    {
      path: "/ui/globals",
      name: "globals",
      component: GlobalDataView,
    },
    {
      path: "/ui/globals/:id",
      name: "globals_detail",
      component: GlobalDataDetailsView,
      props: true,
    },
    {
      path: "/ui/experiments",
      name: "experiments",
      component: ExperimentView,
    },
    {
      path: "/ui/experiments/:id",
      name: "experiments_detail",
      component: ExperimentDetailsView,
      props: true,
    },
    {
      path: "/ui/experiments/:id/run",
      name: "experiments_run_detail",
      component: ExperimentRunDetailsView,
      props: true,
    },
  ],
});

export default router;
