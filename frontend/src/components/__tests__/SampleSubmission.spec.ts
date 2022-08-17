import SampleSubmission from "../SampleSubmission.vue";
import { describe, it, expect } from "vitest";

import { mount } from "@vue/test-utils";
import { Quasar } from "quasar";

describe("SampleSubmission", () => {
  it("renders properly", () => {
    const wrapper = mount(SampleSubmission, {
      props: {},
      global: {
        plugins: [Quasar],
      },
    });
    //expect(wrapper.text()).toContain("Hello world!");
  });
});
