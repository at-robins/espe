import SampleSubmission from "@/components/SampleSubmission.vue";
import { describe, it, expect, vi } from "vitest";

import { VueWrapper, flushPromises, mount } from "@vue/test-utils";
import { Quasar } from "quasar";
import { type ComponentPublicInstance } from "vue";
import axios from "axios";
import vitest from "vitest";

// ----------------------------
// Convenience functions
// ----------------------------

function mountWrapper(): VueWrapper<ComponentPublicInstance> {
  return mount(SampleSubmission, {
    props: {},
    global: {
      plugins: [Quasar],
    },
  });
}

function findSubmitButton(wrapper: VueWrapper<ComponentPublicInstance>) {
  return wrapper.find("button");
}

function findFileInput(wrapper: VueWrapper<ComponentPublicInstance>) {
  return wrapper.find("#sample-submission-input-file");
}

function findMailInput(wrapper: VueWrapper<ComponentPublicInstance>) {
  return wrapper.find("#sample-submission-input-mail");
}

function findNameInput(wrapper: VueWrapper<ComponentPublicInstance>) {
  return wrapper.find("#sample-submission-input-name");
}

function findAllErrors(wrapper: VueWrapper<ComponentPublicInstance>) {
  return wrapper.findAll(".q-field--error");
}

function triggerSubmit(
  wrapper: VueWrapper<ComponentPublicInstance>
): Promise<void> {
  return findSubmitButton(wrapper).trigger("submit");
}

// ----------------------------
// Mocks
// ----------------------------

vi.mock("axios");

// ----------------------------
// Test constants
// ----------------------------

const TEST_NAME = "sample name";
const TEST_FILE = new File(["content"], "file.fastq.gz");
const TEST_MAIL_VALID = "fritz.strassner@brandner.by";
const TEST_PIPELINE = {
  id: 0,
  name: "test",
  comment: "",
};

// ----------------------------
// Tests
// ----------------------------

describe("correct upload if", () => {
  it("has no errors", async () => {
    (axios.get as vitest.Mock).mockResolvedValue({
      status: 200,
      data: [TEST_PIPELINE],
    });
    (axios.post as vitest.Mock).mockResolvedValue({
      status: 201,
    });

    const wrapper = mountWrapper();

    expect(findAllErrors(wrapper)).toHaveLength(0);

    expect(axios.get).toHaveBeenCalledTimes(1);
    expect(axios.get).toHaveBeenCalledWith("/api/pipeline/blueprint");
    expect(axios.post).toHaveBeenCalledTimes(0);

    await flushPromises();

    const inputName = findNameInput(wrapper);
    await inputName.setValue(TEST_NAME);

    const inputMail = findMailInput(wrapper);
    await inputMail.setValue(TEST_MAIL_VALID);

    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (wrapper.vm as any).pipeline = TEST_PIPELINE;
    await flushPromises();

    const inputFile = findFileInput(wrapper);
    const fileElement = inputFile.element as HTMLInputElement;

    Object.defineProperty(fileElement, "files", {
      value: [TEST_FILE],
    });
    await inputFile.trigger("change");

    await triggerSubmit(wrapper);
    await flushPromises();

    expect(findAllErrors(wrapper)).toHaveLength(0);
    expect(axios.post).toHaveBeenCalledTimes(1);

    const formData = new FormData();
    formData.append("file", TEST_FILE);
    formData.append(
      "form",
      JSON.stringify({
        name: TEST_NAME,
        mail: TEST_MAIL_VALID,
        pipelineId: TEST_PIPELINE.id,
      })
    );
    const config = {
      headers: {
        "content-type": "multipart/form-data",
      },
    };

    expect(axios.post).toHaveBeenCalledWith(
      "/api/experiment",
      formData,
      config
    );
  });
});

describe("error if", () => {
  it("has an empty name field", async () => {
    const wrapper = mountWrapper();
    expect(findAllErrors(wrapper)).toHaveLength(0);
    await triggerSubmit(wrapper);
    expect(
      wrapper.findAll("label").at(0)?.classes().includes("q-field--error")
    ).toBe(true);
  });
  it("has an empty pipeline field", async () => {
    const wrapper = mountWrapper();
    expect(findAllErrors(wrapper)).toHaveLength(0);
    const inputName = findNameInput(wrapper);
    await inputName.setValue(TEST_NAME);
    await triggerSubmit(wrapper);
    expect(
      wrapper.findAll("label").at(3)?.classes().includes("q-field--error")
    ).toBe(true);
  });
  it("has an erroneous mail field", async () => {
    const wrapper = mountWrapper();
    expect(findAllErrors(wrapper)).toHaveLength(0);
    const inputName = findNameInput(wrapper);
    await inputName.setValue(TEST_NAME);
    const inputMail = findMailInput(wrapper);
    await inputMail.setValue("not an e-mail address");
    await triggerSubmit(wrapper);
    expect(
      wrapper.findAll("label").at(1)?.classes().includes("q-field--error")
    ).toBe(true);
  });
});
