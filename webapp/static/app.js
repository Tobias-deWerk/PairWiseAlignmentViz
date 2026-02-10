const state = {
  token: null,
  globalExtent: 0,
  viewportStart: 0,
  viewportEnd: 0,
  lastPayload: null,
};

const paramIds = [
  "width",
  "height",
  "dpi",
  "min_gap_size",
  "block_size",
  "min_sequence_identity",
  "window_size",
  "tick_interval",
  "backbone_gap",
  "backbone_thickness",
  "bump_scale",
  "mismatch_line_width",
  "gap_max_height",
  "gap_width",
  "gap_height_scale",
  "indel_height_scale",
  "gap_label_size",
  "annotation_label_size",
  "annotation_thickness",
  "annotation_alpha",
  "reference_annotation_color",
  "query_annotation_color",
  "annotation_label_jitter",
  "annotation_max_layers",
  "annotation_spacing",
];

function $(id) {
  return document.getElementById(id);
}

function setStatus(text, isError = false) {
  const node = $("status");
  node.textContent = text;
  node.style.color = isError ? "#ffd9d9" : "#ecf3f7";
}

function readParams() {
  const params = {};
  for (const id of paramIds) {
    params[id] = $(id).value;
  }
  return params;
}

function buildPayload() {
  return {
    input_path: $("input_path").value.trim(),
    reference_annotation_path: $("reference_annotation_path").value.trim(),
    query_annotation_path: $("query_annotation_path").value.trim(),
    viewport_start: state.viewportStart,
    viewport_end: state.viewportEnd,
    params: readParams(),
  };
}

function updateViewportRange() {
  const viewportWidth = Math.max(1, Number($("viewport_width").value) || 1);
  const maxStart = Math.max(0, state.globalExtent - viewportWidth);
  const slider = $("viewport_start");
  slider.max = String(Math.max(0, Math.floor(maxStart)));
  slider.value = String(Math.min(maxStart, state.viewportStart));
}

async function renderCurrent() {
  const payload = buildPayload();
  const viewportWidth = Math.max(1, Number($("viewport_width").value) || 1);
  payload.viewport_start = Math.max(0, Number($("viewport_start").value) || 0);
  payload.viewport_end = payload.viewport_start + viewportWidth;

  state.viewportStart = payload.viewport_start;
  state.viewportEnd = payload.viewport_end;

  setStatus("Rendering...");
  const response = await fetch("/api/render", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "Render failed");
  }

  state.token = data.token;
  state.globalExtent = Number(data.global_extent) || 0;
  state.viewportStart = Number(data.viewport_start) || 0;
  state.viewportEnd = Number(data.viewport_end) || 0;
  state.lastPayload = payload;

  $("viewer").innerHTML = data.svg;
  updateViewportRange();
  setStatus(`Rendered. Query: ${data.query_name} | Reference: ${data.reference_name}`);
  await updateProbe();
}

async function updateProbe() {
  if (!state.token) {
    return;
  }
  const slider = $("probe_slider");
  const sliderVal = Number(slider.value);
  const sliderMax = Number(slider.max) || 1000;
  const ratio = sliderVal / sliderMax;
  const xCoord = state.viewportStart + ratio * (state.viewportEnd - state.viewportStart);

  const response = await fetch("/api/probe", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ token: state.token, x_coord: xCoord }),
  });
  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "Probe failed");
  }

  $("probe_ref_base").textContent = data.reference_base;
  $("probe_query_base").textContent = data.query_base;
  $("probe_ref_pos").textContent = data.reference_local_position;
  $("probe_query_pos").textContent = data.query_local_position;
  $("probe_idx").textContent = data.column_index;

  const flags = [];
  if (data.is_weak) flags.push("weak");
  if (data.is_mismatch) flags.push("mismatch");
  if (data.is_query_gap) flags.push("query-gap");
  if (data.is_reference_gap) flags.push("reference-gap");
  $("probe_flags").textContent = flags.length ? flags.join(", ") : "none";
  drawProbeLine(data.global_x);
}

function drawProbeLine(xValue) {
  const svg = $("viewer").querySelector("svg");
  if (!svg) {
    return;
  }
  let line = svg.querySelector("#probe-indicator-line");
  if (!line) {
    line = document.createElementNS("http://www.w3.org/2000/svg", "line");
    line.setAttribute("id", "probe-indicator-line");
    line.setAttribute("stroke", "#0ea5d7");
    line.setAttribute("stroke-width", "1.2");
    line.setAttribute("stroke-dasharray", "6 4");
    line.setAttribute("pointer-events", "none");
    svg.appendChild(line);
  }

  const viewBoxAttr = svg.getAttribute("viewBox");
  if (!viewBoxAttr) {
    return;
  }
  const parts = viewBoxAttr.split(/\s+/).map(Number);
  if (parts.length !== 4 || parts.some((v) => Number.isNaN(v))) {
    return;
  }
  const [, y, , h] = parts;
  line.setAttribute("x1", String(xValue));
  line.setAttribute("x2", String(xValue));
  line.setAttribute("y1", String(y));
  line.setAttribute("y2", String(y + h));
}

async function exportFigure(format) {
  if (!state.token) {
    throw new Error("Render before exporting.");
  }

  const response = await fetch("/api/export", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      token: state.token,
      format,
      viewport_start: state.viewportStart,
      viewport_end: state.viewportEnd,
    }),
  });

  if (!response.ok) {
    const data = await response.json();
    throw new Error(data.error || "Export failed");
  }

  const blob = await response.blob();
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = `alignment_export.${format}`;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function setupFileBrowseButtons() {
  const picker = $("file_picker");
  let targetId = null;

  document.querySelectorAll("[data-browse]").forEach((btn) => {
    btn.addEventListener("click", () => {
      targetId = btn.getAttribute("data-browse");
      picker.click();
    });
  });

  picker.addEventListener("change", () => {
    if (!targetId || !picker.files || picker.files.length === 0) {
      return;
    }
    const file = picker.files[0];
    $(targetId).value = file.name;
    setStatus(
      "Browser security only exposes filename. For local server access, paste the full absolute path.",
      false,
    );
  });
}

function bindEvents() {
  $("render_btn").addEventListener("click", async () => {
    try {
      await renderCurrent();
    } catch (err) {
      setStatus(String(err.message || err), true);
    }
  });

  $("viewport_start").addEventListener("input", async (event) => {
    state.viewportStart = Number(event.target.value) || 0;
    state.viewportEnd = state.viewportStart + Math.max(1, Number($("viewport_width").value) || 1);
    if (state.token) {
      try {
        await renderCurrent();
      } catch (err) {
        setStatus(String(err.message || err), true);
      }
    }
  });

  $("viewport_width").addEventListener("change", async () => {
    updateViewportRange();
    if (state.token) {
      try {
        await renderCurrent();
      } catch (err) {
        setStatus(String(err.message || err), true);
      }
    }
  });

  $("probe_slider").addEventListener("input", async () => {
    try {
      await updateProbe();
    } catch (err) {
      setStatus(String(err.message || err), true);
    }
  });

  $("export_svg").addEventListener("click", async () => {
    try {
      await exportFigure("svg");
      setStatus("Exported SVG.");
    } catch (err) {
      setStatus(String(err.message || err), true);
    }
  });

  $("export_png").addEventListener("click", async () => {
    try {
      await exportFigure("png");
      setStatus("Exported PNG.");
    } catch (err) {
      setStatus(String(err.message || err), true);
    }
  });
}

function bootstrap() {
  setupFileBrowseButtons();
  bindEvents();
  updateViewportRange();
}

bootstrap();
