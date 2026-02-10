const state = {
  token: null,
  globalExtent: 0,
  metadata: null,
  probeLocked: false,
  probeGlobalX: null,
  probeSvgX: null,
  probeTimer: null,
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
    params: readParams(),
  };
}

function clamp(value, minValue, maxValue) {
  return Math.min(maxValue, Math.max(minValue, value));
}

function getViewerSvg() {
  return $("viewer").querySelector("svg");
}

function updateScrollReadout() {
  const viewer = $("viewer");
  const maxScroll = Math.max(0, viewer.scrollWidth - viewer.clientWidth);
  const probeText = state.probeGlobalX == null ? "-" : state.probeGlobalX.toFixed(2);
  $("scroll_readout").textContent = `Scroll: ${viewer.scrollLeft.toFixed(0)} / ${maxScroll.toFixed(0)} | Probe x: ${probeText}`;
}

function extractSvgPoint(evt, svg) {
  if (!svg || !svg.createSVGPoint) {
    return null;
  }
  const ctm = svg.getScreenCTM();
  if (!ctm) {
    return null;
  }
  const pt = svg.createSVGPoint();
  pt.x = evt.clientX;
  pt.y = evt.clientY;
  return pt.matrixTransform(ctm.inverse());
}

function mapSvgXToDataX(svgX) {
  const m = state.metadata;
  if (!m) {
    return null;
  }
  const denom = m.axes_right_px - m.axes_left_px;
  if (denom <= 0) {
    return null;
  }
  const ratio = clamp((svgX - m.axes_left_px) / denom, 0, 1);
  return m.x_data_min + ratio * (m.x_data_max - m.x_data_min);
}

async function fetchProbeAt(xCoord, svgX) {
  const response = await fetch("/api/probe", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ token: state.token, x_coord: xCoord }),
  });
  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "Probe failed");
  }

  state.probeGlobalX = Number(data.global_x);
  state.probeSvgX = Number(svgX);

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

  drawProbeLine(svgX);
  updateScrollReadout();
}

function queueProbe(xCoord, svgX) {
  if (state.probeTimer) {
    clearTimeout(state.probeTimer);
  }
  state.probeTimer = setTimeout(async () => {
    state.probeTimer = null;
    try {
      await fetchProbeAt(xCoord, svgX);
    } catch (err) {
      setStatus(String(err.message || err), true);
    }
  }, 25);
}

function drawProbeLine(xValue) {
  const svg = getViewerSvg();
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

function attachProbeInteractions() {
  const viewer = $("viewer");
  const svg = getViewerSvg();
  if (!svg) {
    return;
  }

  viewer.addEventListener("scroll", () => {
    updateScrollReadout();
  });

  svg.addEventListener("mousemove", (evt) => {
    if (!state.token || state.probeLocked) {
      return;
    }
    const svgPoint = extractSvgPoint(evt, svg);
    if (!svgPoint) {
      return;
    }
    const dataX = mapSvgXToDataX(svgPoint.x);
    if (dataX == null) {
      return;
    }
    queueProbe(dataX, svgPoint.x);
  });

  svg.addEventListener("click", (evt) => {
    if (!state.token) {
      return;
    }
    if (state.probeLocked) {
      state.probeLocked = false;
      setStatus("Probe unlocked (hover active).", false);
      return;
    }
    const svgPoint = extractSvgPoint(evt, svg);
    if (!svgPoint) {
      return;
    }
    const dataX = mapSvgXToDataX(svgPoint.x);
    if (dataX == null) {
      return;
    }
    state.probeLocked = true;
    queueProbe(dataX, svgPoint.x);
    setStatus("Probe locked. Click again to unlock.", false);
  });

  document.addEventListener("keydown", (evt) => {
    if (evt.key === "Escape" && state.probeLocked) {
      state.probeLocked = false;
      setStatus("Probe unlocked (hover active).", false);
    }
  });

  updateScrollReadout();
}

async function renderCurrent() {
  const payload = buildPayload();

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
  state.metadata = {
    x_data_min: Number(data.x_data_min),
    x_data_max: Number(data.x_data_max),
    axes_left_px: Number(data.axes_left_px),
    axes_right_px: Number(data.axes_right_px),
    svg_width_px: Number(data.svg_width_px),
    svg_height_px: Number(data.svg_height_px),
  };
  state.probeLocked = false;
  state.probeGlobalX = null;

  $("viewer").innerHTML = data.svg;
  attachProbeInteractions();
  setStatus(`Rendered. Query: ${data.query_name} | Reference: ${data.reference_name}`);
}

async function exportFigure(format) {
  if (!state.token) {
    throw new Error("Render before exporting.");
  }

  const response = await fetch("/api/export", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ token: state.token, format }),
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
  updateScrollReadout();
}

bootstrap();
