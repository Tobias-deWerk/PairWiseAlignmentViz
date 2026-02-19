const state = {
  token: null,
  globalExtent: 0,
  metadata: null,
  probeGlobalX: null,
  probeSvgX: null,
  probeTimer: null,
  probeData: null,
  selection: null,
  selectionAnchor: null,
  dragAnchor: null,
  selectionDraft: null,
  dragActive: false,
  dragMoved: false,
  suppressNextClick: false,
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

const PRESET_VALUES = Object.freeze({
  short_lt_5000: Object.freeze({
    width: "15",
    height: "4",
    min_gap_size: "1",
    block_size: "50",
    min_sequence_identity: "0.3",
    backbone_gap: "0.2",
    bump_scale: "0.1",
    gap_max_height: "4",
    gap_height_scale: "0.03",
    gap_label_size: "10",
    backbone_thickness: "3",
    gap_width: "30",
    indel_height_scale: "0.01",
  }),
  long_ge_5000: Object.freeze({
    width: "auto",
    height: "6.0",
    min_gap_size: "10",
    block_size: "100",
    min_sequence_identity: "0.7",
    backbone_gap: "0.5",
    bump_scale: "0.3",
    gap_max_height: "2.0",
    gap_height_scale: "0.001",
    gap_label_size: "NULL",
    backbone_thickness: "3.0",
    gap_width: "50",
    indel_height_scale: "0.01",
  }),
});

const PRESET_LABELS = Object.freeze({
  short_lt_5000: "Short alignment (< 5,000 bases)",
  long_ge_5000: "Long alignment (>= 5,000 bases / current defaults)",
});

function $(id) {
  return document.getElementById(id);
}

function setStatus(text, isError = false) {
  const node = $("status");
  node.textContent = text;
  node.style.color = isError ? "#9a1328" : "#243a4c";
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
    swap_roles: $("swap_roles").checked,
    params: readParams(),
  };
}

function clamp(value, minValue, maxValue) {
  return Math.min(maxValue, Math.max(minValue, value));
}

function getViewerSvg() {
  return $("viewer").querySelector("svg");
}

function parseViewBox(svg) {
  if (!svg) {
    return null;
  }
  const viewBoxAttr = svg.getAttribute("viewBox");
  if (!viewBoxAttr) {
    return null;
  }
  const parts = viewBoxAttr.split(/\s+/).map(Number);
  if (parts.length !== 4 || parts.some((v) => Number.isNaN(v))) {
    return null;
  }
  const [x, y, width, height] = parts;
  return { x, y, width, height };
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

function mapDataXToSvgX(dataX) {
  const m = state.metadata;
  if (!m) {
    return null;
  }
  const denom = m.x_data_max - m.x_data_min;
  if (denom <= 0) {
    return null;
  }
  const ratio = clamp((dataX - m.x_data_min) / denom, 0, 1);
  return m.axes_left_px + ratio * (m.axes_right_px - m.axes_left_px);
}

function updateExtractButtons() {
  const enabled = Boolean(state.token && state.selection);
  $("extract_query_btn").disabled = !enabled;
  $("extract_reference_btn").disabled = !enabled;
}

function updateCopyButton() {
  $("copy_extract_btn").disabled = $("extract_output").value.trim().length === 0;
}

async function requestProbe(xCoord) {
  const response = await fetch("/api/probe", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ token: state.token, x_coord: xCoord }),
  });
  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "Probe failed");
  }
  return data;
}

function updateProbeSummary(data) {
  $("probe_ref_base").textContent = data.reference_base;
  $("probe_query_base").textContent = data.query_base;
  $("probe_ref_pos").textContent = data.reference_local_position;
  $("probe_query_pos").textContent = data.query_local_position;
}

function resetProbeSummary() {
  $("probe_ref_base").textContent = "-";
  $("probe_query_base").textContent = "-";
  $("probe_ref_pos").textContent = "-";
  $("probe_query_pos").textContent = "-";
}

function setPresetHelper(text) {
  $("preset_helper").textContent = text;
}

function applyPreset(presetKey, helperText = null) {
  const values = PRESET_VALUES[presetKey];
  if (!values) {
    return;
  }
  for (const [id, value] of Object.entries(values)) {
    const node = $(id);
    if (node) {
      node.value = String(value);
    }
  }
  setPresetHelper(helperText || `Preset: ${PRESET_LABELS[presetKey]}`);
}

async function fetchPresetRecommendation() {
  const inputPath = $("input_path").value.trim();
  if (!inputPath) {
    return null;
  }
  const response = await fetch("/api/preset_recommendation", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      input_path: inputPath,
      swap_roles: $("swap_roles").checked,
    }),
  });
  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "Preset recommendation failed");
  }
  return data;
}

async function resolvePresetBeforeRender() {
  const mode = $("render_preset").value;
  if (mode === "short_lt_5000" || mode === "long_ge_5000") {
    applyPreset(mode);
    return;
  }
  setPresetHelper("Preset: Auto (resolving from sequence length)");
  try {
    const recommendation = await fetchPresetRecommendation();
    if (!recommendation) {
      setPresetHelper("Preset: Auto (awaiting input path)");
      return;
    }
    const recommendedPreset = recommendation.recommended_preset;
    applyPreset(
      recommendedPreset,
      `Preset: Auto selected ${PRESET_LABELS[recommendedPreset]} (basis length ${recommendation.basis_length.toLocaleString()})`,
    );
  } catch (err) {
    setPresetHelper("Preset: Auto (could not resolve recommendation)");
    setStatus(`Preset recommendation warning: ${String(err.message || err)}`, false);
  }
}

async function fetchProbeAt(xCoord, svgX) {
  const data = await requestProbe(xCoord);
  state.probeData = data;
  state.probeGlobalX = Number(data.global_x);
  state.probeSvgX = Number(svgX);
  updateProbeSummary(data);
  drawProbeLine(svgX);
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

function ensureProbeLine(svg) {
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
  return line;
}

function drawProbeLine(xValue) {
  const svg = getViewerSvg();
  if (!svg) {
    return;
  }
  const viewBox = parseViewBox(svg);
  if (!viewBox) {
    return;
  }
  const line = ensureProbeLine(svg);
  line.setAttribute("x1", String(xValue));
  line.setAttribute("x2", String(xValue));
  line.setAttribute("y1", String(viewBox.y));
  line.setAttribute("y2", String(viewBox.y + viewBox.height));
}

function ensureSelectionOverlay(svg) {
  let group = svg.querySelector("#selection-overlay");
  if (!group) {
    group = document.createElementNS("http://www.w3.org/2000/svg", "g");
    group.setAttribute("id", "selection-overlay");
    group.setAttribute("pointer-events", "none");

    const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
    rect.setAttribute("id", "selection-overlay-rect");
    rect.setAttribute("fill", "#2ea8d266");
    rect.setAttribute("stroke", "#0b87b3");
    rect.setAttribute("stroke-width", "1.2");

    const startLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
    startLine.setAttribute("id", "selection-overlay-start");
    startLine.setAttribute("stroke", "#085f7f");
    startLine.setAttribute("stroke-width", "1.2");
    startLine.setAttribute("stroke-dasharray", "4 3");

    const endLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
    endLine.setAttribute("id", "selection-overlay-end");
    endLine.setAttribute("stroke", "#085f7f");
    endLine.setAttribute("stroke-width", "1.2");
    endLine.setAttribute("stroke-dasharray", "4 3");

    group.appendChild(rect);
    group.appendChild(startLine);
    group.appendChild(endLine);
    svg.appendChild(group);
  }
  return group;
}

function hideSelectionOverlay() {
  const svg = getViewerSvg();
  if (!svg) {
    return;
  }
  const group = ensureSelectionOverlay(svg);
  group.setAttribute("display", "none");
}

function drawSelectionOverlay(startSvgX, endSvgX) {
  const svg = getViewerSvg();
  if (!svg) {
    return;
  }
  const viewBox = parseViewBox(svg);
  if (!viewBox) {
    return;
  }

  const group = ensureSelectionOverlay(svg);
  group.removeAttribute("display");

  const startX = Math.min(startSvgX, endSvgX);
  const endX = Math.max(startSvgX, endSvgX);
  const width = Math.max(0.8, endX - startX);
  const y1 = viewBox.y;
  const y2 = viewBox.y + viewBox.height;

  const rect = group.querySelector("#selection-overlay-rect");
  const startLine = group.querySelector("#selection-overlay-start");
  const endLine = group.querySelector("#selection-overlay-end");

  rect.setAttribute("x", String(startX));
  rect.setAttribute("y", String(viewBox.y));
  rect.setAttribute("width", String(width));
  rect.setAttribute("height", String(viewBox.height));

  startLine.setAttribute("x1", String(startX));
  startLine.setAttribute("x2", String(startX));
  startLine.setAttribute("y1", String(y1));
  startLine.setAttribute("y2", String(y2));

  endLine.setAttribute("x1", String(endX));
  endLine.setAttribute("x2", String(endX));
  endLine.setAttribute("y1", String(y1));
  endLine.setAttribute("y2", String(y2));
}

function renderSelectionVisuals() {
  if (state.selectionDraft) {
    drawSelectionOverlay(state.selectionDraft.startSvgX, state.selectionDraft.endSvgX);
    return;
  }
  if (state.selection) {
    drawSelectionOverlay(state.selection.startSvgX, state.selection.endSvgX);
    return;
  }
  hideSelectionOverlay();
}

function clearSelectionOutput() {
  $("extract_output").value = "";
  updateCopyButton();
}

function clearSelection() {
  state.selection = null;
  state.selectionAnchor = null;
  state.dragAnchor = null;
  state.selectionDraft = null;
  state.dragActive = false;
  state.dragMoved = false;
  state.suppressNextClick = false;
  updateExtractButtons();
  renderSelectionVisuals();
}

async function commitSelection(startDataX, endDataX) {
  const [startProbe, endProbe] = await Promise.all([requestProbe(startDataX), requestProbe(endDataX)]);

  let left = startProbe;
  let right = endProbe;
  if (right.column_index < left.column_index) {
    left = endProbe;
    right = startProbe;
  }

  const startSvgX = mapDataXToSvgX(Number(left.global_x));
  const endSvgX = mapDataXToSvgX(Number(right.global_x));
  if (startSvgX == null || endSvgX == null) {
    throw new Error("Could not map selected range to viewport coordinates.");
  }

  state.selection = {
    startColumn: Number(left.column_index),
    endColumn: Number(right.column_index),
    startGlobalX: Number(left.global_x),
    endGlobalX: Number(right.global_x),
    startSvgX,
    endSvgX,
    startQueryLocal: Number(left.query_local_position),
    endQueryLocal: Number(right.query_local_position),
    startReferenceLocal: Number(left.reference_local_position),
    endReferenceLocal: Number(right.reference_local_position),
  };

  state.selectionDraft = null;
  state.selectionAnchor = null;
  state.dragAnchor = null;
  updateExtractButtons();
  renderSelectionVisuals();
  setStatus(
    `Selected columns ${state.selection.startColumn + 1}-${state.selection.endColumn + 1}.`,
    false,
  );
}

function attachProbeInteractions() {
  const svg = getViewerSvg();
  if (!svg) {
    return;
  }

  svg.addEventListener("mousemove", (evt) => {
    if (!state.token) {
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

    if (state.dragActive && state.dragAnchor) {
      if (Math.abs(svgPoint.x - state.dragAnchor.svgX) > 1.0) {
        state.dragMoved = true;
      }
      state.selectionDraft = {
        startSvgX: state.dragAnchor.svgX,
        endSvgX: svgPoint.x,
      };
      renderSelectionVisuals();
    }
  });

  svg.addEventListener("mousedown", (evt) => {
    if (!state.token || evt.button !== 0) {
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

    state.dragActive = true;
    state.dragMoved = false;
    state.dragAnchor = { dataX, svgX: svgPoint.x };
    state.selectionDraft = {
      startSvgX: svgPoint.x,
      endSvgX: svgPoint.x,
    };
    renderSelectionVisuals();
  });

  svg.addEventListener("mouseup", async (evt) => {
    if (!state.dragActive || !state.dragAnchor || evt.button !== 0) {
      return;
    }

    const svgPoint = extractSvgPoint(evt, svg);
    state.dragActive = false;

    if (!svgPoint) {
      state.selectionDraft = null;
      state.dragAnchor = null;
      renderSelectionVisuals();
      return;
    }

    const dataX = mapSvgXToDataX(svgPoint.x);
    if (dataX == null) {
      state.selectionDraft = null;
      state.dragAnchor = null;
      renderSelectionVisuals();
      return;
    }

    if (state.dragMoved) {
      state.suppressNextClick = true;
      try {
        await commitSelection(state.dragAnchor.dataX, dataX);
      } catch (err) {
        state.selectionDraft = null;
        renderSelectionVisuals();
        setStatus(String(err.message || err), true);
      }
      state.dragAnchor = null;
      return;
    }

    state.dragAnchor = null;
    state.selectionDraft = null;
    renderSelectionVisuals();
  });

  svg.addEventListener("click", async (evt) => {
    if (!state.token) {
      return;
    }
    if (state.suppressNextClick) {
      state.suppressNextClick = false;
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

    if (!state.selectionAnchor) {
      state.selectionAnchor = { dataX, svgX: svgPoint.x };
      state.selectionDraft = { startSvgX: svgPoint.x, endSvgX: svgPoint.x };
      renderSelectionVisuals();
      setStatus("Selection start set. Click end point to complete range.", false);
      return;
    }

    try {
      await commitSelection(state.selectionAnchor.dataX, dataX);
    } catch (err) {
      state.selectionDraft = null;
      renderSelectionVisuals();
      setStatus(String(err.message || err), true);
    }
  });
  renderSelectionVisuals();
}

async function renderCurrent() {
  setStatus("Rendering...");
  await resolvePresetBeforeRender();
  const payload = buildPayload();
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
  state.probeGlobalX = null;
  state.probeSvgX = null;
  state.probeData = null;
  resetProbeSummary();
  clearSelection();
  clearSelectionOutput();

  $("viewer").innerHTML = data.svg;
  attachProbeInteractions();
  updateExtractButtons();
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

async function extractSequence(stream) {
  if (!state.token || !state.selection) {
    throw new Error("Select a range before extracting sequence.");
  }

  const response = await fetch("/api/extract_sequence", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      token: state.token,
      stream,
      start_x: state.selection.startGlobalX,
      end_x: state.selection.endGlobalX,
    }),
  });

  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "Sequence extraction failed");
  }

  $("extract_output").value = data.fasta || "";
  updateCopyButton();
  setStatus(
    `Extracted ${stream} sequence (${data.sequence_length} bp) from columns ${data.start_column_index + 1}-${data.end_column_index + 1}.`,
    false,
  );
}

async function copyExtractOutput() {
  const text = $("extract_output").value;
  if (!text.trim()) {
    return;
  }

  if (navigator.clipboard && navigator.clipboard.writeText) {
    await navigator.clipboard.writeText(text);
    setStatus("Extract output copied to clipboard.", false);
    return;
  }

  $("extract_output").focus();
  $("extract_output").select();
  document.execCommand("copy");
  setStatus("Extract output copied to clipboard.", false);
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
  document.addEventListener("keydown", (evt) => {
    if (evt.key !== "Escape") {
      return;
    }
    if (state.selectionAnchor) {
      state.selectionAnchor = null;
      state.dragAnchor = null;
      state.selectionDraft = null;
      renderSelectionVisuals();
      setStatus("Selection anchor cleared.", false);
      return;
    }
    if (state.selection) {
      clearSelection();
      setStatus("Selection cleared.", false);
    }
  });

  $("render_btn").addEventListener("click", async () => {
    try {
      await renderCurrent();
    } catch (err) {
      setStatus(String(err.message || err), true);
    }
  });

  $("render_preset").addEventListener("change", () => {
    const selected = $("render_preset").value;
    if (selected === "short_lt_5000" || selected === "long_ge_5000") {
      applyPreset(selected);
      return;
    }
    setPresetHelper("Preset: Auto (resolve on render)");
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

  $("extract_query_btn").addEventListener("click", async () => {
    try {
      await extractSequence("query");
    } catch (err) {
      setStatus(String(err.message || err), true);
    }
  });

  $("extract_reference_btn").addEventListener("click", async () => {
    try {
      await extractSequence("reference");
    } catch (err) {
      setStatus(String(err.message || err), true);
    }
  });

  $("clear_selection_btn").addEventListener("click", () => {
    clearSelection();
    setStatus("Selection cleared.", false);
  });

  $("copy_extract_btn").addEventListener("click", async () => {
    try {
      await copyExtractOutput();
    } catch (err) {
      setStatus(String(err.message || err), true);
    }
  });
}

function bootstrap() {
  setupFileBrowseButtons();
  bindEvents();
  updateExtractButtons();
  clearSelectionOutput();
  setPresetHelper("Preset: Auto (resolve on render)");
}

bootstrap();
