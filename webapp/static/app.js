const state = {
  token: null,
  globalExtent: 0,
  metadata: null,
  featureTags: [],
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
  genome: {
    uploadId: null,
    dotplot: null,
    dotplotJobId: null,
    alignJobId: null,
    selectedBlocks: new Set(),
    segmentCache: [],
    hoveredBlockId: null,
    canvasBoundHandlers: false,
  },
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
  if (!node) {
    return;
  }
  node.textContent = text;
  node.style.color = isError ? "#9a1328" : "#243a4c";
}

function setGenomeStatus(text, isError = false) {
  const node = $("genome_status");
  if (!node) {
    return;
  }
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

function setProbeFeatureNote(text, isDuplication = false) {
  const node = $("probe_feature_note");
  if (!node) {
    return;
  }
  node.textContent = text;
  node.classList.toggle("dup-active", Boolean(isDuplication));
}

function updateProbeFeaturePointer(globalX) {
  if (!Number.isFinite(globalX)) {
    setProbeFeatureNote("-", false);
    return;
  }
  const duplicationTags = state.featureTags.filter((tag) => {
    if (!tag || String(tag.kind || "").toLowerCase() !== "duplication") {
      return false;
    }
    const start = Number(tag.start_x);
    const end = Number(tag.end_x);
    if (!Number.isFinite(start) || !Number.isFinite(end)) {
      return false;
    }
    const low = Math.min(start, end);
    const high = Math.max(start, end);
    return globalX >= low && globalX <= high;
  });
  if (!duplicationTags.length) {
    setProbeFeatureNote("-", false);
    return;
  }
  const first = duplicationTags[0];
  const pointerText = String(first.pointer_text || "").trim() || "Duplicated sequence";
  const extra = duplicationTags.length > 1 ? ` (+${duplicationTags.length - 1} more)` : "";
  setProbeFeatureNote(`${pointerText}${extra}`, true);
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
  updateProbeFeaturePointer(Number(data.global_x));
}

function resetProbeSummary() {
  $("probe_ref_base").textContent = "-";
  $("probe_query_base").textContent = "-";
  $("probe_ref_pos").textContent = "-";
  $("probe_query_pos").textContent = "-";
  setProbeFeatureNote("-", false);
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
  setStatus(`Selected columns ${state.selection.startColumn + 1}-${state.selection.endColumn + 1}.`, false);
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

function applyRenderPayload(data, statusText = null) {
  state.token = data.token;
  state.globalExtent = Number(data.global_extent) || 0;
  state.featureTags = Array.isArray(data.feature_tags) ? data.feature_tags : [];
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

  const viewerNode = $("viewer");
  viewerNode.innerHTML = data.svg;
  const svgNode = viewerNode.querySelector("svg");
  if (svgNode) {
    const metadataWidth = Number(state.metadata.svg_width_px) || 0;
    if (metadataWidth > 0) {
      svgNode.style.width = `${metadataWidth}px`;
      svgNode.style.maxWidth = "none";
      svgNode.style.height = "auto";
    } else {
      svgNode.style.width = "max-content";
      svgNode.style.maxWidth = "none";
      svgNode.style.height = "auto";
    }
  }
  attachProbeInteractions();
  updateExtractButtons();
  const invCount = Array.isArray(data.inversion_regions) ? data.inversion_regions.length : 0;
  const dupCount = state.featureTags.filter((tag) => String(tag?.kind || "").toLowerCase() === "duplication").length;
  const featureParts = [];
  if (invCount > 0) {
    featureParts.push(`Inversions: ${invCount}`);
  }
  if (dupCount > 0) {
    featureParts.push(`Duplications: ${dupCount}`);
  }
  const featureSuffix = featureParts.length ? ` | ${featureParts.join(" | ")}` : "";
  if (statusText) {
    const details = featureParts.length ? ` (${featureParts.join(", ")})` : "";
    setStatus(`${statusText}${details}`, false);
    return;
  }
  setStatus(`Rendered. Query: ${data.query_name} | Reference: ${data.reference_name}${featureSuffix}`, false);
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
  applyRenderPayload(data);
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

function setActiveTab(tabName) {
  document.querySelectorAll(".tab-btn").forEach((btn) => {
    btn.classList.toggle("active", btn.getAttribute("data-tab-target") === tabName);
  });
  document.querySelectorAll(".tab-panel").forEach((panel) => {
    panel.classList.toggle("active", panel.id === `tab-${tabName}`);
  });
  if (tabName === "genome") {
    requestAnimationFrame(() => {
      redrawGenomeDotplot();
    });
  }
}

const DOTPLOT_HIT_RADIUS_PX = 6;

function getGenomeBlocks() {
  return Array.isArray(state.genome.dotplot?.dotplot?.blocks) ? state.genome.dotplot.dotplot.blocks : [];
}

function updateGenomeAlignButton() {
  $("genome_align_btn").disabled = state.genome.selectedBlocks.size === 0;
}

function syncGenomeSelectionUI() {
  renderBlockList(getGenomeBlocks());
  redrawGenomeDotplot();
  updateGenomeAlignButton();
}

function ensureDotplotCanvasSize(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const cssWidth = Math.max(360, Math.floor(canvas.clientWidth || 1000));
  const cssHeight = Math.max(420, Math.floor(canvas.clientHeight || 620));
  const targetWidth = Math.floor(cssWidth * dpr);
  const targetHeight = Math.floor(cssHeight * dpr);

  if (canvas.width !== targetWidth || canvas.height !== targetHeight) {
    canvas.width = targetWidth;
    canvas.height = targetHeight;
  }

  const ctx = canvas.getContext("2d");
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.scale(dpr, dpr);
  return { ctx, width: cssWidth, height: cssHeight };
}

function getDotplotHoverCandidate(evt) {
  const canvas = $("genome_dotplot_canvas");
  const rect = canvas.getBoundingClientRect();
  const x = evt.clientX - rect.left;
  const y = evt.clientY - rect.top;
  let best = null;
  let bestDistance = Infinity;

  for (const segment of state.genome.segmentCache) {
    const ax = segment.x1;
    const ay = segment.y1;
    const bx = segment.x2;
    const by = segment.y2;
    const abx = bx - ax;
    const aby = by - ay;
    const apx = x - ax;
    const apy = y - ay;
    const ab2 = abx * abx + aby * aby;
    const t = ab2 <= 1e-9 ? 0 : clamp((apx * abx + apy * aby) / ab2, 0, 1);
    const cx = ax + t * abx;
    const cy = ay + t * aby;
    const dx = x - cx;
    const dy = y - cy;
    const dist = Math.hypot(dx, dy);
    if (dist < bestDistance) {
      bestDistance = dist;
      best = segment;
    }
  }

  if (best && bestDistance <= DOTPLOT_HIT_RADIUS_PX) {
    return best;
  }
  return null;
}

function redrawGenomeDotplot() {
  drawDotplot(state.genome.dotplot);
}

function bindDotplotInteractions() {
  const canvas = $("genome_dotplot_canvas");
  if (!canvas || state.genome.canvasBoundHandlers) {
    return;
  }
  state.genome.canvasBoundHandlers = true;

  canvas.addEventListener("mousemove", (evt) => {
    if (!state.genome.dotplot) {
      return;
    }
    const candidate = getDotplotHoverCandidate(evt);
    const hoveredBlockId = candidate ? candidate.blockId : null;
    canvas.style.cursor = candidate ? "pointer" : "crosshair";
    if (hoveredBlockId === state.genome.hoveredBlockId) {
      return;
    }
    state.genome.hoveredBlockId = hoveredBlockId;
    redrawGenomeDotplot();
  });

  canvas.addEventListener("mouseleave", () => {
    if (state.genome.hoveredBlockId == null) {
      return;
    }
    state.genome.hoveredBlockId = null;
    canvas.style.cursor = "crosshair";
    redrawGenomeDotplot();
  });

  canvas.addEventListener("click", (evt) => {
    if (!state.genome.dotplot) {
      return;
    }
    const candidate = getDotplotHoverCandidate(evt);
    if (!candidate) {
      return;
    }
    const blockId = candidate.blockId;
    if (evt.shiftKey) {
      if (state.genome.selectedBlocks.has(blockId)) {
        state.genome.selectedBlocks.delete(blockId);
      } else {
        state.genome.selectedBlocks.add(blockId);
      }
    } else {
      state.genome.selectedBlocks = new Set([blockId]);
    }
    syncGenomeSelectionUI();
  });
}

function drawDotplot(dotplotPayload) {
  const canvas = $("genome_dotplot_canvas");
  const { ctx, width, height } = ensureDotplotCanvasSize(canvas);
  const pad = 48;

  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  ctx.strokeStyle = "#92a8bf";
  ctx.lineWidth = 1;
  ctx.strokeRect(pad, pad, width - 2 * pad, height - 2 * pad);

  if (!dotplotPayload) {
    state.genome.segmentCache = [];
    return;
  }

  const qLen = Number(dotplotPayload.query_length) || 1;
  const rLen = Number(dotplotPayload.reference_length) || 1;
  const blocks = getGenomeBlocks();
  const rangeWidth = Math.max(1, width - 2 * pad);
  const rangeHeight = Math.max(1, height - 2 * pad);

  ctx.fillStyle = "#2a4155";
  ctx.font = "12px IBM Plex Sans, sans-serif";
  ctx.fillText("Reference", 8, pad - 12);
  ctx.save();
  ctx.translate(width - 10, height - 8);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("Query", 0, 0);
  ctx.restore();

  const segments = [];
  for (const block of blocks) {
    const x1 = pad + (Number(block.r_start) / rLen) * rangeWidth;
    const x2 = pad + (Number(block.r_end) / rLen) * rangeWidth;
    const y1 = height - pad - (Number(block.q_start) / qLen) * rangeHeight;
    const y2 = height - pad - (Number(block.q_end) / qLen) * rangeHeight;
    const selected = state.genome.selectedBlocks.has(block.block_id);
    const hovered = state.genome.hoveredBlockId === block.block_id;
    const reverse = block.orientation === "reverse";

    let strokeStyle = reverse ? "rgba(191, 45, 55, 0.32)" : "rgba(35, 110, 171, 0.32)";
    let lineWidth = 1.0;
    if (selected) {
      strokeStyle = reverse ? "rgba(191, 45, 55, 0.9)" : "rgba(24, 111, 187, 0.9)";
      lineWidth = 2.2;
    }
    if (hovered) {
      strokeStyle = reverse ? "rgba(147, 23, 31, 0.95)" : "rgba(11, 87, 153, 0.95)";
      lineWidth += 1.4;
    }

    ctx.beginPath();
    ctx.strokeStyle = strokeStyle;
    ctx.lineWidth = lineWidth;
    ctx.moveTo(x1, y1);
    ctx.lineTo(x2, y2);
    ctx.stroke();

    if (selected || hovered) {
      ctx.fillStyle = reverse ? "#932129" : "#0b5799";
      ctx.beginPath();
      ctx.arc(x1, y1, 1.9, 0, Math.PI * 2);
      ctx.fill();
      ctx.beginPath();
      ctx.arc(x2, y2, 1.9, 0, Math.PI * 2);
      ctx.fill();
    }

    segments.push({
      blockId: block.block_id,
      x1,
      y1,
      x2,
      y2,
    });
  }

  state.genome.segmentCache = segments;
}

function renderBlockList(blocks) {
  const node = $("genome_block_list");
  node.innerHTML = "";

  if (!Array.isArray(blocks) || blocks.length === 0) {
    node.textContent = "No blocks found. You can still run whole-sequence alignment by selecting all once available.";
    return;
  }

  const frag = document.createDocumentFragment();
  for (const block of blocks) {
    const selected = state.genome.selectedBlocks.has(block.block_id);
    const row = document.createElement("div");
    row.className = `block-row ${block.orientation === "reverse" ? "reverse" : "forward"}${selected ? " selected" : ""}`;

    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.checked = selected;
    checkbox.addEventListener("change", () => {
      if (checkbox.checked) {
        state.genome.selectedBlocks.add(block.block_id);
      } else {
        state.genome.selectedBlocks.delete(block.block_id);
      }
      syncGenomeSelectionUI();
    });

    const text = document.createElement("div");
    text.textContent = `${block.block_id} | q:${block.q_start}-${block.q_end} | r:${block.r_start}-${block.r_end} | id:${Number(block.score).toFixed(2)}%`;

    const orient = document.createElement("span");
    orient.className = "block-orientation";
    orient.textContent = block.orientation === "reverse" ? "REV" : "FWD";

    row.appendChild(checkbox);
    row.appendChild(text);
    row.appendChild(orient);
    frag.appendChild(row);
  }
  node.appendChild(frag);
}

function readNucmerOptions() {
  const matchMode = $("genome_nucmer_match_mode").value || "maxmatch";
  const deltaFilterMode = $("genome_delta_filter_mode").value || "none";

  const parseOptionalPositiveInt = (raw, label) => {
    const text = String(raw || "").trim();
    if (!text) {
      return null;
    }
    const parsed = Number.parseInt(text, 10);
    if (!Number.isFinite(parsed) || parsed <= 0) {
      throw new Error(`${label} must be a positive integer.`);
    }
    return parsed;
  };

  const parseOptionalPositiveFloat = (raw, label) => {
    const text = String(raw || "").trim();
    if (!text) {
      return null;
    }
    const parsed = Number.parseFloat(text);
    if (!Number.isFinite(parsed) || parsed <= 0) {
      throw new Error(`${label} must be a positive number.`);
    }
    return parsed;
  };

  const mincluster = parseOptionalPositiveInt($("genome_nucmer_mincluster").value, "Min cluster (-c)");
  const diagfactor = parseOptionalPositiveFloat($("genome_nucmer_diagfactor").value, "Diag factor (-D)");
  const breaklen = parseOptionalPositiveInt($("genome_nucmer_breaklen").value, "Break len (-b)");

  const options = { match_mode: matchMode, delta_filter_mode: deltaFilterMode };
  if (mincluster != null) {
    options.mincluster = mincluster;
  }
  if (diagfactor != null) {
    options.diagfactor = diagfactor;
  }
  if (breaklen != null) {
    options.breaklen = breaklen;
  }
  return options;
}

async function pollJob(jobId, onProgress) {
  while (true) {
    const response = await fetch(`/api/genome/jobs/${jobId}`);
    const data = await response.json();
    if (!response.ok) {
      throw new Error(data.error || "Job polling failed");
    }

    if (typeof onProgress === "function") {
      onProgress(data);
    }

    if (data.status === "done") {
      return data;
    }
    if (data.status === "failed" || data.status === "cancelled") {
      throw new Error(data.error || data.message || `Job ended with status ${data.status}`);
    }

    await new Promise((resolve) => setTimeout(resolve, 800));
  }
}

async function uploadGenomeFastas() {
  const qFile = $("genome_query_fasta").files?.[0];
  const rFile = $("genome_reference_fasta").files?.[0];
  if (!qFile || !rFile) {
    throw new Error("Choose both query and reference FASTA files.");
  }

  const form = new FormData();
  form.append("query_fasta", qFile);
  form.append("reference_fasta", rFile);

  const response = await fetch("/api/genome/upload", {
    method: "POST",
    body: form,
  });
  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "Genome upload failed");
  }

  state.genome.uploadId = data.upload_id;
  state.genome.dotplot = null;
  state.genome.dotplotJobId = null;
  state.genome.alignJobId = null;
  state.genome.selectedBlocks = new Set();
  state.genome.segmentCache = [];
  state.genome.hoveredBlockId = null;

  $("genome_dotplot_btn").disabled = false;
  $("genome_align_btn").disabled = true;
  $("genome_send_viewer_btn").disabled = true;
  $("genome_dotplot_meta").textContent = "Upload and run dotplot to populate this panel.";
  renderBlockList([]);
  redrawGenomeDotplot();

  const warnings = Array.isArray(data.warnings) && data.warnings.length > 0 ? ` Warnings: ${data.warnings.join(" | ")}` : "";
  setGenomeStatus(
    `Uploaded. Query: ${data.query_name} (${Number(data.query_length).toLocaleString()} bp), Reference: ${data.reference_name} (${Number(data.reference_length).toLocaleString()} bp).${warnings}`,
    false,
  );
}

async function runDotplot() {
  if (!state.genome.uploadId) {
    throw new Error("Upload FASTA files first.");
  }
  const nucmerOptions = readNucmerOptions();

  const response = await fetch("/api/genome/dotplot/start", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ upload_id: state.genome.uploadId, nucmer_options: nucmerOptions }),
  });
  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "Failed to start dotplot job");
  }

  state.genome.dotplotJobId = data.job_id;
  setGenomeStatus("Dotplot started.", false);

  await pollJob(data.job_id, (job) => {
    const pct = Math.round(Number(job.progress || 0) * 100);
    let extra = "";
    if (job.stage === "done" && Array.isArray(job.result?.nucmer_command)) {
      extra = ` | cmd: ${job.result.nucmer_command.join(" ")}`;
    }
    setGenomeStatus(`Dotplot job: ${job.stage} (${pct}%) ${job.message || ""}${extra}`.trim(), false);
  });

  const dotplotResp = await fetch(`/api/genome/dotplot/${state.genome.uploadId}`);
  const dotplotData = await dotplotResp.json();
  if (!dotplotResp.ok) {
    throw new Error(dotplotData.error || "Failed to load dotplot");
  }

  state.genome.dotplot = dotplotData;
  const blocks = Array.isArray(dotplotData.dotplot?.blocks) ? dotplotData.dotplot.blocks : [];
  state.genome.selectedBlocks = new Set(blocks.map((item) => item.block_id));
  state.genome.hoveredBlockId = null;
  syncGenomeSelectionUI();

  $("genome_dotplot_meta").textContent = `Blocks: ${blocks.length.toLocaleString()} | Points: ${(dotplotData.dotplot?.points || []).length.toLocaleString()} | Query: ${Number(dotplotData.query_length).toLocaleString()} bp | Reference: ${Number(dotplotData.reference_length).toLocaleString()} bp`;

  const jobStatusResp = await fetch(`/api/genome/jobs/${state.genome.dotplotJobId}`);
  const jobStatus = await jobStatusResp.json();
  let statusMsg = "Dotplot complete. Select blocks and run MAFFT.";
  if (jobStatusResp.ok && Array.isArray(jobStatus.result?.nucmer_command)) {
    statusMsg += ` Command: ${jobStatus.result.nucmer_command.join(" ")}`;
  }
  if (blocks.length === 0) {
    statusMsg += " No blocks found. Try switching match mode or lowering mincluster.";
  }
  setGenomeStatus(statusMsg, false);
}

async function runGenomeAlignment() {
  if (!state.genome.uploadId || !state.genome.dotplot) {
    throw new Error("Run the dotplot first.");
  }

  const blocks = Array.isArray(state.genome.dotplot.dotplot?.blocks) ? state.genome.dotplot.dotplot.blocks : [];
  const selectedPayload = blocks.map((block) => ({
    block_id: block.block_id,
    include: state.genome.selectedBlocks.has(block.block_id),
  }));
  const includeIntervals = $("genome_include_intervals").checked;

  const response = await fetch("/api/genome/align/start", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      upload_id: state.genome.uploadId,
      selected_blocks: selectedPayload,
      align_options: {
        include_inter_block_intervals: includeIntervals,
      },
    }),
  });
  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "Failed to start alignment job");
  }

  state.genome.alignJobId = data.job_id;
  $("genome_send_viewer_btn").disabled = true;

  await pollJob(data.job_id, (job) => {
    const pct = Math.round(Number(job.progress || 0) * 100);
    setGenomeStatus(`Alignment job: ${job.stage} (${pct}%) ${job.message || ""}`.trim(), false);
  });

  $("genome_send_viewer_btn").disabled = false;
  setGenomeStatus("Alignment complete. Send result to Viewer when ready.", false);
}

async function sendGenomeResultToViewer() {
  if (!state.genome.alignJobId) {
    throw new Error("Run MAFFT alignment first.");
  }

  await resolvePresetBeforeRender();
  const response = await fetch("/api/genome/send_to_viewer", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      alignment_job_id: state.genome.alignJobId,
      params: readParams(),
      query_annotation_path: $("query_annotation_path").value.trim(),
      reference_annotation_path: $("reference_annotation_path").value.trim(),
      swap_roles: $("swap_roles").checked,
    }),
  });
  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "Failed to send genome alignment to viewer");
  }

  applyRenderPayload(data, "Genome alignment loaded into viewer.");
  setActiveTab("viewer");
  setGenomeStatus("Sent to viewer. Inversion and duplication regions are shaded and tagged.", false);
}

function bindViewerEvents() {
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
      setStatus("Exported SVG.", false);
    } catch (err) {
      setStatus(String(err.message || err), true);
    }
  });

  $("export_png").addEventListener("click", async () => {
    try {
      await exportFigure("png");
      setStatus("Exported PNG.", false);
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

function bindGenomeEvents() {
  bindDotplotInteractions();

  let resizeTimer = null;
  window.addEventListener("resize", () => {
    if (!state.genome.dotplot) {
      return;
    }
    if (resizeTimer) {
      clearTimeout(resizeTimer);
    }
    resizeTimer = setTimeout(() => {
      redrawGenomeDotplot();
      resizeTimer = null;
    }, 100);
  });

  $("genome_upload_btn").addEventListener("click", async () => {
    try {
      await uploadGenomeFastas();
    } catch (err) {
      setGenomeStatus(String(err.message || err), true);
    }
  });

  $("genome_dotplot_btn").addEventListener("click", async () => {
    try {
      await runDotplot();
    } catch (err) {
      setGenomeStatus(String(err.message || err), true);
    }
  });

  $("genome_align_btn").addEventListener("click", async () => {
    try {
      await runGenomeAlignment();
    } catch (err) {
      setGenomeStatus(String(err.message || err), true);
    }
  });

  $("genome_send_viewer_btn").addEventListener("click", async () => {
    try {
      await sendGenomeResultToViewer();
    } catch (err) {
      setGenomeStatus(String(err.message || err), true);
    }
  });
}

function bindTabEvents() {
  document.querySelectorAll(".tab-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      const target = btn.getAttribute("data-tab-target");
      setActiveTab(target);
    });
  });
}

function bootstrap() {
  setupFileBrowseButtons();
  bindTabEvents();
  bindViewerEvents();
  bindGenomeEvents();
  updateExtractButtons();
  clearSelectionOutput();
  setPresetHelper("Preset: Auto (resolve on render)");
  $("genome_dotplot_btn").disabled = true;
  $("genome_align_btn").disabled = true;
  $("genome_send_viewer_btn").disabled = true;
  setActiveTab("viewer");
}

bootstrap();
