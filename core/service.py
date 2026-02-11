from __future__ import annotations

import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

from .compute import compute_alignment_data, construct_stream_paths
from .inputs import parse_alignment_pair, parse_annotation_file
from .params import RenderParams
from .render import (
    plot_alignment,
    render_bytes,
    render_svg_with_metadata,
    resolve_figure_width,
    write_stream_debug_tables,
)


@dataclass
class RenderSession:
    token: str
    params: RenderParams
    input_path: Path
    query_annotation_path: Optional[Path]
    reference_annotation_path: Optional[Path]
    data: object
    query_annotations: list
    reference_annotations: list
    global_x: np.ndarray
    query_x: np.ndarray
    reference_x: np.ndarray
    query_positions: np.ndarray
    reference_positions: np.ndarray
    global_extent: float
    gap_labels: list
    resolved_width: float


@dataclass
class RenderResult:
    svg: str
    token: str
    global_extent: float
    x_data_min: float
    x_data_max: float
    axes_left_px: float
    axes_right_px: float
    svg_width_px: float
    svg_height_px: float


SESSION_CACHE: Dict[str, RenderSession] = {}
MAX_SESSIONS = 12


def _trim_cache() -> None:
    while len(SESSION_CACHE) > MAX_SESSIONS:
        first_key = next(iter(SESSION_CACHE))
        del SESSION_CACHE[first_key]


def _to_optional_path(path_text: Optional[str]) -> Optional[Path]:
    if path_text is None:
        return None
    text = str(path_text).strip()
    if not text:
        return None
    return Path(text)


def prepare_session(
    *,
    input_path: Path,
    params: RenderParams,
    query_annotation_path: Optional[Path] = None,
    reference_annotation_path: Optional[Path] = None,
) -> RenderSession:
    query, reference, query_name, reference_name = parse_alignment_pair(input_path)
    query_annotations = parse_annotation_file(query_annotation_path)
    reference_annotations = parse_annotation_file(reference_annotation_path)

    data = compute_alignment_data(
        query,
        reference,
        query_name,
        reference_name,
        params.min_gap_size,
        params.block_size,
        params.min_sequence_identity,
        params.window_size,
    )

    (
        global_x,
        query_x,
        reference_x,
        query_positions,
        reference_positions,
        global_extent,
        gap_labels,
    ) = construct_stream_paths(
        data,
        data.weak_regions,
        params.backbone_gap,
        params.bump_scale,
        params.gap_max_height,
        params.gap_width,
        params.min_gap_size,
        params.gap_height_scale,
        params.indel_height_scale,
    )

    resolved_width = resolve_figure_width(params.width, global_extent)
    token = uuid.uuid4().hex
    session = RenderSession(
        token=token,
        params=params,
        input_path=input_path,
        query_annotation_path=query_annotation_path,
        reference_annotation_path=reference_annotation_path,
        data=data,
        query_annotations=query_annotations,
        reference_annotations=reference_annotations,
        global_x=global_x,
        query_x=query_x,
        reference_x=reference_x,
        query_positions=query_positions,
        reference_positions=reference_positions,
        global_extent=global_extent,
        gap_labels=gap_labels,
        resolved_width=resolved_width,
    )
    SESSION_CACHE[token] = session
    _trim_cache()
    return session


def render_alignment(
    session: RenderSession,
    viewport: Optional[Tuple[float, float]] = None,
) -> RenderResult:
    del viewport
    svg_bytes, metadata = render_svg_with_metadata(
        data=session.data,
        global_x=session.global_x,
        query_x=session.query_x,
        reference_x=session.reference_x,
        query_positions=session.query_positions,
        reference_positions=session.reference_positions,
        global_extent=session.global_extent,
        gap_labels=session.gap_labels,
        width=session.resolved_width,
        height=session.params.height,
        dpi=session.params.dpi,
        tick_interval=session.params.tick_interval,
        backbone_thickness=session.params.backbone_thickness,
        mismatch_line_width=session.params.mismatch_line_width,
        gap_label_size=session.params.gap_label_size,
        query_annotations=session.query_annotations,
        reference_annotations=session.reference_annotations,
        annotation_label_size=session.params.annotation_label_size,
        annotation_thickness=session.params.annotation_thickness,
        annotation_alpha=session.params.annotation_alpha,
        reference_annotation_color=session.params.reference_annotation_color,
        query_annotation_color=session.params.query_annotation_color,
        annotation_label_jitter=session.params.annotation_label_jitter,
        annotation_max_layers=session.params.annotation_max_layers,
        annotation_spacing=session.params.annotation_spacing,
    )
    return RenderResult(
        svg=svg_bytes.decode("utf-8"),
        token=session.token,
        global_extent=session.global_extent,
        x_data_min=float(metadata.get("x_data_min", 0.0)),
        x_data_max=float(metadata.get("x_data_max", max(1.0, session.global_extent))),
        axes_left_px=float(metadata.get("axes_left_px", 0.0)),
        axes_right_px=float(metadata.get("axes_right_px", 0.0)),
        svg_width_px=float(metadata.get("svg_width_px", 0.0)),
        svg_height_px=float(metadata.get("svg_height_px", 0.0)),
    )


def export_alignment(
    session: RenderSession,
    fmt: str,
    viewport: Optional[Tuple[float, float]] = None,
) -> bytes:
    if fmt not in {"svg", "png"}:
        raise ValueError("Export format must be 'svg' or 'png'")
    x_window = normalize_viewport(viewport, session.global_extent)
    return render_bytes(
        fmt=fmt,
        data=session.data,
        global_x=session.global_x,
        query_x=session.query_x,
        reference_x=session.reference_x,
        query_positions=session.query_positions,
        reference_positions=session.reference_positions,
        global_extent=session.global_extent,
        gap_labels=session.gap_labels,
        width=session.resolved_width,
        height=session.params.height,
        dpi=session.params.dpi,
        tick_interval=session.params.tick_interval,
        backbone_thickness=session.params.backbone_thickness,
        mismatch_line_width=session.params.mismatch_line_width,
        gap_label_size=session.params.gap_label_size,
        query_annotations=session.query_annotations,
        reference_annotations=session.reference_annotations,
        annotation_label_size=session.params.annotation_label_size,
        annotation_thickness=session.params.annotation_thickness,
        annotation_alpha=session.params.annotation_alpha,
        reference_annotation_color=session.params.reference_annotation_color,
        query_annotation_color=session.params.query_annotation_color,
        annotation_label_jitter=session.params.annotation_label_jitter,
        annotation_max_layers=session.params.annotation_max_layers,
        annotation_spacing=session.params.annotation_spacing,
        x_window=x_window,
    )


def probe_alignment(session: RenderSession, x_coord: float) -> Dict[str, object]:
    if session.global_x.size == 0:
        raise ValueError("No alignment coordinates are available")

    idx = int(np.searchsorted(session.global_x, x_coord, side="left"))
    if idx <= 0:
        nearest_idx = 0
    elif idx >= len(session.global_x):
        nearest_idx = len(session.global_x) - 1
    else:
        left = idx - 1
        right = idx
        nearest_idx = left if abs(session.global_x[left] - x_coord) <= abs(session.global_x[right] - x_coord) else right

    q_base = session.data.query[nearest_idx]
    r_base = session.data.reference[nearest_idx]
    is_mismatch = (q_base != r_base) and not session.data.is_query_gap[nearest_idx] and not session.data.is_reference_gap[nearest_idx]

    return {
        "column_index": nearest_idx,
        "global_x": float(session.global_x[nearest_idx]),
        "query_base": q_base,
        "reference_base": r_base,
        "query_local_position": int(session.data.query_local_positions[nearest_idx]),
        "reference_local_position": int(session.data.reference_local_positions[nearest_idx]),
        "is_weak": bool(session.data.is_weak[nearest_idx]),
        "is_query_gap": bool(session.data.is_query_gap[nearest_idx]),
        "is_reference_gap": bool(session.data.is_reference_gap[nearest_idx]),
        "is_mismatch": bool(is_mismatch),
    }


def normalize_viewport(
    viewport: Optional[Tuple[float, float]],
    global_extent: float,
) -> Tuple[float, float]:
    if viewport is None:
        return (0.0, max(1.0, global_extent))

    start, end = viewport
    if end < start:
        start, end = end, start

    start = max(0.0, float(start))
    hard_max = max(1.0, float(global_extent))
    end = min(hard_max, float(end))

    if end <= start:
        end = min(hard_max, start + 1.0)
        start = max(0.0, end - 1.0)

    return (start, end)


def session_from_payload(payload: Dict[str, object]) -> RenderSession:
    input_text = str(payload.get("input_path", "")).strip()
    if not input_text:
        raise ValueError("input_path is required")

    params = RenderParams.from_payload(payload.get("params", {}))
    session = prepare_session(
        input_path=Path(input_text),
        params=params,
        query_annotation_path=_to_optional_path(payload.get("query_annotation_path")),
        reference_annotation_path=_to_optional_path(payload.get("reference_annotation_path")),
    )
    return session


def get_session(token: str) -> RenderSession:
    try:
        return SESSION_CACHE[token]
    except KeyError as exc:
        raise ValueError("Unknown or expired render token") from exc


def write_debug_tables_for_cli(session: RenderSession, output: Path) -> None:
    write_stream_debug_tables(
        output,
        session.data,
        session.global_x,
        session.query_x,
        session.reference_x,
        session.query_positions,
        session.reference_positions,
    )


def export_to_file_for_cli(session: RenderSession, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    plot_alignment(
        data=session.data,
        global_x=session.global_x,
        query_x=session.query_x,
        reference_x=session.reference_x,
        query_positions=session.query_positions,
        reference_positions=session.reference_positions,
        global_extent=session.global_extent,
        gap_labels=session.gap_labels,
        width=session.resolved_width,
        height=session.params.height,
        dpi=session.params.dpi,
        output=output,
        tick_interval=session.params.tick_interval,
        backbone_thickness=session.params.backbone_thickness,
        mismatch_line_width=session.params.mismatch_line_width,
        gap_label_size=session.params.gap_label_size,
        query_annotations=session.query_annotations,
        reference_annotations=session.reference_annotations,
        annotation_label_size=session.params.annotation_label_size,
        annotation_thickness=session.params.annotation_thickness,
        annotation_alpha=session.params.annotation_alpha,
        reference_annotation_color=session.params.reference_annotation_color,
        query_annotation_color=session.params.query_annotation_color,
        annotation_label_jitter=session.params.annotation_label_jitter,
        annotation_max_layers=session.params.annotation_max_layers,
        annotation_spacing=session.params.annotation_spacing,
        x_window=None,
    )
