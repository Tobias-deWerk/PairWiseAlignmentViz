from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Optional, Sequence, Tuple


def resolve_figure_width(width_arg, global_extent: float) -> float:
    from align_viz import resolve_figure_width as _resolve_figure_width

    return _resolve_figure_width(width_arg, global_extent)


def write_stream_debug_tables(
    output: Path,
    data,
    global_x,
    query_x,
    reference_x,
    query_positions,
    reference_positions,
) -> None:
    from align_viz import write_stream_debug_tables as _write_stream_debug_tables

    _write_stream_debug_tables(
        output,
        data,
        global_x,
        query_x,
        reference_x,
        query_positions,
        reference_positions,
    )


def plot_alignment(
    data,
    global_x,
    query_x,
    reference_x,
    query_positions,
    reference_positions,
    global_extent: float,
    gap_labels,
    width: float,
    height: float,
    dpi: int,
    output: Path,
    tick_interval: int,
    backbone_thickness: float,
    mismatch_line_width: float,
    gap_label_size,
    query_annotations,
    reference_annotations,
    annotation_label_size,
    annotation_thickness: float,
    annotation_alpha: float,
    reference_annotation_color: str,
    query_annotation_color: str,
    annotation_label_jitter: float,
    annotation_max_layers: int,
    annotation_spacing: float,
    x_window: Optional[Tuple[float, float]] = None,
) -> None:
    from align_viz import plot_alignment as _plot_alignment

    _plot_alignment(
        data,
        global_x,
        query_x,
        reference_x,
        query_positions,
        reference_positions,
        global_extent,
        gap_labels,
        width,
        height,
        dpi,
        output,
        tick_interval,
        backbone_thickness,
        mismatch_line_width,
        gap_label_size,
        query_annotations,
        reference_annotations,
        annotation_label_size,
        annotation_thickness,
        annotation_alpha,
        reference_annotation_color,
        query_annotation_color,
        annotation_label_jitter,
        annotation_max_layers,
        annotation_spacing,
        x_window,
    )


def render_bytes(
    *,
    fmt: str,
    data,
    global_x,
    query_x,
    reference_x,
    query_positions,
    reference_positions,
    global_extent: float,
    gap_labels,
    width: float,
    height: float,
    dpi: int,
    tick_interval: int,
    backbone_thickness: float,
    mismatch_line_width: float,
    gap_label_size,
    query_annotations,
    reference_annotations,
    annotation_label_size,
    annotation_thickness: float,
    annotation_alpha: float,
    reference_annotation_color: str,
    query_annotation_color: str,
    annotation_label_jitter: float,
    annotation_max_layers: int,
    annotation_spacing: float,
    x_window: Optional[Tuple[float, float]] = None,
) -> bytes:
    suffix = ".svg" if fmt == "svg" else ".png"
    with NamedTemporaryFile(suffix=suffix, delete=True) as handle:
        plot_alignment(
            data,
            global_x,
            query_x,
            reference_x,
            query_positions,
            reference_positions,
            global_extent,
            gap_labels,
            width,
            height,
            dpi,
            Path(handle.name),
            tick_interval,
            backbone_thickness,
            mismatch_line_width,
            gap_label_size,
            query_annotations,
            reference_annotations,
            annotation_label_size,
            annotation_thickness,
            annotation_alpha,
            reference_annotation_color,
            query_annotation_color,
            annotation_label_jitter,
            annotation_max_layers,
            annotation_spacing,
            x_window,
        )
        handle.seek(0)
        return handle.read()
