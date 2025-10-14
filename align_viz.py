#!/usr/bin/env python3
"""Command line visualization tool for pairwise DNA alignment outputs.

The tool expects a FASTA file that contains two aligned sequences. The
visualization renders two horizontal streams (query on top, reference on
bottom) and decorates them with glyphs describing indels, weakly aligned
regions, and mismatch ladder rungs.
"""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib import transforms
from matplotlib.path import Path as MplPath


NUCLEOTIDE_COLORS: Dict[str, str] = {
    "A": "#4daf4a",  # green
    "C": "#377eb8",  # blue
    "G": "#000000",  # black
    "T": "#e41a1c",  # red
    "U": "#e41a1c",  # treat U like T
}
@dataclass
class GapRun:
    start: int
    end: int
    length: int
    stream: str  # "query" or "reference"


@dataclass
class WeakRegion:
    start: int
    end: int
    identity: float


@dataclass
class GapCircle:
    run: GapRun
    anchor_x: float
    anchor_y: float
    center_x: float
    center_y: float
    radius: float
    jitter: float


@dataclass
class AlignmentData:
    query: str
    reference: str
    matches: List[bool]
    is_query_gap: List[bool]
    is_reference_gap: List[bool]
    query_local_positions: List[int]
    reference_local_positions: List[int]
    is_weak: List[bool]
    gap_runs: List[GapRun]
    weak_regions: List[WeakRegion]


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize pairwise DNA sequence alignments with loop/bump glyphs."
    )
    parser.add_argument("input", type=Path, help="FASTA file containing the alignment")
    parser.add_argument("width", type=float, help="Figure width (inches)")
    parser.add_argument("height", type=float, help="Figure height (inches)")
    parser.add_argument("dpi", type=int, help="Figure resolution in dots per inch")
    parser.add_argument(
        "output",
        type=Path,
        help="Output file (.svg, .pdf, .png, .jpg, etc.)",
    )
    parser.add_argument(
        "--min-gap-size",
        type=int,
        default=10,
        help="Gaps shorter than this threshold are rendered as beaks instead of loops",
    )
    parser.add_argument(
        "--block-size",
        type=int,
        default=100,
        help="Block length (in columns) that separates the per-block and sliding window heuristics",
    )
    parser.add_argument(
        "--min-sequence-identity",
        type=float,
        default=0.7,
        help="Identity threshold used to flag weakly aligned sequence segments",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=20,
        help="Window size for scanning long regions for weak alignment",
    )
    parser.add_argument(
        "--tick-interval",
        type=int,
        default=10_000,
        help="Spacing for local coordinate tick marks (set to 0 to disable)",
    )
    parser.add_argument(
        "--backbone-gap",
        type=float,
        default=1.0,
        help="Vertical distance between the query and reference backbones",
    )
    parser.add_argument(
        "--backbone-thickness",
        type=float,
        default=2.0,
        help="Line width used when drawing the query and reference backbones",
    )
    parser.add_argument(
        "--bump-scale",
        type=float,
        default=1.0,
        help="Multiplier applied to weak-alignment bump heights",
    )
    parser.add_argument(
        "--gap-circle-scale",
        type=float,
        default=0.02,
        help=(
            "Scale factor that converts gap length (bp) to circle diameter units; "
            "center height tracks half the gap length after scaling"
        ),
    )
    parser.add_argument(
        "--mismatch-line-width",
        type=float,
        default=1.2,
        help="Line width used for mismatch ladder rungs",
    )
    return parser.parse_args(argv)


def parse_fasta_pair(path: Path) -> Tuple[str, str]:
    sequences: List[List[str]] = []
    current_seq: List[str] | None = None

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_seq = []
                sequences.append(current_seq)
            else:
                if current_seq is None:
                    raise ValueError("FASTA file must start with a header line beginning with '>'")
                current_seq.append(line.strip())

    if len(sequences) != 2:
        raise ValueError(
            f"Expected exactly two sequences in the FASTA file, found {len(sequences)}"
        )

    query = "".join(sequences[1]).upper()
    reference = "".join(sequences[0]).upper()

    if len(query) != len(reference):
        raise ValueError(
            "Aligned sequences must have the same length (including gap characters)."
        )

    return query, reference


def identify_gap_runs(query: str, reference: str) -> List[GapRun]:
    runs: List[GapRun] = []
    n = len(query)

    def scan(sequence: str, stream_name: str) -> None:
        idx = 0
        while idx < n:
            if sequence[idx] == "-" and not (query[idx] == "-" and reference[idx] == "-"):
                start = idx
                while idx < n and sequence[idx] == "-":
                    idx += 1
                end = idx - 1
                runs.append(GapRun(start, end, end - start + 1, stream_name))
            else:
                idx += 1

    scan(query, "query")
    scan(reference, "reference")
    return sorted(runs, key=lambda r: r.start)


def compute_alignment_data(
    query: str,
    reference: str,
    min_gap_size: int,
    block_size: int,
    min_sequence_identity: float,
    window_size: int,
) -> AlignmentData:
    n = len(query)
    matches: List[bool] = []
    is_query_gap: List[bool] = []
    is_reference_gap: List[bool] = []
    query_local_positions: List[int] = []
    reference_local_positions: List[int] = []
    query_local = 0
    reference_local = 0

    for q_char, r_char in zip(query, reference):
        is_gap_q = q_char == "-"
        is_gap_r = r_char == "-"
        is_match = (q_char == r_char) and not (is_gap_q or is_gap_r)
        matches.append(is_match)
        is_query_gap.append(is_gap_q)
        is_reference_gap.append(is_gap_r)
        if not is_gap_q:
            query_local += 1
        query_local_positions.append(query_local)
        if not is_gap_r:
            reference_local += 1
        reference_local_positions.append(reference_local)

    gap_runs = identify_gap_runs(query, reference)
    is_weak = identify_weak_regions(
        query,
        reference,
        matches,
        is_query_gap,
        is_reference_gap,
        gap_runs,
        min_gap_size,
        block_size,
        min_sequence_identity,
        window_size,
    )
    weak_regions = build_weak_region_list(is_weak, matches, min_sequence_identity)

    return AlignmentData(
        query=query,
        reference=reference,
        matches=matches,
        is_query_gap=is_query_gap,
        is_reference_gap=is_reference_gap,
        query_local_positions=query_local_positions,
        reference_local_positions=reference_local_positions,
        is_weak=is_weak,
        gap_runs=gap_runs,
        weak_regions=weak_regions,
    )


def identify_weak_regions(
    query: str,
    reference: str,
    matches: Sequence[bool],
    is_query_gap: Sequence[bool],
    is_reference_gap: Sequence[bool],
    gap_runs: Sequence[GapRun],
    min_gap_size: int,
    block_size: int,
    min_sequence_identity: float,
    window_size: int,
) -> List[bool]:
    n = len(query)
    is_weak = [False] * n

    gap_end_lengths: Dict[int, int] = {}
    gap_start_lengths: Dict[int, int] = {}
    for run in gap_runs:
        gap_start_lengths[run.start] = max(gap_start_lengths.get(run.start, 0), run.length)
        gap_end_lengths[run.end] = max(gap_end_lengths.get(run.end, 0), run.length)

    segments: List[Tuple[int, int]] = []
    idx = 0
    while idx < n:
        if not is_query_gap[idx] and not is_reference_gap[idx]:
            start = idx
            while idx < n and not is_query_gap[idx] and not is_reference_gap[idx]:
                idx += 1
            segments.append((start, idx - 1))
        else:
            idx += 1

    for start, end in segments:
        length = end - start + 1
        if length <= 0:
            continue
        segment_matches = matches[start : end + 1]
        match_count = sum(1 for m in segment_matches if m)
        identity = match_count / length

        left_gap_len = gap_end_lengths.get(start - 1, 0)
        right_gap_len = gap_start_lengths.get(end + 1, 0)

        if (
            length < block_size
            and left_gap_len >= min_gap_size
            and right_gap_len >= min_gap_size
            and identity < min_sequence_identity
        ):
            for idx in range(start, end + 1):
                is_weak[idx] = True
            continue

        if length < window_size:
            if identity < min_sequence_identity:
                for idx in range(start, end + 1):
                    is_weak[idx] = True
            continue

        window_identity_scan(
            start,
            end,
            segment_matches,
            matches,
            is_weak,
            min_sequence_identity,
            window_size,
        )

    return is_weak


def window_identity_scan(
    segment_start: int,
    segment_end: int,
    segment_matches: Sequence[bool],
    matches: Sequence[bool],
    is_weak: List[bool],
    min_sequence_identity: float,
    window_size: int,
) -> None:
    length = segment_end - segment_start + 1
    # Precompute cumulative counts for fast window identity calculation
    cumulative: List[int] = [0]
    for m in segment_matches:
        cumulative.append(cumulative[-1] + (1 if m else 0))

    for offset in range(0, length - window_size + 1):
        window_match_count = cumulative[offset + window_size] - cumulative[offset]
        identity = window_match_count / window_size
        if identity < min_sequence_identity:
            for pos in range(offset, offset + window_size):
                is_weak[segment_start + pos] = True


def build_weak_region_list(
    is_weak: Sequence[bool],
    matches: Sequence[bool],
    min_sequence_identity: float,
) -> List[WeakRegion]:
    weak_regions: List[WeakRegion] = []
    idx = 0
    n = len(is_weak)
    while idx < n:
        if is_weak[idx]:
            start = idx
            match_count = 0
            length = 0
            while idx < n and is_weak[idx]:
                if matches[idx]:
                    match_count += 1
                length += 1
                idx += 1
            identity = match_count / length if length > 0 else 0.0
            weak_regions.append(WeakRegion(start, idx - 1, identity))
        else:
            idx += 1
    return weak_regions


def nucleotide_color(base: str) -> str:
    return NUCLEOTIDE_COLORS.get(base.upper(), "#888888")


def resolve_circle_position(
    base_x: float, radius: float, existing: List[Tuple[float, float]]
) -> float:
    """Return a jittered x position that prevents circle overlap."""

    spacing = 0.05
    if all(abs(base_x - x) >= radius + r + spacing for x, r in existing):
        return base_x

    step = max(radius * 0.5, 0.2)
    for attempt in range(1, 256):
        magnitude = math.ceil(attempt / 2) * step
        direction = -1 if attempt % 2 == 0 else 1
        candidate = base_x + magnitude * direction
        if all(abs(candidate - x) >= radius + r + spacing for x, r in existing):
            return candidate

    return base_x


def construct_stream_paths(
    data: AlignmentData,
    min_gap_size: int,
    weak_regions: Sequence[WeakRegion],
    backbone_gap: float,
    bump_scale: float,
    gap_circle_scale: float,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    float,
    List[GapCircle],
    List[GapCircle],
]:
    n = len(data.query)
    global_x = np.zeros(n)
    current_global = 0.0
    for idx in range(n):
        global_x[idx] = current_global
        if not data.is_query_gap[idx] and not data.is_reference_gap[idx]:
            current_global += 1.0

    query_baseline = backbone_gap
    reference_baseline = 0.0
    query_offsets = np.zeros(n)
    reference_offsets = np.zeros(n)
    query_x = global_x.copy()
    reference_x = global_x.copy()

    gap_height_scale = 0.04
    max_beak_width = 1.2
    max_beak_height = 0.8

    for run in data.gap_runs:
        indices = range(run.start, run.end + 1)
        length = run.length
        if length <= 0:
            continue

        base_x = global_x[run.start]
        if length < min_gap_size:
            amplitude = min(max_beak_height, gap_height_scale * length)
            width = min(max_beak_width, 0.4 + 0.1 * length)
            denom = max(length - 1, 1)
            for idx, pos in enumerate(indices):
                t = idx / denom if denom > 0 else 0.5
                vertical_shape = 1.0 - abs(2.0 * t - 1.0)
                horizontal_shape = math.sin(math.pi * t)
                if run.stream == "reference":
                    query_offsets[pos] = max(
                        query_offsets[pos], amplitude * vertical_shape
                    )
                    query_x[pos] = base_x + width * horizontal_shape
                else:
                    reference_offsets[pos] = min(
                        reference_offsets[pos], -amplitude * vertical_shape
                    )
                    reference_x[pos] = base_x + width * horizontal_shape

    for region in weak_regions:
        length = region.end - region.start + 1
        if length <= 0:
            continue
        base_amplitude = min(0.5, 0.12 + 0.3 * (1.0 - region.identity))
        amplitude = max(0.0, base_amplitude * max(bump_scale, 0.0))
        for idx, pos in enumerate(range(region.start, region.end + 1)):
            t = (idx + 1) / (length + 1)
            shape = math.sin(math.pi * t)
            query_offsets[pos] = max(query_offsets[pos], amplitude * shape)
            reference_offsets[pos] = min(
                reference_offsets[pos], -amplitude * shape
            )

    query_positions = query_baseline + query_offsets
    reference_positions = reference_baseline + reference_offsets

    circle_scale = max(gap_circle_scale, 0.0)
    min_diameter = 0.08
    query_circles: List[GapCircle] = []
    reference_circles: List[GapCircle] = []
    circle_registry = {"query": [], "reference": []}

    for run in data.gap_runs:
        length = run.length
        if length <= 0 or length < min_gap_size:
            continue

        if circle_scale <= 0.0:
            continue

        anchor_idx = run.start + run.length // 2
        anchor_idx = min(max(anchor_idx, run.start), run.end)

        if run.stream == "reference":
            anchor_x = float(query_x[anchor_idx])
            anchor_y = float(query_positions[anchor_idx])
            stream_key = "query"
        else:
            anchor_x = float(reference_x[anchor_idx])
            anchor_y = float(reference_positions[anchor_idx])
            stream_key = "reference"

        diameter = max(length * circle_scale, min_diameter)
        radius = diameter / 2.0
        stem_length = max(radius * 0.25, 0.05)
        sign = 1.0 if run.stream == "reference" else -1.0

        target_list = circle_registry[stream_key]
        center_x = resolve_circle_position(anchor_x, radius, target_list)
        jitter = center_x - anchor_x
        center_y = anchor_y + sign * (radius + stem_length)

        circle = GapCircle(
            run=run,
            anchor_x=anchor_x,
            anchor_y=anchor_y,
            center_x=center_x,
            center_y=center_y,
            radius=radius,
            jitter=jitter,
        )

        target_list.append((center_x, radius))
        if stream_key == "query":
            query_circles.append(circle)
        else:
            reference_circles.append(circle)

    return (
        global_x,
        query_x,
        reference_x,
        query_positions,
        reference_positions,
        current_global,
        query_circles,
        reference_circles,
    )


def compute_tick_marks(
    local_positions: Sequence[int],
    stream_x: np.ndarray,
    stream_y: np.ndarray,
    is_gap: Sequence[bool],
    interval: int,
) -> List[Tuple[float, float, int]]:
    ticks: List[Tuple[float, float, int]] = []
    if interval <= 0:
        return ticks

    next_tick = interval
    max_local = local_positions[-1] if local_positions else 0

    for idx, local in enumerate(local_positions):
        if is_gap[idx]:
            continue
        while local >= next_tick:
            ticks.append((float(stream_x[idx]), float(stream_y[idx]), next_tick))
            next_tick += interval
        if next_tick > max_local:
            break

    return ticks


def write_stream_debug_tables(
    output: Path,
    data: AlignmentData,
    global_x: np.ndarray,
    query_x: np.ndarray,
    reference_x: np.ndarray,
    query_positions: np.ndarray,
    reference_positions: np.ndarray,
    query_circles: Sequence[GapCircle],
    reference_circles: Sequence[GapCircle],
) -> None:
    """Write TSV tables describing the computed stream coordinates."""

    def fmt_float(value: float) -> str:
        return format(value, ".8f")

    output.parent.mkdir(parents=True, exist_ok=True)
    base_name = output.stem
    query_path = output.parent / f"{base_name}_query_stream.tsv"
    reference_path = output.parent / f"{base_name}_reference_stream.tsv"

    n = len(data.query)
    global_values = [float(val) for val in global_x.tolist()]
    global_deltas: List[float] = []
    prev_value: float | None = None
    for value in global_values:
        if prev_value is None:
            global_deltas.append(0.0)
        else:
            global_deltas.append(value - prev_value)
        prev_value = value

    query_gap_lengths: Dict[int, int] = {}
    reference_gap_lengths: Dict[int, int] = {}
    for run in data.gap_runs:
        target = query_gap_lengths if run.stream == "query" else reference_gap_lengths
        for pos in range(run.start, run.end + 1):
            target[pos] = run.length

    query_circle_lookup: Dict[int, GapCircle] = {}
    for circle in query_circles:
        for pos in range(circle.run.start, circle.run.end + 1):
            query_circle_lookup[pos] = circle

    reference_circle_lookup: Dict[int, GapCircle] = {}
    for circle in reference_circles:
        for pos in range(circle.run.start, circle.run.end + 1):
            reference_circle_lookup[pos] = circle

    header = (
        "column_index\tglobal_x\tstream_x\tstream_y\tlocal_position\tcharacter"
        "\tpartner_character\tis_gap\tis_match\tis_weak\tgap_run_length"
        "\tglobal_delta\tcircle_gap_length\tcircle_center_x\tcircle_center_y"
        "\tcircle_radius\tcircle_anchor_x\tcircle_anchor_y\tcircle_jitter\n"
    )

    with query_path.open("w", encoding="utf-8") as handle:
        handle.write(header)
        for idx in range(n):
            global_value = global_values[idx]
            circle = query_circle_lookup.get(idx)
            circle_length = circle.run.length if circle else ""
            circle_center_x = fmt_float(circle.center_x) if circle else ""
            circle_center_y = fmt_float(circle.center_y) if circle else ""
            circle_radius = fmt_float(circle.radius) if circle else ""
            circle_anchor_x = fmt_float(circle.anchor_x) if circle else ""
            circle_anchor_y = fmt_float(circle.anchor_y) if circle else ""
            circle_jitter = fmt_float(circle.jitter) if circle else ""
            handle.write(
                "\t".join(
                    [
                        str(idx),
                        fmt_float(global_value),
                        fmt_float(float(query_x[idx])),
                        fmt_float(float(query_positions[idx])),
                        str(data.query_local_positions[idx]),
                        data.query[idx],
                        data.reference[idx],
                        "true" if data.is_query_gap[idx] else "false",
                        "true" if data.matches[idx] else "false",
                        "true" if data.is_weak[idx] else "false",
                        str(query_gap_lengths.get(idx, 0)),
                        fmt_float(global_deltas[idx]),
                        str(circle_length),
                        circle_center_x,
                        circle_center_y,
                        circle_radius,
                        circle_anchor_x,
                        circle_anchor_y,
                        circle_jitter,
                    ]
                )
                + "\n"
            )

    with reference_path.open("w", encoding="utf-8") as handle:
        handle.write(header)
        for idx in range(n):
            global_value = global_values[idx]
            circle = reference_circle_lookup.get(idx)
            circle_length = circle.run.length if circle else ""
            circle_center_x = fmt_float(circle.center_x) if circle else ""
            circle_center_y = fmt_float(circle.center_y) if circle else ""
            circle_radius = fmt_float(circle.radius) if circle else ""
            circle_anchor_x = fmt_float(circle.anchor_x) if circle else ""
            circle_anchor_y = fmt_float(circle.anchor_y) if circle else ""
            circle_jitter = fmt_float(circle.jitter) if circle else ""
            handle.write(
                "\t".join(
                    [
                        str(idx),
                        fmt_float(global_value),
                        fmt_float(float(reference_x[idx])),
                        fmt_float(float(reference_positions[idx])),
                        str(data.reference_local_positions[idx]),
                        data.reference[idx],
                        data.query[idx],
                        "true" if data.is_reference_gap[idx] else "false",
                        "true" if data.matches[idx] else "false",
                        "true" if data.is_weak[idx] else "false",
                        str(reference_gap_lengths.get(idx, 0)),
                        fmt_float(global_deltas[idx]),
                        str(circle_length),
                        circle_center_x,
                        circle_center_y,
                        circle_radius,
                        circle_anchor_x,
                        circle_anchor_y,
                        circle_jitter,
                    ]
                )
                + "\n"
            )


def plot_alignment(
    data: AlignmentData,
    global_x: np.ndarray,
    query_x: np.ndarray,
    reference_x: np.ndarray,
    query_positions: np.ndarray,
    reference_positions: np.ndarray,
    global_extent: float,
    width: float,
    height: float,
    dpi: int,
    output: Path,
    tick_interval: int,
    query_circles: Sequence[GapCircle],
    reference_circles: Sequence[GapCircle],
    backbone_thickness: float,
    mismatch_line_width: float,
) -> None:
    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)

    line_width = max(backbone_thickness, 0.0)
    ax.plot(
        query_x,
        query_positions,
        color="#222222",
        linewidth=line_width if line_width > 0 else 0.0,
        label="Query",
        zorder=5,
    )
    ax.plot(
        reference_x,
        reference_positions,
        color="#222222",
        linewidth=line_width if line_width > 0 else 0.0,
        label="Reference",
        zorder=5,
    )

    tick_length = 0.06
    if tick_interval > 0:
        query_ticks = compute_tick_marks(
            data.query_local_positions,
            query_x,
            query_positions,
            data.is_query_gap,
            tick_interval,
        )
        reference_ticks = compute_tick_marks(
            data.reference_local_positions,
            reference_x,
            reference_positions,
            data.is_reference_gap,
            tick_interval,
        )

        for x_pos, y_pos, value in query_ticks:
            ax.plot([x_pos, x_pos], [y_pos, y_pos + tick_length], color="#555555", linewidth=1.0)
            ax.text(
                x_pos,
                y_pos + tick_length * 1.4,
                f"{value:,}",
                fontsize=6,
                ha="center",
                va="bottom",
                color="#333333",
            )

        for x_pos, y_pos, value in reference_ticks:
            ax.plot([x_pos, x_pos], [y_pos, y_pos - tick_length], color="#555555", linewidth=1.0)
            ax.text(
                x_pos,
                y_pos - tick_length * 1.4,
                f"{value:,}",
                fontsize=6,
                ha="center",
                va="top",
                color="#333333",
            )

    weak_spans = [(region.start, region.end) for region in data.weak_regions]
    for start, end in weak_spans:
        span_start = global_x[start]
        span_end = global_x[end]
        if math.isclose(span_start, span_end):
            span_end = span_start + 1e-3
        ax.axvspan(span_start, span_end, color="#f0f0f0", alpha=0.7, zorder=0)

    for idx, (x_pos, q_char, r_char) in enumerate(
        zip(global_x, data.query, data.reference)
    ):
        if data.is_query_gap[idx] or data.is_reference_gap[idx]:
            continue
        if any(region.start <= idx <= region.end for region in data.weak_regions):
            continue
        if q_char == r_char:
            continue
        top = query_positions[idx]
        bottom = reference_positions[idx]
        middle = (top + bottom) / 2.0
        rung_width = max(mismatch_line_width, 0.0)
        if rung_width > 0:
            ax.plot(
                [x_pos, x_pos],
                [top, middle],
                color=nucleotide_color(q_char),
                linewidth=rung_width,
                zorder=3,
            )
            ax.plot(
                [x_pos, x_pos],
                [middle, bottom],
                color=nucleotide_color(r_char),
                linewidth=rung_width,
                zorder=3,
            )

    connector_color = "#555555"
    circle_edge_color = "#222222"
    circle_face_color = "#ffffff"

    x_min = global_x[0] if len(global_x) > 0 else 0.0
    x_max = max(global_extent, x_min + 1.0)

    top_candidates: List[float] = [float(np.max(query_positions))] if len(query_positions) else []
    top_candidates.extend(circle.center_y + circle.radius for circle in query_circles)
    bottom_candidates: List[float] = [float(np.min(reference_positions))] if len(reference_positions) else []
    bottom_candidates.extend(circle.center_y - circle.radius for circle in reference_circles)

    y_min = min(bottom_candidates) - 0.2 if bottom_candidates else -0.5
    y_max = max(top_candidates) + 0.2 if top_candidates else 1.5

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Ensure the renderer is ready so coordinate transforms are accurate.
    fig.canvas.draw()

    x0y0 = ax.transData.transform((0.0, 0.0))
    x1y0 = ax.transData.transform((1.0, 0.0))
    x0y1 = ax.transData.transform((0.0, 1.0))
    px_per_x = x1y0[0] - x0y0[0]
    px_per_y = x0y1[1] - x0y0[1]
    y_scale = px_per_x / px_per_y if px_per_y != 0 else 1.0

    for circle in query_circles:
        contact_y = circle.center_y - circle.radius
        control_point = (
            (circle.anchor_x + circle.center_x) / 2.0,
            circle.anchor_y + (circle.center_y - circle.anchor_y) * 0.45,
        )
        path = MplPath(
            [
                (circle.anchor_x, circle.anchor_y),
                control_point,
                (circle.center_x, contact_y),
            ],
            [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3],
        )
        ax.add_patch(
            patches.PathPatch(
                path,
                fill=False,
                linewidth=1.0,
                color=connector_color,
                zorder=4,
            )
        )
        circle_patch = patches.Ellipse(
            (circle.center_x, circle.center_y),
            width=2 * circle.radius,
            height=2 * circle.radius * y_scale,
            facecolor=circle_face_color,
            edgecolor=circle_edge_color,
            linewidth=1.2,
            zorder=4.5,
        )
        ax.add_patch(circle_patch)
        ax.text(
            circle.center_x,
            circle.center_y + circle.radius + 0.05,
            f"+{circle.run.length:,}",
            fontsize=6,
            ha="center",
            va="bottom",
            color="#222222",
            zorder=6,
        )

    for circle in reference_circles:
        contact_y = circle.center_y + circle.radius
        control_point = (
            (circle.anchor_x + circle.center_x) / 2.0,
            circle.anchor_y + (circle.center_y - circle.anchor_y) * 0.45,
        )
        path = MplPath(
            [
                (circle.anchor_x, circle.anchor_y),
                control_point,
                (circle.center_x, contact_y),
            ],
            [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3],
        )
        ax.add_patch(
            patches.PathPatch(
                path,
                fill=False,
                linewidth=1.0,
                color=connector_color,
                zorder=4,
            )
        )
        circle_patch = patches.Ellipse(
            (circle.center_x, circle.center_y),
            width=2 * circle.radius,
            height=2 * circle.radius * y_scale,
            facecolor=circle_face_color,
            edgecolor=circle_edge_color,
            linewidth=1.2,
            zorder=4.5,
        )
        ax.add_patch(circle_patch)
        ax.text(
            circle.center_x,
            circle.center_y - circle.radius - 0.05,
            f"+{circle.run.length:,}",
            fontsize=6,
            ha="center",
            va="top",
            color="#222222",
            zorder=6,
        )
    ax.axis("off")

    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    try:
        query, reference = parse_fasta_pair(args.input)
    except Exception as exc:  # pragma: no cover - user input validation
        print(f"Error while parsing FASTA file: {exc}", file=sys.stderr)
        return 1

    try:
        data = compute_alignment_data(
            query,
            reference,
            args.min_gap_size,
            args.block_size,
            args.min_sequence_identity,
            args.window_size,
        )
        (
            global_x,
            query_x,
            reference_x,
            query_positions,
            reference_positions,
            global_extent,
            query_circles,
            reference_circles,
        ) = construct_stream_paths(
            data,
            args.min_gap_size,
            data.weak_regions,
            args.backbone_gap,
            args.bump_scale,
            args.gap_circle_scale,
        )
        write_stream_debug_tables(
            args.output,
            data,
            global_x,
            query_x,
            reference_x,
            query_positions,
            reference_positions,
            query_circles,
            reference_circles,
        )
        plot_alignment(
            data,
            global_x,
            query_x,
            reference_x,
            query_positions,
            reference_positions,
            global_extent,
            args.width,
            args.height,
            args.dpi,
            args.output,
            args.tick_interval,
            query_circles,
            reference_circles,
            args.backbone_thickness,
            args.mismatch_line_width,
        )
    except Exception as exc:  # pragma: no cover - runtime safety
        print(f"Error while creating visualization: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main(sys.argv[1:]))
