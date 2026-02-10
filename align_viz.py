#!/usr/bin/env python3
"""Command line visualization tool for pairwise DNA alignment outputs.

The tool expects a FASTA file that contains two aligned sequences. The
visualization renders two horizontal streams (query on top, reference on
bottom) and decorates them with glyphs describing indels, weakly aligned
regions, and mismatch ladder rungs.
"""
from __future__ import annotations

import argparse
import hashlib
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from core.params import RenderParams
from core.service import export_to_file_for_cli, prepare_session, write_debug_tables_for_cli


NUCLEOTIDE_COLORS: Dict[str, str] = {
    "A": "#4daf4a",  # green
    "C": "#377eb8",  # blue
    "G": "#000000",  # black
    "T": "#e41a1c",  # red
    "U": "#e41a1c",  # treat U like T
}


DEFAULT_ANNOTATION_THICKNESS = 4.0
DEFAULT_ANNOTATION_ALPHA = 0.85
DEFAULT_REF_ANNOTATION_COLOR = "#984ea3"
DEFAULT_QUERY_ANNOTATION_COLOR = "#ff7f00"
DEFAULT_ANNOTATION_MAX_LAYERS = 3
DEFAULT_ANNOTATION_SPACING = 0.8
ANNOTATION_LABEL_OFFSET = 0.18
DEFAULT_ANNOTATION_LABEL_JITTER = 0.35
AUTO_WIDTH_TOKEN = "auto"
AUTO_WIDTH_REFERENCE_GLOBAL_X = 80_278.0
AUTO_WIDTH_REFERENCE_INCHES = 125.0
AUTO_WIDTH_MIN_INCHES = 8.0


def parse_optional_positive_float(value: str) -> Optional[float]:
    normalized = value.strip().lower()
    if normalized in {"na", "null"}:
        return None
    try:
        parsed = float(value)
    except ValueError as exc:  # pragma: no cover - argparse error propagation
        raise argparse.ArgumentTypeError(
            "value must be a floating-point number or 'NA'/'NULL'"
        ) from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def parse_nonnegative_float(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:  # pragma: no cover - argparse error propagation
        raise argparse.ArgumentTypeError("value must be a floating-point number") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return parsed


def parse_width_arg(value: str) -> Union[float, str]:
    normalized = value.strip().lower()
    if normalized == AUTO_WIDTH_TOKEN:
        return AUTO_WIDTH_TOKEN
    try:
        parsed = float(value)
    except ValueError as exc:  # pragma: no cover - argparse error propagation
        raise argparse.ArgumentTypeError(
            "width must be a floating-point number or 'auto'"
        ) from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("width must be positive")
    return parsed


def resolve_figure_width(width_arg: Union[float, str], global_extent: float) -> float:
    if isinstance(width_arg, float):
        return width_arg
    if width_arg != AUTO_WIDTH_TOKEN:
        raise ValueError(f"Unsupported width argument: {width_arg!r}")
    if global_extent <= 0:
        return AUTO_WIDTH_MIN_INCHES
    scaled = AUTO_WIDTH_REFERENCE_INCHES * (global_extent / AUTO_WIDTH_REFERENCE_GLOBAL_X)
    return max(AUTO_WIDTH_MIN_INCHES, scaled)


def parse_positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:  # pragma: no cover - argparse error propagation
        raise argparse.ArgumentTypeError("value must be an integer") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def parse_alpha(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:  # pragma: no cover - argparse error propagation
        raise argparse.ArgumentTypeError("value must be a floating-point number") from exc
    if not (0.0 <= parsed <= 1.0):
        raise argparse.ArgumentTypeError("alpha must fall between 0 and 1")
    return parsed


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
class GapLabel:
    x: float
    y: float
    direction: int  # +1 for upward labels, -1 for downward
    length: int
    amplitude: float


@dataclass
class AlignmentData:
    query: str
    reference: str
    query_name: str
    reference_name: str
    matches: List[bool]
    is_query_gap: List[bool]
    is_reference_gap: List[bool]
    query_local_positions: List[int]
    reference_local_positions: List[int]
    is_weak: List[bool]
    gap_runs: List[GapRun]
    weak_regions: List[WeakRegion]


@dataclass
class GeneFeature:
    name: str
    start: int
    end: int
    kind: str


@dataclass
class GeneAnnotation:
    gene_id: str
    start: int
    end: int
    direction: str
    features: List[GeneFeature] = field(default_factory=list)


@dataclass
class GeneLabelInfo:
    x: float
    y: float
    annotation: GeneAnnotation
    layer: int


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize pairwise DNA sequence alignments with loop/bump glyphs."
    )
    parser.add_argument("input", type=Path, help="FASTA file containing the alignment")
    parser.add_argument(
        "width",
        type=parse_width_arg,
        help=(
            "Figure width (inches) or 'auto' to scale from alignment span "
            f"(calibrated so global_x≈{int(AUTO_WIDTH_REFERENCE_GLOBAL_X):,} maps to "
            f"{AUTO_WIDTH_REFERENCE_INCHES:g} inches)"
        ),
    )
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
        "--mismatch-line-width",
        type=float,
        default=1.2,
        help="Line width used for mismatch ladder rungs",
    )
    parser.add_argument(
        "--gap-max-height",
        type=float,
        default=0.8,
        help="Maximum amplitude (in data units) for gap beak glyphs",
    )
    parser.add_argument(
        "--gap-width",
        type=float,
        default=0.0,
        help="Horizontal width (in nucleotides) assigned to each gap column",
    )
    parser.add_argument(
        "--gap-height-scale",
        type=float,
        default=0.04,
        help="Multiplier applied to gap amplitudes before clamping to the maximum height",
    )
    parser.add_argument(
        "--indel-height-scale",
        type=float,
        default=0.04,
        help="Multiplier applied to short indel amplitudes (length < min gap size)",
    )
    parser.add_argument(
        "--gap-label-size",
        type=parse_optional_positive_float,
        default=8.0,
        help=(
            "Font size used for '+x bp' gap labels; set to NA or NULL to disable label rendering"
        ),
    )
    parser.add_argument(
        "--query-annotation",
        type=Path,
        help="Path to a query stream gene annotation file",
    )
    parser.add_argument(
        "--reference-annotation",
        type=Path,
        help="Path to a reference stream gene annotation file",
    )
    parser.add_argument(
        "--annotation-label-size",
        type=parse_optional_positive_float,
        default=10.0,
        help=(
            "Font size for gene annotation labels; set to NA or NULL to hide the labels"
        ),
    )
    parser.add_argument(
        "--annotation-thickness",
        type=parse_nonnegative_float,
        default=DEFAULT_ANNOTATION_THICKNESS,
        help="Line width used for coding segments in gene annotations",
    )
    parser.add_argument(
        "--annotation-alpha",
        type=parse_alpha,
        default=DEFAULT_ANNOTATION_ALPHA,
        help="Alpha transparency applied to annotation overlays (0-1)",
    )
    parser.add_argument(
        "--reference-annotation-color",
        default=DEFAULT_REF_ANNOTATION_COLOR,
        help="Stroke color used for reference annotations (any Matplotlib color)",
    )
    parser.add_argument(
        "--query-annotation-color",
        default=DEFAULT_QUERY_ANNOTATION_COLOR,
        help="Stroke color used for query annotations (any Matplotlib color)",
    )
    parser.add_argument(
        "--annotation-label-jitter",
        type=parse_nonnegative_float,
        default=DEFAULT_ANNOTATION_LABEL_JITTER,
        help="Horizontal jitter amplitude for annotation labels (set to 0 to disable)",
    )
    parser.add_argument(
        "--annotation-max-layers",
        type=parse_positive_int,
        default=DEFAULT_ANNOTATION_MAX_LAYERS,
        help="Maximum number of parallel annotation layers rendered per stream",
    )
    parser.add_argument(
        "--annotation-spacing",
        type=parse_nonnegative_float,
        default=DEFAULT_ANNOTATION_SPACING,
        help="Vertical spacing (in data units) between stacked annotation layers",
    )
    return parser.parse_args(argv)


def parse_fasta_pair(path: Path) -> Tuple[str, str, str, str]:
    records: List[Tuple[str, List[str]]] = []
    current_header: str | None = None
    current_chunks: List[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, current_chunks))
                current_header = line[1:].strip()
                current_chunks = []
            else:
                if current_header is None:
                    raise ValueError("FASTA file must start with a header line beginning with '>'")
                current_chunks.append(line)

    if current_header is not None:
        records.append((current_header, current_chunks))

    if len(records) != 2:
        raise ValueError(
            f"Expected exactly two sequences in the FASTA file, found {len(records)}"
        )

    (header1, chunks1), (header2, chunks2) = records
    seq1 = "".join(chunks1).upper()
    seq2 = "".join(chunks2).upper()

    if len(seq1) != len(seq2):
        raise ValueError(
            "Aligned sequences must have the same length (including gap characters)."
        )

    def is_reference(name: str) -> bool:
        lowered = name.lower()
        return "reference" in lowered or "ref" in lowered

    def is_query(name: str) -> bool:
        lowered = name.lower()
        return "query" in lowered

    if is_reference(header1) and not is_reference(header2):
        reference_seq, reference_name = seq1, header1
        query_seq, query_name = seq2, header2
    elif is_reference(header2) and not is_reference(header1):
        reference_seq, reference_name = seq2, header2
        query_seq, query_name = seq1, header1
    elif is_query(header1) and not is_query(header2):
        query_seq, query_name = seq1, header1
        reference_seq, reference_name = seq2, header2
    elif is_query(header2) and not is_query(header1):
        query_seq, query_name = seq2, header2
        reference_seq, reference_name = seq1, header1
    else:
        query_seq, query_name = seq1, header1
        reference_seq, reference_name = seq2, header2

    return query_seq, reference_seq, query_name, reference_name


def normalize_feature_kind(name: str) -> str:
    lowered = name.strip().lower()
    if lowered == "exon":
        return "exon"
    if lowered == "utr":
        return "utr"
    return lowered or name


def finalize_gene_annotation(
    gene_id: str,
    start: int,
    end: int,
    direction: str,
    features: List[GeneFeature],
) -> GeneAnnotation:
    if start > end:
        raise ValueError(f"Gene '{gene_id}' has start greater than end")
    direction_symbol = direction.strip()
    if direction_symbol not in {"<", ">"}:
        raise ValueError(
            f"Gene '{gene_id}' must specify direction '<' or '>', found '{direction}'"
        )

    sorted_features = sorted(features, key=lambda feature: (feature.start, feature.end))
    merged: List[GeneFeature] = []
    current = start - 1
    for feature in sorted_features:
        if feature.start > feature.end:
            raise ValueError(
                f"Feature '{feature.name}' in gene '{gene_id}' has start greater than end"
            )
        if feature.start < start or feature.end > end:
            raise ValueError(
                f"Feature '{feature.name}' in gene '{gene_id}' falls outside the gene bounds"
            )
        if feature.start <= current:
            raise ValueError(
                f"Features in gene '{gene_id}' overlap or are unsorted around position {feature.start}"
            )
        if feature.start > current + 1:
            merged.append(
                GeneFeature(
                    name="intron",
                    start=current + 1,
                    end=feature.start - 1,
                    kind="intron",
                )
            )
        merged.append(feature)
        current = feature.end

    if current < end:
        merged.append(
            GeneFeature(name="intron", start=current + 1, end=end, kind="intron")
        )

    return GeneAnnotation(
        gene_id=gene_id,
        start=start,
        end=end,
        direction=direction_symbol,
        features=merged,
    )


def parse_annotation_file(path: Path) -> List[GeneAnnotation]:
    annotations: List[GeneAnnotation] = []
    if path is None:
        return annotations

    if not path.exists():
        raise FileNotFoundError(f"Annotation file '{path}' was not found")

    current_gene_id: Optional[str] = None
    current_start: Optional[int] = None
    current_end: Optional[int] = None
    current_direction: Optional[str] = None
    current_features: List[GeneFeature] = []

    def flush_current() -> None:
        nonlocal current_gene_id, current_start, current_end, current_direction, current_features
        if current_gene_id is None:
            return
        if current_start is None or current_end is None or current_direction is None:
            raise ValueError(
                f"Gene '{current_gene_id}' is missing coordinate or direction information"
            )
        annotations.append(
            finalize_gene_annotation(
                current_gene_id,
                current_start,
                current_end,
                current_direction,
                current_features,
            )
        )
        current_gene_id = None
        current_start = None
        current_end = None
        current_direction = None
        current_features = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                flush_current()
                continue
            parts = line.split()
            if current_gene_id is None:
                if len(parts) != 4:
                    raise ValueError(
                        f"Gene header lines must have four columns, found: '{line}'"
                    )
                gene_id, start_text, end_text, direction = parts
                try:
                    start_val = int(start_text)
                    end_val = int(end_text)
                except ValueError as exc:
                    raise ValueError(
                        f"Gene '{gene_id}' has non-integer bounds ('{start_text}', '{end_text}')"
                    ) from exc
                current_gene_id = gene_id
                current_start = start_val
                current_end = end_val
                current_direction = direction
                current_features = []
            else:
                if len(parts) != 3:
                    raise ValueError(
                        f"Feature rows must have three columns, found: '{line}'"
                    )
                feature_name, start_text, end_text = parts
                try:
                    feature_start = int(start_text)
                    feature_end = int(end_text)
                except ValueError as exc:
                    raise ValueError(
                        f"Feature '{feature_name}' in gene '{current_gene_id}' has non-integer bounds"
                    ) from exc
                current_features.append(
                    GeneFeature(
                        name=feature_name,
                        start=feature_start,
                        end=feature_end,
                        kind=normalize_feature_kind(feature_name),
                    )
                )

    flush_current()
    return annotations


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
    query_name: str,
    reference_name: str,
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
        query_name=query_name,
        reference_name=reference_name,
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


def construct_stream_paths(
    data: AlignmentData,
    weak_regions: Sequence[WeakRegion],
    backbone_gap: float,
    bump_scale: float,
    gap_max_height: float,
    gap_column_width: float,
    min_gap_size: int,
    gap_height_scale: float,
    indel_height_scale: float,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    float,
    List[GapLabel],
]:
    n = len(data.query)
    global_x = np.zeros(n)
    query_x = np.zeros(n)
    reference_x = np.zeros(n)
    query_offsets = np.zeros(n)
    reference_offsets = np.zeros(n)

    query_baseline = backbone_gap
    reference_baseline = 0.0

    gap_width = max(gap_column_width, 0.0)
    height_scale = max(gap_height_scale, 0.0)
    indel_scale = max(indel_height_scale, 0.0)
    max_beak_height = max(gap_max_height, 0.0)
    gap_labels: List[GapLabel] = []

    gap_runs_by_start = {run.start: run for run in data.gap_runs}

    def compute_label_jitter(run: GapRun, span: float) -> float:
        seed = (run.start << 16) ^ (run.end << 1) ^ (1 if run.stream == "reference" else 0)
        raw = ((seed * 1103515245 + 12345) & 0x7FFFFFFF) / 0x7FFFFFFF
        jitter_extent = max(span * 0.35, 0.5)
        return (raw * 2.0 - 1.0) * jitter_extent

    idx = 0
    current_global = 0.0
    while idx < n:
        run = gap_runs_by_start.get(idx)
        if run and run.length > 0:
            is_large_gap = run.length >= max(min_gap_size, 1)
            span = gap_width if is_large_gap else 0.0
            start_x = current_global
            length = run.length
            denom = max(length - 1, 1)
            scale = height_scale if is_large_gap else indel_scale
            amplitude = min(max_beak_height, scale * length) if length > 0 else 0.0
            for offset, pos in enumerate(range(run.start, run.end + 1)):
                if length == 1:
                    edge_t = 0.5
                    midpoint_t = 0.5
                else:
                    edge_t = offset / denom
                    midpoint_t = (offset + 0.5) / length
                weight = 4.0 * midpoint_t * (1.0 - midpoint_t)
                if length > 1 and (offset == 0 or offset == length - 1):
                    weight = 0.0
                x_value = start_x + span * edge_t if span > 0 else start_x
                global_x[pos] = x_value
                query_x[pos] = x_value
                reference_x[pos] = x_value
                if run.stream == "reference":
                    query_offsets[pos] = max(query_offsets[pos], amplitude * weight)
                else:
                    reference_offsets[pos] = min(
                        reference_offsets[pos], -amplitude * weight
                    )
            if is_large_gap and amplitude > 0:
                base_x = start_x + (span * 0.5 if span > 0 else 0.0)
                apex_x = base_x + compute_label_jitter(run, span)
                if run.stream == "reference":
                    apex_y = query_baseline + amplitude
                    direction = 1
                else:
                    apex_y = reference_baseline - amplitude
                    direction = -1
                gap_labels.append(
                    GapLabel(
                        x=apex_x,
                        y=apex_y,
                        direction=direction,
                        length=run.length,
                        amplitude=amplitude,
                    )
                )
            current_global = start_x + span
            idx = run.end + 1
            continue

        x_value = current_global
        global_x[idx] = x_value
        query_x[idx] = x_value
        reference_x[idx] = x_value
        current_global += 1.0
        idx += 1

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

    if n:
        last_x = float(global_x[-1])
        global_extent = max(current_global, last_x + 1.0)
    else:
        global_extent = 0.0

    return (
        global_x,
        query_x,
        reference_x,
        query_positions,
        reference_positions,
        global_extent,
        gap_labels,
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

    header = (
        "column_index\tglobal_x\tstream_x\tstream_y\tlocal_position\tcharacter"
        "\tpartner_character\tis_gap\tis_match\tis_weak\tgap_run_length"
        "\tglobal_delta\n"
    )

    with query_path.open("w", encoding="utf-8") as handle:
        handle.write(header)
        for idx in range(n):
            global_value = global_values[idx]
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
                    ]
                )
                + "\n"
            )

    with reference_path.open("w", encoding="utf-8") as handle:
        handle.write(header)
        for idx in range(n):
            global_value = global_values[idx]
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
                    ]
                )
                + "\n"
            )


def build_local_index_map(
    local_positions: Sequence[int], is_gap: Sequence[bool]
) -> Dict[int, List[int]]:
    mapping: Dict[int, List[int]] = {}
    for idx, (local, gap) in enumerate(zip(local_positions, is_gap)):
        if gap:
            continue
        mapping.setdefault(local, []).append(idx)
    return mapping


def collect_indices_for_range(
    mapping: Dict[int, List[int]], start: int, end: int
) -> List[int]:
    if start > end:
        return []
    collected: List[int] = []
    for value in range(start, end + 1):
        collected.extend(mapping.get(value, []))
    return sorted(set(collected))


def split_indices_into_runs(indices: Sequence[int]) -> List[List[int]]:
    if not indices:
        return []
    runs: List[List[int]] = []
    current_run: List[int] = [indices[0]]
    for idx in indices[1:]:
        if idx == current_run[-1] + 1:
            current_run.append(idx)
        else:
            runs.append(current_run)
            current_run = [idx]
    runs.append(current_run)
    return runs


def indices_to_segment(
    run: Sequence[int], stream_x: np.ndarray, stream_y: np.ndarray
) -> Optional[np.ndarray]:
    if not run:
        return None
    if len(run) == 1:
        idx = run[0]
        x_pos = float(stream_x[idx])
        y_pos = float(stream_y[idx])
        epsilon = 0.35
        return np.array([[x_pos - epsilon, y_pos], [x_pos + epsilon, y_pos]])
    x_values = stream_x[list(run)]
    y_values = stream_y[list(run)]
    return np.column_stack((x_values, y_values))


def format_gene_label(annotation: GeneAnnotation) -> str:
    if annotation.direction == "<":
        return f"< {annotation.gene_id}"
    return f"{annotation.gene_id} >"


def compute_annotation_label_jitter(
    annotation: GeneAnnotation, amplitude: float
) -> float:
    if amplitude <= 0:
        return 0.0
    key = f"{annotation.gene_id}|{annotation.start}|{annotation.end}|{annotation.direction}"
    digest = hashlib.blake2b(key.encode("utf-8"), digest_size=8).digest()
    value = int.from_bytes(digest, "big") / float(1 << 64)
    return (value * 2.0 - 1.0) * amplitude


def assign_annotation_layers(
    annotations: Sequence[GeneAnnotation], max_layers: int
) -> Dict[int, int]:
    if max_layers <= 0:
        return {}
    layer_end: List[float] = [-math.inf] * max_layers
    assignments: Dict[int, int] = {}
    sorted_indices = sorted(
        range(len(annotations)), key=lambda idx: (annotations[idx].start, annotations[idx].end, idx)
    )
    for original_idx in sorted_indices:
        annotation = annotations[original_idx]
        assigned_layer: Optional[int] = None
        for layer in range(max_layers):
            if annotation.start > layer_end[layer]:
                assigned_layer = layer
                layer_end[layer] = annotation.end
                break
        if assigned_layer is not None:
            assignments[original_idx] = assigned_layer
    return assignments


def prepare_annotation_drawables(
    annotations: Sequence[GeneAnnotation],
    local_positions: Sequence[int],
    is_gap: Sequence[bool],
    stream_x: np.ndarray,
    stream_y: np.ndarray,
    max_layers: int,
    spacing: float,
    vertical_direction: int,
) -> Tuple[List[np.ndarray], List[np.ndarray], List[GeneLabelInfo], int]:
    coding_segments: List[np.ndarray] = []
    noncoding_segments: List[np.ndarray] = []
    labels: List[GeneLabelInfo] = []
    max_used_layer = -1

    index_map = build_local_index_map(local_positions, is_gap)
    assignments = assign_annotation_layers(annotations, max_layers)

    for idx, annotation in enumerate(annotations):
        layer = assignments.get(idx)
        if layer is None:
            continue
        max_used_layer = max(max_used_layer, layer)
        offset = vertical_direction * spacing * layer if spacing > 0 else 0.0
        gene_indices = collect_indices_for_range(index_map, annotation.start, annotation.end)
        if gene_indices:
            mid_index = gene_indices[len(gene_indices) // 2]
            labels.append(
                GeneLabelInfo(
                    x=float(stream_x[mid_index]),
                    y=float(stream_y[mid_index] + offset),
                    annotation=annotation,
                    layer=layer,
                )
            )
        for feature in annotation.features:
            feature_indices = collect_indices_for_range(index_map, feature.start, feature.end)
            for run in split_indices_into_runs(feature_indices):
                segment = indices_to_segment(run, stream_x, stream_y)
                if segment is None:
                    continue
                if offset != 0.0:
                    segment = segment.copy()
                    segment[:, 1] = segment[:, 1] + offset
                if feature.kind == "exon":
                    coding_segments.append(segment)
                else:
                    noncoding_segments.append(segment)

    return coding_segments, noncoding_segments, labels, max_used_layer


def plot_alignment(
    data: AlignmentData,
    global_x: np.ndarray,
    query_x: np.ndarray,
    reference_x: np.ndarray,
    query_positions: np.ndarray,
    reference_positions: np.ndarray,
    global_extent: float,
    gap_labels: Sequence[GapLabel],
    width: float,
    height: float,
    dpi: int,
    output: Path,
    tick_interval: int,
    backbone_thickness: float,
    mismatch_line_width: float,
    gap_label_size: Optional[float],
    query_annotations: Sequence[GeneAnnotation],
    reference_annotations: Sequence[GeneAnnotation],
    annotation_label_size: Optional[float],
    annotation_thickness: float,
    annotation_alpha: float,
    reference_annotation_color: str,
    query_annotation_color: str,
    annotation_label_jitter: float,
    annotation_max_layers: int,
    annotation_spacing: float,
    x_window: Optional[Tuple[float, float]] = None,
) -> None:
    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)

    line_width = max(backbone_thickness, 0.0)
    ax.plot(
        query_x,
        query_positions,
        color="#222222",
        linewidth=line_width if line_width > 0 else 0.0,
        label=data.query_name or "Query",
        zorder=5,
    )
    ax.plot(
        reference_x,
        reference_positions,
        color="#222222",
        linewidth=line_width if line_width > 0 else 0.0,
        label=data.reference_name or "Reference",
        zorder=5,
    )

    (
        query_coding_segments,
        query_noncoding_segments,
        query_label_info,
        query_max_layer,
    ) = prepare_annotation_drawables(
        query_annotations,
        data.query_local_positions,
        data.is_query_gap,
        query_x,
        query_positions,
        annotation_max_layers,
        annotation_spacing,
        1,
    )
    (
        reference_coding_segments,
        reference_noncoding_segments,
        reference_label_info,
        reference_max_layer,
    ) = prepare_annotation_drawables(
        reference_annotations,
        data.reference_local_positions,
        data.is_reference_gap,
        reference_x,
        reference_positions,
        annotation_max_layers,
        annotation_spacing,
        -1,
    )

    def add_annotation_segments(
        segments: Sequence[np.ndarray], color: str, width: float
    ) -> None:
        if not segments:
            return
        effective_width = max(width, 0.0)
        if effective_width <= 0:
            return
        collection = LineCollection(
            segments,
            colors=[color],
            linewidths=effective_width,
            alpha=annotation_alpha,
            zorder=6,
        )
        ax.add_collection(collection)

    add_annotation_segments(
        query_coding_segments, query_annotation_color, annotation_thickness
    )
    add_annotation_segments(
        query_noncoding_segments, query_annotation_color, backbone_thickness
    )
    add_annotation_segments(
        reference_coding_segments, reference_annotation_color, annotation_thickness
    )
    add_annotation_segments(
        reference_noncoding_segments, reference_annotation_color, backbone_thickness
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

    if gap_label_size is not None:
        label_offset_base = 0.12
        for label in gap_labels:
            text = f"+{label.length:,} bp"
            offset = label_offset_base + 0.05 * label.amplitude
            text_y = label.y + label.direction * offset
            vertical_alignment = "bottom" if label.direction > 0 else "top"
            ax.text(
                label.x,
                text_y,
                text,
                fontsize=gap_label_size,
                ha="center",
                va=vertical_alignment,
                color="#222222",
                zorder=6,
            )

    if annotation_label_size is not None:
        def draw_labels(
            entries: Sequence[GeneLabelInfo],
            vertical_direction: int,
            color: str,
        ) -> None:
            for info in entries:
                jitter = compute_annotation_label_jitter(
                    info.annotation, annotation_label_jitter
                )
                ax.text(
                    info.x + jitter,
                    info.y + vertical_direction * ANNOTATION_LABEL_OFFSET,
                    format_gene_label(info.annotation),
                    fontsize=annotation_label_size,
                    ha="center",
                    va="bottom" if vertical_direction > 0 else "top",
                    color=color,
                    zorder=7,
                )

        draw_labels(query_label_info, 1, query_annotation_color)
        draw_labels(reference_label_info, -1, reference_annotation_color)

    x_min = global_x[0] if len(global_x) > 0 else 0.0
    x_max = max(global_extent, x_min + 1.0)

    y_min_candidates: List[float] = []
    if len(reference_positions):
        y_min_candidates.append(float(np.min(reference_positions)))
    if len(query_positions):
        y_min_candidates.append(float(np.min(query_positions)))
    y_max_candidates: List[float] = []
    if len(query_positions):
        y_max_candidates.append(float(np.max(query_positions)))
    if len(reference_positions):
        y_max_candidates.append(float(np.max(reference_positions)))

    y_min = min(y_min_candidates) - 0.2 if y_min_candidates else -0.5
    y_max = max(y_max_candidates) + 0.2 if y_max_candidates else 1.5

    if annotation_spacing > 0:
        if query_max_layer >= 0:
            y_max += annotation_spacing * query_max_layer
        if reference_max_layer >= 0:
            y_min -= annotation_spacing * reference_max_layer

    if annotation_label_size is not None:
        if query_max_layer >= 0:
            y_max += ANNOTATION_LABEL_OFFSET
        if reference_max_layer >= 0:
            y_min -= ANNOTATION_LABEL_OFFSET

    if x_window is not None:
        window_start, window_end = x_window
        x_min = max(x_min, min(window_start, window_end))
        x_max = min(x_max, max(window_start, window_end))
        if x_max <= x_min:
            x_max = x_min + 1.0

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    ax.axis("off")

    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    try:
        params = RenderParams.from_cli_args(args)
        session = prepare_session(
            input_path=args.input,
            params=params,
            query_annotation_path=args.query_annotation,
            reference_annotation_path=args.reference_annotation,
        )
    except Exception as exc:  # pragma: no cover - user input validation
        print(f"Error while parsing input files: {exc}", file=sys.stderr)
        return 1

    try:
        write_debug_tables_for_cli(session, args.output)
        export_to_file_for_cli(session, args.output)
    except Exception as exc:  # pragma: no cover - runtime safety
        print(f"Error while creating visualization: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main(sys.argv[1:]))
