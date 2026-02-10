from __future__ import annotations


def compute_alignment_data(
    query: str,
    reference: str,
    query_name: str,
    reference_name: str,
    min_gap_size: int,
    block_size: int,
    min_sequence_identity: float,
    window_size: int,
):
    from align_viz import compute_alignment_data as _compute_alignment_data

    return _compute_alignment_data(
        query,
        reference,
        query_name,
        reference_name,
        min_gap_size,
        block_size,
        min_sequence_identity,
        window_size,
    )


def construct_stream_paths(
    data,
    weak_regions,
    backbone_gap: float,
    bump_scale: float,
    gap_max_height: float,
    gap_width: float,
    min_gap_size: int,
    gap_height_scale: float,
    indel_height_scale: float,
):
    from align_viz import construct_stream_paths as _construct_stream_paths

    return _construct_stream_paths(
        data,
        weak_regions,
        backbone_gap,
        bump_scale,
        gap_max_height,
        gap_width,
        min_gap_size,
        gap_height_scale,
        indel_height_scale,
    )
