from __future__ import annotations

import re
from pathlib import Path
from typing import List, Tuple


_BLAST_SCORE_RE = re.compile(r"Score:\s*([0-9]+(?:\.[0-9]+)?)\s*bits", re.IGNORECASE)


def _first_nonempty_line(path: Path) -> str:
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line:
                return line
    return ""


def parse_fasta_pair(path: Path) -> Tuple[str, str, str, str]:
    from core.engine import parse_fasta_pair as _parse_fasta_pair

    return _parse_fasta_pair(path)


def _parse_alignment_row(line: str, expected_label: str) -> str:
    parts = line.split()
    if len(parts) < 4 or parts[0] != expected_label:
        raise ValueError(f"Malformed BLAST alignment row: '{line}'")
    return parts[2]


def parse_blast_pairwise(path: Path) -> Tuple[str, str, str, str]:
    with path.open("r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle]

    query_name = ""
    reference_name = ""
    for raw_line in lines:
        line = raw_line.strip()
        if not query_name and line.startswith("Query:"):
            query_name = line[len("Query:") :].strip()
        elif not query_name and line.startswith("Query="):
            query_name = line[len("Query=") :].strip()
        if not reference_name and line.startswith(">"):
            reference_name = line[1:].strip()
        if query_name and reference_name:
            break

    if not query_name:
        raise ValueError("Invalid BLAST output: missing query header line ('Query:' or 'Query=')")
    if not reference_name:
        raise ValueError("Invalid BLAST output: missing subject header line beginning with '>'")

    score_indices: List[Tuple[int, float]] = []
    for idx, raw_line in enumerate(lines):
        match = _BLAST_SCORE_RE.search(raw_line)
        if match:
            score_indices.append((idx, float(match.group(1))))

    if not score_indices:
        raise ValueError("Invalid BLAST output: no HSP score lines found")

    best_start_idx, _ = max(score_indices, key=lambda item: item[1])
    sorted_score_starts = sorted(idx for idx, _ in score_indices)
    best_pos = sorted_score_starts.index(best_start_idx)
    best_end_idx = (
        sorted_score_starts[best_pos + 1] if best_pos + 1 < len(sorted_score_starts) else len(lines)
    )

    query_chunks: List[str] = []
    reference_chunks: List[str] = []
    idx = best_start_idx + 1
    while idx < best_end_idx:
        line = lines[idx].strip()
        if not line:
            idx += 1
            continue
        if line.startswith("Query "):
            query_segment = _parse_alignment_row(line, "Query")
            sbjct_line: str | None = None
            scan = idx + 1
            while scan < best_end_idx:
                candidate = lines[scan].strip()
                if not candidate:
                    scan += 1
                    continue
                if candidate.startswith("Sbjct "):
                    sbjct_line = candidate
                    break
                if candidate.startswith("Query ") or candidate.startswith("Score:") or candidate.startswith(">"):
                    break
                scan += 1
            if sbjct_line is None:
                raise ValueError("Invalid BLAST output: Query row without a matching Sbjct row in selected HSP")
            reference_segment = _parse_alignment_row(sbjct_line, "Sbjct")
            if len(query_segment) != len(reference_segment):
                raise ValueError(
                    "Invalid BLAST output: Query/Sbjct segment length mismatch in selected HSP"
                )
            expanded_reference = []
            for q_char, r_char in zip(query_segment, reference_segment):
                if r_char == ".":
                    expanded_reference.append(q_char)
                else:
                    expanded_reference.append(r_char)
            query_chunks.append(query_segment.upper())
            reference_chunks.append("".join(expanded_reference).upper())
            idx = scan + 1
            continue
        idx += 1

    if not query_chunks or not reference_chunks:
        raise ValueError("Invalid BLAST output: selected top-scoring HSP contains no Query/Sbjct alignment rows")

    query = "".join(query_chunks)
    reference = "".join(reference_chunks)
    if len(query) != len(reference):
        raise ValueError("Invalid BLAST output: reconstructed query/reference alignment lengths differ")

    return query, reference, query_name, reference_name


def parse_alignment_pair(path: Path) -> Tuple[str, str, str, str]:
    first_line = _first_nonempty_line(path)
    fasta_error: Exception | None = None
    blast_error: Exception | None = None

    if first_line.startswith(">"):
        try:
            return parse_fasta_pair(path)
        except Exception as exc:
            fasta_error = exc
        try:
            return parse_blast_pairwise(path)
        except Exception as exc:
            blast_error = exc
    else:
        try:
            return parse_blast_pairwise(path)
        except Exception as exc:
            blast_error = exc

    details: List[str] = []
    if fasta_error is not None:
        details.append(f"FASTA parse failed: {fasta_error}")
    if blast_error is not None:
        details.append(f"BLAST parse failed: {blast_error}")
    suffix = f" Details: {' | '.join(details)}" if details else ""
    raise ValueError(
        "Unsupported alignment input format. Expected either a 2-sequence aligned FASTA "
        "or standard BLAST pairwise output with 'Query'/'Sbjct' rows."
        f"{suffix}"
    )


def parse_annotation_file(path: Path | None):
    from core.engine import parse_annotation_file as _parse_annotation_file

    return _parse_annotation_file(path)
