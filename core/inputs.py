from __future__ import annotations

from pathlib import Path
from typing import List, Tuple


def parse_fasta_pair(path: Path) -> Tuple[str, str, str, str]:
    from align_viz import parse_fasta_pair as _parse_fasta_pair

    return _parse_fasta_pair(path)


def parse_annotation_file(path: Path | None):
    from align_viz import parse_annotation_file as _parse_annotation_file

    return _parse_annotation_file(path)
