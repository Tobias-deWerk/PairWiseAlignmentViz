from __future__ import annotations

from pathlib import Path
from typing import Tuple


def parse_fasta_pair(path: Path) -> Tuple[str, str, str, str]:
    from core.engine import parse_fasta_pair as _parse_fasta_pair

    return _parse_fasta_pair(path)


def parse_annotation_file(path: Path | None):
    from core.engine import parse_annotation_file as _parse_annotation_file

    return _parse_annotation_file(path)
