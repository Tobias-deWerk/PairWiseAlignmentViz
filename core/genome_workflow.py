from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import BinaryIO, Dict, List, Optional, Sequence, Tuple

from core.params import RenderParams
from core.service import RenderSession, prepare_session, render_alignment


MAX_UPLOAD_BYTES = 250 * 1024 * 1024
MIN_GENOME_WARN_BP = 100_000
WORK_ROOT = Path(tempfile.gettempdir()) / "pwav_genome_jobs"
UPLOAD_TTL_SECONDS = 6 * 60 * 60
JOB_TTL_SECONDS = 6 * 60 * 60


@dataclass
class UploadRecord:
    upload_id: str
    dir_path: Path
    query_path: Path
    reference_path: Path
    query_name: str
    reference_name: str
    query_sequence: str
    reference_sequence: str
    query_length: int
    reference_length: int
    warnings: List[str] = field(default_factory=list)
    dotplot: Optional[Dict[str, object]] = None
    created_at: float = field(default_factory=time.time)


@dataclass
class JobRecord:
    job_id: str
    job_type: str
    upload_id: str
    status: str = "queued"
    progress: float = 0.0
    stage: str = "queued"
    message: str = ""
    error: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    result: Dict[str, object] = field(default_factory=dict)
    cancel_requested: bool = False
    created_at: float = field(default_factory=time.time)
    updated_at: float = field(default_factory=time.time)


UPLOADS: Dict[str, UploadRecord] = {}
JOBS: Dict[str, JobRecord] = {}
_LOCK = threading.RLock()


_COMPLEMENT = str.maketrans({"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"})
NUCMER_MATCH_MODES = {"maxmatch", "mumreference", "mum"}
DELTA_FILTER_MODES = {
    "none": None,
    "one_to_one": "-1",
    "many_to_many": "-m",
    "best_query": "-q",
    "best_reference": "-r",
}


def _ensure_work_root() -> None:
    WORK_ROOT.mkdir(parents=True, exist_ok=True)


def cleanup_expired_records() -> None:
    now = time.time()
    with _LOCK:
        upload_ids_to_delete = [
            upload_id
            for upload_id, record in UPLOADS.items()
            if now - record.created_at > UPLOAD_TTL_SECONDS
        ]
        for upload_id in upload_ids_to_delete:
            record = UPLOADS.pop(upload_id)
            shutil.rmtree(record.dir_path, ignore_errors=True)

        job_ids_to_delete = [
            job_id
            for job_id, record in JOBS.items()
            if now - record.created_at > JOB_TTL_SECONDS
        ]
        for job_id in job_ids_to_delete:
            JOBS.pop(job_id, None)


def _safe_name(name: str, fallback: str) -> str:
    text = (name or "").strip()
    if not text:
        return fallback
    return "".join(ch for ch in text if ch.isalnum() or ch in {"_", "-", "."}) or fallback


def _read_single_fasta(path: Path) -> Tuple[str, str, int]:
    header: Optional[str] = None
    chunks: List[str] = []
    records = 0

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                records += 1
                if records == 1:
                    header = line[1:].strip() or "sequence"
                elif records == 2:
                    break
                continue
            if records == 0:
                raise ValueError(f"FASTA file '{path.name}' must start with a header line beginning with '>'")
            if records == 1:
                chunks.append(line.upper())

    if records == 0:
        raise ValueError(f"FASTA file '{path.name}' does not contain any FASTA records")

    sequence = "".join(chunks).replace(" ", "")
    if not sequence:
        raise ValueError(f"FASTA file '{path.name}' contains an empty first sequence")
    return header or "sequence", sequence, records


def _save_upload_file(upload_dir: Path, filename: str, stream: BinaryIO, target_name: str) -> Path:
    safe = _safe_name(filename, target_name)
    suffix = Path(safe).suffix or ".fa"
    target = upload_dir / f"{target_name}{suffix}"
    total = 0
    with target.open("wb") as out:
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            total += len(chunk)
            if total > MAX_UPLOAD_BYTES:
                raise ValueError(f"Uploaded file '{filename}' exceeds {MAX_UPLOAD_BYTES} bytes")
            out.write(chunk)
    if total == 0:
        raise ValueError(f"Uploaded file '{filename}' is empty")
    return target


def create_upload_from_streams(
    query_filename: str,
    query_stream: BinaryIO,
    reference_filename: str,
    reference_stream: BinaryIO,
) -> Dict[str, object]:
    _ensure_work_root()
    cleanup_expired_records()

    upload_id = uuid.uuid4().hex
    upload_dir = WORK_ROOT / upload_id
    upload_dir.mkdir(parents=True, exist_ok=True)

    query_path = _save_upload_file(upload_dir, query_filename, query_stream, "query")
    reference_path = _save_upload_file(upload_dir, reference_filename, reference_stream, "reference")

    query_name, query_seq, query_records = _read_single_fasta(query_path)
    reference_name, reference_seq, reference_records = _read_single_fasta(reference_path)

    warnings: List[str] = []
    if query_records > 1:
        warnings.append("Query FASTA contains multiple records; only the first record is used.")
    if reference_records > 1:
        warnings.append("Reference FASTA contains multiple records; only the first record is used.")
    if max(len(query_seq), len(reference_seq)) < MIN_GENOME_WARN_BP:
        warnings.append(
            f"Input appears shorter than genome-scale ({MIN_GENOME_WARN_BP:,} bp). Workflow still allowed."
        )

    record = UploadRecord(
        upload_id=upload_id,
        dir_path=upload_dir,
        query_path=query_path,
        reference_path=reference_path,
        query_name=query_name,
        reference_name=reference_name,
        query_sequence=query_seq,
        reference_sequence=reference_seq,
        query_length=len(query_seq),
        reference_length=len(reference_seq),
        warnings=warnings,
    )
    with _LOCK:
        UPLOADS[upload_id] = record

    return {
        "upload_id": upload_id,
        "query_name": query_name,
        "reference_name": reference_name,
        "query_length": len(query_seq),
        "reference_length": len(reference_seq),
        "warnings": warnings,
    }


def _require_upload(upload_id: str) -> UploadRecord:
    with _LOCK:
        record = UPLOADS.get(upload_id)
    if record is None:
        raise ValueError("Unknown upload_id")
    return record


def _new_job(job_type: str, upload_id: str) -> JobRecord:
    job = JobRecord(job_id=uuid.uuid4().hex, job_type=job_type, upload_id=upload_id)
    with _LOCK:
        JOBS[job.job_id] = job
    return job


def _set_job_state(job: JobRecord, *, status: Optional[str] = None, progress: Optional[float] = None, stage: Optional[str] = None, message: Optional[str] = None, error: Optional[str] = None, result: Optional[Dict[str, object]] = None) -> None:
    with _LOCK:
        if status is not None:
            job.status = status
        if progress is not None:
            job.progress = float(max(0.0, min(1.0, progress)))
        if stage is not None:
            job.stage = stage
        if message is not None:
            job.message = message
        if error is not None:
            job.error = error
        if result is not None:
            job.result = result
        job.updated_at = time.time()


def _require_command(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise RuntimeError(f"Required command '{name}' is not installed or not on PATH")
    return path


def _run_command(args: Sequence[str], *, cwd: Path, timeout_seconds: int = 7200) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(args),
        cwd=str(cwd),
        check=True,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )


def _check_cancel(job: JobRecord) -> None:
    if job.cancel_requested:
        raise RuntimeError("Job cancelled")


def _to_positive_int(value: object, name: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if parsed <= 0:
        raise ValueError(f"{name} must be > 0")
    return parsed


def _to_positive_float(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a number") from exc
    if parsed <= 0:
        raise ValueError(f"{name} must be > 0")
    return parsed


def _to_bool(value: object, *, default: bool = False, name: str = "value") -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        if value in {0, 1}:
            return bool(value)
        raise ValueError(f"{name} must be a boolean")
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "y", "on"}:
            return True
        if normalized in {"false", "0", "no", "n", "off", ""}:
            return False
    raise ValueError(f"{name} must be a boolean")


def _parse_dotplot_options(payload_options: Optional[Dict[str, object]]) -> Dict[str, object]:
    options = payload_options or {}
    if not isinstance(options, dict):
        raise ValueError("nucmer_options must be an object")

    match_mode = str(options.get("match_mode", "maxmatch")).strip().lower() or "maxmatch"
    if match_mode not in NUCMER_MATCH_MODES:
        raise ValueError("nucmer_options.match_mode must be one of: maxmatch, mumreference, mum")

    parsed: Dict[str, object] = {"match_mode": match_mode}

    mincluster_raw = options.get("mincluster")
    if mincluster_raw not in (None, ""):
        parsed["mincluster"] = _to_positive_int(mincluster_raw, "nucmer_options.mincluster")

    diagfactor_raw = options.get("diagfactor")
    if diagfactor_raw not in (None, ""):
        parsed["diagfactor"] = _to_positive_float(diagfactor_raw, "nucmer_options.diagfactor")

    breaklen_raw = options.get("breaklen")
    if breaklen_raw not in (None, ""):
        parsed["breaklen"] = _to_positive_int(breaklen_raw, "nucmer_options.breaklen")

    delta_filter_mode = str(options.get("delta_filter_mode", "none")).strip().lower() or "none"
    if delta_filter_mode not in DELTA_FILTER_MODES:
        raise ValueError(
            "nucmer_options.delta_filter_mode must be one of: "
            "none, one_to_one, many_to_many, best_query, best_reference"
        )
    parsed["delta_filter_mode"] = delta_filter_mode

    return parsed


def _parse_align_options(payload_options: Optional[Dict[str, object]]) -> Dict[str, object]:
    options = payload_options or {}
    if not isinstance(options, dict):
        raise ValueError("align_options must be an object")
    include_intervals = _to_bool(
        options.get("include_inter_block_intervals"),
        default=False,
        name="align_options.include_inter_block_intervals",
    )
    return {"include_inter_block_intervals": include_intervals}


def _build_nucmer_command(
    nucmer_bin: str,
    *,
    options: Dict[str, object],
    prefix: Path,
    reference_path: Path,
    query_path: Path,
) -> List[str]:
    mode_flag = f"--{str(options.get('match_mode', 'maxmatch'))}"
    cmd: List[str] = [nucmer_bin, mode_flag]

    if "mincluster" in options:
        cmd.extend(["-c", str(options["mincluster"])])
    if "diagfactor" in options:
        cmd.extend(["-D", str(options["diagfactor"])])
    if "breaklen" in options:
        cmd.extend(["-b", str(options["breaklen"])])

    cmd.extend(["-p", str(prefix), str(reference_path), str(query_path)])
    return cmd


def _build_delta_filter_command(
    delta_filter_bin: str,
    *,
    mode: str,
    input_delta_path: Path,
) -> Optional[List[str]]:
    mode_flag = DELTA_FILTER_MODES.get(mode)
    if mode_flag is None:
        return None
    return [delta_filter_bin, mode_flag, str(input_delta_path)]


def _parse_show_coords(text: str) -> Dict[str, object]:
    points: List[Dict[str, object]] = []
    blocks: List[Dict[str, object]] = []
    block_idx = 0
    parsed_rows = 0
    skipped_rows = 0

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("=") or line.startswith("[") or line.startswith("NUCMER") or line.startswith("/"):
            continue

        # Prefer tabular parsing (show-coords -T), but remain compatible with legacy
        # output that may include visual separators like '|'.
        tab_parts = [part.strip() for part in line.split("\t")]
        if len(tab_parts) >= 7:
            parts = tab_parts
        else:
            parts = [part for part in line.replace("|", " ").split() if part]

        if len(parts) < 7:
            skipped_rows += 1
            continue

        try:
            r_start = int(parts[0])
            r_end = int(parts[1])
            q_start = int(parts[2])
            q_end = int(parts[3])
            len_r = int(parts[4])
            len_q = int(parts[5])
            identity = float(parts[6])
        except ValueError:
            skipped_rows += 1
            continue

        parsed_rows += 1
        orientation = "forward" if q_end >= q_start else "reverse"
        block_id = f"b{block_idx}"
        block_idx += 1

        block = {
            "block_id": block_id,
            "q_start": q_start,
            "q_end": q_end,
            "r_start": r_start,
            "r_end": r_end,
            "orientation": orientation,
            "score": identity,
            "length_query": len_q,
            "length_reference": len_r,
        }
        blocks.append(block)

        points.append({"q_pos": q_start, "r_pos": r_start, "orientation": orientation, "identity": identity})
        points.append({"q_pos": q_end, "r_pos": r_end, "orientation": orientation, "identity": identity})
        points.append(
            {
                "q_pos": int((q_start + q_end) / 2),
                "r_pos": int((r_start + r_end) / 2),
                "orientation": orientation,
                "identity": identity,
            }
        )

    return {
        "points": points,
        "blocks": blocks,
        "parsed_rows_count": parsed_rows,
        "skipped_rows_count": skipped_rows,
    }


def _run_dotplot_job(job: JobRecord) -> None:
    try:
        upload = _require_upload(job.upload_id)
        _set_job_state(job, status="running", stage="preflight", progress=0.05, message="Checking external tools")
        nucmer = _require_command("nucmer")
        show_coords = _require_command("show-coords")
        dotplot_options = _parse_dotplot_options(job.result.get("nucmer_options"))
        delta_filter_mode = str(dotplot_options.get("delta_filter_mode", "none"))
        delta_filter = _require_command("delta-filter") if delta_filter_mode != "none" else None

        _check_cancel(job)
        mode = str(dotplot_options.get("match_mode", "maxmatch"))
        mode_message = f"Running nucmer ({mode})"
        _set_job_state(job, stage="nucmer", progress=0.20, message=mode_message)
        prefix = upload.dir_path / "dotplot"
        nucmer_cmd = _build_nucmer_command(
            nucmer,
            options=dotplot_options,
            prefix=prefix,
            reference_path=upload.reference_path,
            query_path=upload.query_path,
        )
        _run_command(nucmer_cmd, cwd=upload.dir_path)

        _check_cancel(job)
        _set_job_state(job, stage="coords", progress=0.65, message="Extracting synteny blocks")
        delta_path = upload.dir_path / "dotplot.delta"
        coords_delta_path = delta_path
        delta_filter_cmd: Optional[List[str]] = None

        if delta_filter is not None:
            _set_job_state(job, stage="delta_filter", progress=0.58, message=f"Running delta-filter ({delta_filter_mode})")
            delta_filter_cmd = _build_delta_filter_command(
                delta_filter,
                mode=delta_filter_mode,
                input_delta_path=delta_path,
            )
            if delta_filter_cmd is None:
                raise RuntimeError("delta-filter command could not be constructed")
            filtered = _run_command(delta_filter_cmd, cwd=upload.dir_path)
            filtered_delta_path = upload.dir_path / "dotplot.filtered.delta"
            filtered_delta_path.write_text(filtered.stdout, encoding="utf-8")
            coords_delta_path = filtered_delta_path

        coords_cmd = [show_coords, "-T", "-rcl", str(coords_delta_path)]
        coords = _run_command(coords_cmd, cwd=upload.dir_path)

        dotplot = _parse_show_coords(coords.stdout)
        with _LOCK:
            upload.dotplot = dotplot

        dotplot_json = upload.dir_path / "dotplot.json"
        dotplot_json.write_text(json.dumps(dotplot, indent=2), encoding="utf-8")

        blocks_count = len(dotplot["blocks"])
        points_count = len(dotplot["points"])
        parsed_rows_count = int(dotplot.get("parsed_rows_count", 0))
        skipped_rows_count = int(dotplot.get("skipped_rows_count", 0))
        done_message = "Dotplot complete"
        warnings: List[str] = []
        if blocks_count == 0:
            warning = "No synteny blocks found. Try switching match mode or lowering mincluster."
            if parsed_rows_count == 0:
                warning += " If this persists, inspect parser/format compatibility."
            warnings.append(warning)
            done_message = warning

        _set_job_state(
            job,
            status="done",
            stage="done",
            progress=1.0,
            message=done_message,
            result={
                **dict(job.result),
                "upload_id": upload.upload_id,
                "dotplot_path": str(dotplot_json),
                "delta_path": str(delta_path),
                "coords_delta_path": str(coords_delta_path),
                "delta_filter_command": delta_filter_cmd,
                "nucmer_command": nucmer_cmd,
                "coords_command": coords_cmd,
                "points_count": points_count,
                "blocks_count": blocks_count,
                "parsed_rows_count": parsed_rows_count,
                "skipped_rows_count": skipped_rows_count,
                "warnings": warnings,
            },
        )
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        text = stderr[-1500:] if stderr else str(exc)
        _set_job_state(
            job,
            status="failed",
            stage="failed",
            error=f"Dotplot command failed: {text}",
            result={
                **dict(job.result),
                "stderr_tail": text,
            },
        )
    except Exception as exc:  # pragma: no cover
        message = str(exc)
        if message == "Job cancelled":
            _set_job_state(job, status="cancelled", stage="cancelled", message="Job cancelled")
        else:
            _set_job_state(job, status="failed", stage="failed", error=message)


def start_dotplot_job(upload_id: str, nucmer_options: Optional[Dict[str, object]] = None) -> str:
    _require_upload(upload_id)
    parsed_options = _parse_dotplot_options(nucmer_options)
    job = _new_job("dotplot", upload_id)
    job.result = {"nucmer_options": parsed_options}
    worker = threading.Thread(target=_run_dotplot_job, args=(job,), daemon=True)
    worker.start()
    return job.job_id


def get_job(job_id: str) -> Dict[str, object]:
    cleanup_expired_records()
    with _LOCK:
        job = JOBS.get(job_id)
        if job is None:
            raise ValueError("Unknown job_id")
        return {
            "job_id": job.job_id,
            "job_type": job.job_type,
            "upload_id": job.upload_id,
            "status": job.status,
            "progress": job.progress,
            "stage": job.stage,
            "message": job.message,
            "error": job.error,
            "warnings": list(job.warnings),
            "result": dict(job.result),
            "created_at": job.created_at,
            "updated_at": job.updated_at,
            "cancel_requested": job.cancel_requested,
        }


def cancel_job(job_id: str) -> Dict[str, object]:
    with _LOCK:
        job = JOBS.get(job_id)
        if job is None:
            raise ValueError("Unknown job_id")
        job.cancel_requested = True
        job.updated_at = time.time()
    return {"job_id": job_id, "status": "cancelling"}


def get_dotplot(upload_id: str) -> Dict[str, object]:
    upload = _require_upload(upload_id)
    if upload.dotplot is None:
        raise ValueError("Dotplot not available yet for this upload")
    return {
        "upload_id": upload.upload_id,
        "query_name": upload.query_name,
        "reference_name": upload.reference_name,
        "query_length": upload.query_length,
        "reference_length": upload.reference_length,
        "warnings": list(upload.warnings),
        "dotplot": upload.dotplot,
    }


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]


def _find_block(blocks: Sequence[Dict[str, object]], block_id: str) -> Dict[str, object]:
    for block in blocks:
        if str(block.get("block_id")) == block_id:
            return block
    raise ValueError(f"Unknown block_id '{block_id}'")


def _normalize_bounds(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a <= b else (b, a)


def _trimmed_bounds(base_start: int, base_end: int, trim_start: Optional[int], trim_end: Optional[int]) -> Tuple[int, int]:
    low, high = _normalize_bounds(base_start, base_end)
    t0 = low if trim_start is None else max(low, min(high, int(trim_start)))
    t1 = high if trim_end is None else max(low, min(high, int(trim_end)))
    return _normalize_bounds(t0, t1)


def _extract_segment(seq: str, start_1: int, end_1: int) -> str:
    start = max(1, int(start_1)) - 1
    end = min(len(seq), int(end_1))
    if end <= start:
        return ""
    return seq[start:end]


def _extract_bridge_segment(seq: str, prev_anchor: int, next_anchor: int) -> str:
    start_anchor = int(prev_anchor)
    end_anchor = int(next_anchor)
    if abs(end_anchor - start_anchor) <= 1:
        return ""
    low = min(start_anchor, end_anchor) + 1
    high = max(start_anchor, end_anchor) - 1
    if high < low:
        return ""
    bridge = _extract_segment(seq, low, high)
    if not bridge:
        return ""
    if start_anchor > end_anchor:
        return reverse_complement(bridge)
    return bridge


def _ranges_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    left_a, right_a = _normalize_bounds(int(a_start), int(a_end))
    left_b, right_b = _normalize_bounds(int(b_start), int(b_end))
    return max(left_a, left_b) <= min(right_a, right_b)


def _map_labeled_local_ranges_to_columns(
    aligned_query: str,
    labeled_ranges: Sequence[Dict[str, object]],
) -> List[Dict[str, object]]:
    mapped: List[Dict[str, object]] = []
    for item in labeled_ranges:
        start_local = int(item.get("start_local", 0))
        end_local = int(item.get("end_local", 0))
        if start_local <= 0 or end_local <= 0:
            continue
        intervals = _map_local_ranges_to_columns(aligned_query, [(start_local, end_local)])
        for start_col, end_col in intervals:
            payload = dict(item)
            payload["start_column"] = int(start_col)
            payload["end_column"] = int(end_col)
            mapped.append(payload)
    return mapped


def _parse_mafft_pair(path: Path) -> Tuple[str, str]:
    seqs: List[str] = []
    chunks: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chunks:
                    seqs.append("".join(chunks).upper())
                    chunks = []
                continue
            chunks.append(line)
    if chunks:
        seqs.append("".join(chunks).upper())

    if len(seqs) != 2:
        raise ValueError(f"Expected two sequences in MAFFT output, found {len(seqs)}")
    if len(seqs[0]) != len(seqs[1]):
        raise ValueError("MAFFT output sequences differ in length")
    return seqs[0], seqs[1]


def _map_local_ranges_to_columns(aligned_query: str, local_ranges: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    result: List[Tuple[int, int]] = []
    if not local_ranges:
        return result

    current_local = 0
    interval_index = 0
    interval_start_col: Optional[int] = None
    sorted_ranges = sorted(local_ranges, key=lambda pair: pair[0])

    for col, base in enumerate(aligned_query, start=1):
        if base != "-":
            current_local += 1

        if interval_index >= len(sorted_ranges):
            break
        start_local, end_local = sorted_ranges[interval_index]

        if current_local >= start_local and current_local <= end_local and base != "-":
            if interval_start_col is None:
                interval_start_col = col

        if current_local > end_local:
            if interval_start_col is not None:
                result.append((interval_start_col, col - 1))
            interval_start_col = None
            interval_index += 1
            while interval_index < len(sorted_ranges) and sorted_ranges[interval_index][1] < current_local:
                interval_index += 1
            if interval_index < len(sorted_ranges):
                start_local, end_local = sorted_ranges[interval_index]
                if current_local >= start_local and current_local <= end_local and base != "-":
                    interval_start_col = col

    if interval_start_col is not None and interval_index < len(sorted_ranges):
        result.append((interval_start_col, len(aligned_query)))
    return result


def _build_selected_sequences(
    upload: UploadRecord,
    selected_blocks: Sequence[Dict[str, object]],
    *,
    include_inter_block_intervals: bool = False,
) -> Tuple[
    str,
    str,
    List[Tuple[int, int]],
    List[Dict[str, object]],
    List[Dict[str, object]],
    int,
    List[str],
]:
    blocks = [] if upload.dotplot is None else list(upload.dotplot.get("blocks", []))

    if not selected_blocks:
        return upload.query_sequence, upload.reference_sequence, [], [], [], 0, []

    query_parts: List[str] = []
    reference_parts: List[str] = []
    inversion_ranges_local: List[Tuple[int, int]] = []
    duplication_local_ranges: List[Dict[str, object]] = []
    selected_summary: List[Dict[str, object]] = []
    stitching_notes: List[str] = []
    inter_block_intervals_count = 0
    query_offset = 0
    previous_block: Optional[Dict[str, object]] = None

    for item in selected_blocks:
        include = bool(item.get("include", True))
        if not include:
            continue

        block_id = str(item.get("block_id", "")).strip()
        if not block_id:
            raise ValueError("Each selected block must include block_id")
        block = _find_block(blocks, block_id)

        q0, q1 = _trimmed_bounds(
            int(block["q_start"]),
            int(block["q_end"]),
            item.get("trim_q_start"),
            item.get("trim_q_end"),
        )
        r0, r1 = _trimmed_bounds(
            int(block["r_start"]),
            int(block["r_end"]),
            item.get("trim_r_start"),
            item.get("trim_r_end"),
        )

        query_segment = _extract_segment(upload.query_sequence, q0, q1)
        reference_segment = _extract_segment(upload.reference_sequence, r0, r1)
        if not query_segment or not reference_segment:
            continue

        orientation = str(block.get("orientation", "forward"))
        query_anchor_start = q1 if orientation == "reverse" else q0
        query_anchor_end = q0 if orientation == "reverse" else q1
        reference_anchor_start = r0
        reference_anchor_end = r1

        if include_inter_block_intervals and previous_block is not None:
            query_bridge = _extract_bridge_segment(
                upload.query_sequence,
                int(previous_block["query_anchor_end"]),
                query_anchor_start,
            )
            reference_bridge = _extract_bridge_segment(
                upload.reference_sequence,
                int(previous_block["reference_anchor_end"]),
                reference_anchor_start,
            )
            if query_bridge:
                query_parts.append(query_bridge)
                query_offset += len(query_bridge)
            if reference_bridge:
                reference_parts.append(reference_bridge)
            if query_bridge or reference_bridge:
                inter_block_intervals_count += 1
            if bool(query_bridge) != bool(reference_bridge):
                stitching_notes.append(
                    "Inter-block bridge between "
                    f"{previous_block['block_id']} and {block_id} produced only one stream interval."
                )

        if orientation == "reverse":
            query_segment = reverse_complement(query_segment)

        start_local = query_offset + 1
        end_local = query_offset + len(query_segment)
        if orientation == "reverse":
            inversion_ranges_local.append((start_local, end_local))

        query_parts.append(query_segment)
        reference_parts.append(reference_segment)
        query_offset += len(query_segment)

        selected_summary.append(
            {
                "block_id": block_id,
                "q_start": q0,
                "q_end": q1,
                "r_start": r0,
                "r_end": r1,
                "orientation": orientation,
                "query_local_start": start_local,
                "query_local_end": end_local,
            }
        )
        previous_block = {
            "block_id": block_id,
            "query_anchor_end": query_anchor_end,
            "reference_anchor_end": reference_anchor_end,
        }

    query_concat = "".join(query_parts)
    reference_concat = "".join(reference_parts)

    if not query_concat or not reference_concat:
        raise ValueError("Selection produced empty sequence set")

    duplicate_ids: set[str] = set()
    for left_idx in range(len(selected_summary)):
        left = selected_summary[left_idx]
        for right_idx in range(left_idx + 1, len(selected_summary)):
            right = selected_summary[right_idx]
            if _ranges_overlap(left["q_start"], left["q_end"], right["q_start"], right["q_end"]) or _ranges_overlap(
                left["r_start"], left["r_end"], right["r_start"], right["r_end"]
            ):
                duplicate_ids.add(str(left["block_id"]))
                duplicate_ids.add(str(right["block_id"]))

    for item in selected_summary:
        block_id = str(item["block_id"])
        is_duplicate = block_id in duplicate_ids
        item["is_duplicate"] = is_duplicate
        if not is_duplicate:
            continue
        duplication_local_ranges.append(
            {
                "start_local": int(item["query_local_start"]),
                "end_local": int(item["query_local_end"]),
                "label": "DUP",
                "block_id": block_id,
                "pointer_text": (
                    f"Duplicated on query {int(item['q_start'])}-{int(item['q_end'])}; "
                    f"reference {int(item['r_start'])}-{int(item['r_end'])}"
                ),
            }
        )

    return (
        query_concat,
        reference_concat,
        inversion_ranges_local,
        duplication_local_ranges,
        selected_summary,
        inter_block_intervals_count,
        stitching_notes,
    )


def _run_alignment_job(
    job: JobRecord,
    selected_blocks: Sequence[Dict[str, object]],
    align_options: Dict[str, object],
) -> None:
    try:
        upload = _require_upload(job.upload_id)
        parsed_align_options = _parse_align_options(align_options)
        include_intervals = bool(parsed_align_options.get("include_inter_block_intervals", False))
        _set_job_state(job, status="running", stage="preflight", progress=0.05, message="Checking MAFFT")
        mafft = _require_command("mafft")

        _check_cancel(job)
        _set_job_state(job, stage="selection", progress=0.20, message="Preparing selected segments")
        (
            query_selected,
            reference_selected,
            inversion_local_ranges,
            duplication_local_ranges,
            selected_summary,
            inter_block_intervals_count,
            stitching_notes,
        ) = _build_selected_sequences(
            upload,
            selected_blocks,
            include_inter_block_intervals=include_intervals,
        )

        work_dir = upload.dir_path / f"align_{job.job_id}"
        work_dir.mkdir(parents=True, exist_ok=True)

        mafft_input = work_dir / "selected_input.fa"
        mafft_input.write_text(
            "\n".join([">query_selected", query_selected, ">reference_selected", reference_selected, ""]),
            encoding="utf-8",
        )

        _check_cancel(job)
        _set_job_state(job, stage="mafft", progress=0.50, message="Running MAFFT global pairwise alignment")
        mafft_out = _run_command(
            [mafft, "--quiet", "--globalpair", "--maxiterate", "1000", str(mafft_input)],
            cwd=work_dir,
            timeout_seconds=3 * 60 * 60,
        )

        mafft_output_path = work_dir / "mafft_output.fa"
        mafft_output_path.write_text(mafft_out.stdout, encoding="utf-8")

        aligned_query, aligned_reference = _parse_mafft_pair(mafft_output_path)

        _set_job_state(job, stage="mapping", progress=0.85, message="Mapping inversion intervals")
        inversion_cols = _map_local_ranges_to_columns(aligned_query, inversion_local_ranges)
        inversion_intervals = [
            {"kind": "inversion", "start_column": int(start_col), "end_column": int(end_col), "label": "INV"}
            for start_col, end_col in inversion_cols
            if end_col >= start_col
        ]
        duplication_cols = _map_labeled_local_ranges_to_columns(aligned_query, duplication_local_ranges)
        duplication_intervals = [
            {
                "kind": "duplication",
                "start_column": int(item["start_column"]),
                "end_column": int(item["end_column"]),
                "label": str(item.get("label", "DUP")),
                "block_id": str(item.get("block_id", "")),
                "pointer_text": str(item.get("pointer_text", "")).strip(),
            }
            for item in duplication_cols
            if int(item.get("end_column", 0)) >= int(item.get("start_column", 0))
        ]

        aligned_path = work_dir / "aligned_pair.fa"
        aligned_path.write_text(
            "\n".join([f">{upload.query_name}", aligned_query, f">{upload.reference_name}", aligned_reference, ""]),
            encoding="utf-8",
        )

        done_message = "Alignment complete"
        if include_intervals:
            done_message += f" (inter-block intervals inserted: {inter_block_intervals_count})"

        _set_job_state(
            job,
            status="done",
            stage="done",
            progress=1.0,
            message=done_message,
            result={
                "align_options_used": parsed_align_options,
                "aligned_path": str(aligned_path),
                "alignment_length": len(aligned_query),
                "selected_blocks": selected_summary,
                "inversion_intervals": inversion_intervals,
                "duplication_intervals": duplication_intervals,
                "duplication_intervals_count": len(duplication_intervals),
                "inter_block_intervals_count": int(inter_block_intervals_count),
                "stitching_notes": stitching_notes,
            },
        )
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        text = stderr[-1500:] if stderr else str(exc)
        _set_job_state(job, status="failed", stage="failed", error=f"Alignment command failed: {text}")
    except Exception as exc:  # pragma: no cover
        message = str(exc)
        if message == "Job cancelled":
            _set_job_state(job, status="cancelled", stage="cancelled", message="Job cancelled")
        else:
            _set_job_state(job, status="failed", stage="failed", error=message)


def start_alignment_job(
    upload_id: str,
    selected_blocks: Sequence[Dict[str, object]],
    *,
    align_options: Optional[Dict[str, object]] = None,
) -> str:
    upload = _require_upload(upload_id)
    if upload.dotplot is None:
        raise ValueError("Dotplot must be completed before starting alignment")
    if selected_blocks is None:
        raise ValueError("selected_blocks is required")

    if not any(bool(item.get("include", True)) for item in selected_blocks):
        raise ValueError("At least one selected block must be included")

    parsed_align_options = _parse_align_options(align_options)
    job = _new_job("align", upload_id)
    job.result = {"align_options_used": parsed_align_options}
    worker = threading.Thread(
        target=_run_alignment_job,
        args=(job, list(selected_blocks), parsed_align_options),
        daemon=True,
    )
    worker.start()
    return job.job_id


def _set_feature_regions(
    session: RenderSession,
    inversion_intervals: Sequence[Dict[str, object]],
    duplication_intervals: Sequence[Dict[str, object]],
) -> None:
    regions: List[Dict[str, object]] = []
    all_intervals = [
        *[("inversion", interval) for interval in inversion_intervals],
        *[("duplication", interval) for interval in duplication_intervals],
    ]
    for kind, interval in all_intervals:
        start_col = int(interval.get("start_column", 0))
        end_col = int(interval.get("end_column", 0))
        if start_col <= 0 or end_col <= 0:
            continue
        if end_col < start_col:
            start_col, end_col = end_col, start_col
        if start_col > len(session.global_x):
            continue
        end_col = min(end_col, len(session.global_x))
        start_idx = start_col - 1
        end_idx = end_col - 1
        tag = str(interval.get("label", "INV" if kind == "inversion" else "DUP"))
        region = {
            "kind": kind,
            "start_column": start_col,
            "end_column": end_col,
            "start_x": float(session.global_x[start_idx]),
            "end_x": float(session.global_x[end_idx]),
            "tag": tag,
        }
        pointer_text = str(interval.get("pointer_text", "")).strip()
        if pointer_text:
            region["pointer_text"] = pointer_text
        regions.append(
            region
        )
    session.feature_regions = regions
    session.inversion_regions = [dict(region) for region in regions if str(region.get("kind")) == "inversion"]


def render_alignment_job_for_viewer(
    alignment_job_id: str,
    *,
    params_payload: Optional[Dict[str, object]] = None,
    query_annotation_path: Optional[str] = None,
    reference_annotation_path: Optional[str] = None,
    swap_roles: bool = False,
) -> Dict[str, object]:
    with _LOCK:
        job = JOBS.get(alignment_job_id)
    if job is None:
        raise ValueError("Unknown alignment_job_id")
    if job.job_type != "align":
        raise ValueError("Provided job is not an alignment job")
    if job.status != "done":
        raise ValueError("Alignment job is not finished")

    aligned_path_text = str(job.result.get("aligned_path", "")).strip()
    if not aligned_path_text:
        raise ValueError("Alignment output path is missing")
    aligned_path = Path(aligned_path_text)
    if not aligned_path.exists():
        raise FileNotFoundError("Alignment output file no longer exists")

    params = RenderParams.from_payload(params_payload or {})
    query_ann = Path(query_annotation_path).expanduser() if query_annotation_path else None
    ref_ann = Path(reference_annotation_path).expanduser() if reference_annotation_path else None

    session = prepare_session(
        input_path=aligned_path,
        params=params,
        query_annotation_path=query_ann,
        reference_annotation_path=ref_ann,
        swap_roles=swap_roles,
    )

    inversion_intervals = list(job.result.get("inversion_intervals", []))
    duplication_intervals = list(job.result.get("duplication_intervals", []))
    _set_feature_regions(session, inversion_intervals, duplication_intervals)

    result = render_alignment(session, viewport=None)
    return {
        "token": result.token,
        "svg": result.svg,
        "global_extent": result.global_extent,
        "x_data_min": result.x_data_min,
        "x_data_max": result.x_data_max,
        "axes_left_px": result.axes_left_px,
        "axes_right_px": result.axes_right_px,
        "svg_width_px": result.svg_width_px,
        "svg_height_px": result.svg_height_px,
        "query_name": session.data.query_name,
        "reference_name": session.data.reference_name,
        "inversion_regions": result.inversion_regions,
        "feature_tags": result.feature_tags,
    }
