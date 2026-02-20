from __future__ import annotations

from pathlib import Path

from flask import Flask, jsonify, make_response, request, send_from_directory

from core.genome_workflow import (
    cancel_job,
    cleanup_expired_records,
    create_upload_from_streams,
    get_dotplot,
    get_job,
    render_alignment_job_for_viewer,
    start_alignment_job,
    start_dotplot_job,
)
from core.mpl_backend import configure_headless_matplotlib
from core.service import (
    extract_sequence_range,
    export_alignment,
    get_session,
    probe_alignment,
    recommend_render_preset,
    render_alignment,
    session_from_payload,
)

configure_headless_matplotlib()

BASE_DIR = Path(__file__).resolve().parent
STATIC_DIR = BASE_DIR / "static"

app = Flask(__name__, static_folder=str(STATIC_DIR), static_url_path="/static")


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


@app.get("/")
def index():
    return send_from_directory(STATIC_DIR, "index.html")


@app.post("/api/render")
def api_render():
    cleanup_expired_records()
    payload = request.get_json(silent=True) or {}
    try:
        session = session_from_payload(payload)
        result = render_alignment(session, viewport=None)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify(
        {
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
    )


@app.post("/api/preset_recommendation")
def api_preset_recommendation():
    cleanup_expired_records()
    payload = request.get_json(silent=True) or {}
    input_text = str(payload.get("input_path", "")).strip()
    if not input_text:
        return jsonify({"error": "input_path is required"}), 400

    try:
        recommendation = recommend_render_preset(
            input_path=Path(input_text),
            swap_roles=_to_bool(payload.get("swap_roles"), default=False, name="swap_roles"),
        )
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify(recommendation)


@app.post("/api/probe")
def api_probe():
    cleanup_expired_records()
    payload = request.get_json(silent=True) or {}
    token = str(payload.get("token", "")).strip()
    x_coord = payload.get("x_coord")

    if not token:
        return jsonify({"error": "token is required"}), 400
    if x_coord is None:
        return jsonify({"error": "x_coord is required"}), 400

    try:
        session = get_session(token)
        probe = probe_alignment(session, float(x_coord))
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify(probe)


@app.post("/api/export")
def api_export():
    cleanup_expired_records()
    payload = request.get_json(silent=True) or {}
    token = str(payload.get("token", "")).strip()
    fmt = str(payload.get("format", "svg")).strip().lower()

    if not token:
        return jsonify({"error": "token is required"}), 400

    try:
        session = get_session(token)
        blob = export_alignment(session, fmt, viewport=None)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400

    mime = "image/svg+xml" if fmt == "svg" else "image/png"
    filename = f"alignment_export.{fmt}"

    response = make_response(blob)
    response.headers["Content-Type"] = mime
    response.headers["Content-Disposition"] = f"attachment; filename={filename}"
    return response


@app.post("/api/extract_sequence")
def api_extract_sequence():
    cleanup_expired_records()
    payload = request.get_json(silent=True) or {}
    token = str(payload.get("token", "")).strip()
    stream = str(payload.get("stream", "")).strip().lower()
    start_x = payload.get("start_x")
    end_x = payload.get("end_x")

    if not token:
        return jsonify({"error": "token is required"}), 400
    if stream not in {"query", "reference"}:
        return jsonify({"error": "stream must be 'query' or 'reference'"}), 400
    if start_x is None:
        return jsonify({"error": "start_x is required"}), 400
    if end_x is None:
        return jsonify({"error": "end_x is required"}), 400

    try:
        session = get_session(token)
        result = extract_sequence_range(
            session,
            start_x=float(start_x),
            end_x=float(end_x),
            stream=stream,
        )
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify(result)


@app.post("/api/genome/upload")
def api_genome_upload():
    cleanup_expired_records()
    query_file = request.files.get("query_fasta")
    reference_file = request.files.get("reference_fasta")
    if query_file is None or reference_file is None:
        return jsonify({"error": "query_fasta and reference_fasta are required"}), 400

    try:
        payload = create_upload_from_streams(
            query_filename=query_file.filename or "query.fa",
            query_stream=query_file.stream,
            reference_filename=reference_file.filename or "reference.fa",
            reference_stream=reference_file.stream,
        )
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify(payload)


@app.post("/api/genome/dotplot/start")
def api_genome_dotplot_start():
    cleanup_expired_records()
    payload = request.get_json(silent=True) or {}
    upload_id = str(payload.get("upload_id", "")).strip()
    nucmer_options = payload.get("nucmer_options", {})
    if not upload_id:
        return jsonify({"error": "upload_id is required"}), 400

    try:
        job_id = start_dotplot_job(upload_id, nucmer_options=nucmer_options)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify({"job_id": job_id})


@app.get("/api/genome/jobs/<job_id>")
def api_genome_job_status(job_id: str):
    cleanup_expired_records()
    try:
        status = get_job(job_id)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify(status)


@app.post("/api/genome/jobs/<job_id>/cancel")
def api_genome_job_cancel(job_id: str):
    cleanup_expired_records()
    try:
        payload = cancel_job(job_id)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify(payload)


@app.get("/api/genome/dotplot/<upload_id>")
def api_genome_dotplot(upload_id: str):
    cleanup_expired_records()
    try:
        payload = get_dotplot(upload_id)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify(payload)


@app.post("/api/genome/align/start")
def api_genome_align_start():
    cleanup_expired_records()
    payload = request.get_json(silent=True) or {}
    upload_id = str(payload.get("upload_id", "")).strip()
    selected_blocks = payload.get("selected_blocks")
    if not upload_id:
        return jsonify({"error": "upload_id is required"}), 400
    if not isinstance(selected_blocks, list):
        return jsonify({"error": "selected_blocks must be a list"}), 400

    try:
        job_id = start_alignment_job(upload_id, selected_blocks)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify({"job_id": job_id})


@app.post("/api/genome/send_to_viewer")
def api_genome_send_to_viewer():
    cleanup_expired_records()
    payload = request.get_json(silent=True) or {}
    alignment_job_id = str(payload.get("alignment_job_id", "")).strip()
    if not alignment_job_id:
        return jsonify({"error": "alignment_job_id is required"}), 400

    try:
        viewer_payload = render_alignment_job_for_viewer(
            alignment_job_id,
            params_payload=payload.get("params", {}),
            query_annotation_path=str(payload.get("query_annotation_path", "")).strip() or None,
            reference_annotation_path=str(payload.get("reference_annotation_path", "")).strip() or None,
            swap_roles=_to_bool(payload.get("swap_roles"), default=False, name="swap_roles"),
        )
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify(viewer_payload)


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True, threaded=False)
