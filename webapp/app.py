from __future__ import annotations

from pathlib import Path

from flask import Flask, jsonify, make_response, request, send_from_directory

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
        }
    )


@app.post("/api/preset_recommendation")
def api_preset_recommendation():
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


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True, threaded=False)
