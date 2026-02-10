from __future__ import annotations

from pathlib import Path

from flask import Flask, jsonify, make_response, request, send_from_directory

from core.service import (
    export_alignment,
    get_session,
    normalize_viewport,
    probe_alignment,
    render_alignment,
    session_from_payload,
)

BASE_DIR = Path(__file__).resolve().parent
STATIC_DIR = BASE_DIR / "static"

app = Flask(__name__, static_folder=str(STATIC_DIR), static_url_path="/static")


def _parse_viewport(payload: dict, extent: float):
    start = payload.get("viewport_start")
    end = payload.get("viewport_end")
    if start is None or end is None:
        return None
    return normalize_viewport((float(start), float(end)), extent)


@app.get("/")
def index():
    return send_from_directory(STATIC_DIR, "index.html")


@app.post("/api/render")
def api_render():
    payload = request.get_json(silent=True) or {}
    try:
        session = session_from_payload(payload)
        viewport = _parse_viewport(payload, session.global_extent)
        result = render_alignment(session, viewport=viewport)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400

    return jsonify(
        {
            "token": result.token,
            "svg": result.svg,
            "global_extent": result.global_extent,
            "viewport_start": result.viewport_start,
            "viewport_end": result.viewport_end,
            "query_name": session.data.query_name,
            "reference_name": session.data.reference_name,
        }
    )


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
        viewport = _parse_viewport(payload, session.global_extent)
        blob = export_alignment(session, fmt, viewport=viewport)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400

    mime = "image/svg+xml" if fmt == "svg" else "image/png"
    filename = f"alignment_export.{fmt}"

    response = make_response(blob)
    response.headers["Content-Type"] = mime
    response.headers["Content-Disposition"] = f"attachment; filename={filename}"
    return response


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True)
