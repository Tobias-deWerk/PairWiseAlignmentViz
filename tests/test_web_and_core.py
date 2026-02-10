from __future__ import annotations

import unittest
from pathlib import Path

from core.params import RenderParams
from core.service import prepare_session, probe_alignment, render_alignment
from webapp.app import app


ROOT = Path(__file__).resolve().parents[1]
ALIGN = ROOT / "test_alignment.fa"
REF_ANNOT = ROOT / "test_ref_annotation.txt"


class CoreTests(unittest.TestCase):
    def test_params_payload_parsing(self):
        params = RenderParams.from_payload(
            {
                "width": "auto",
                "height": "6.0",
                "dpi": "150",
                "gap_label_size": "NULL",
                "annotation_label_size": "NA",
            }
        )
        self.assertEqual(params.width, "auto")
        self.assertIsNone(params.gap_label_size)
        self.assertIsNone(params.annotation_label_size)

    def test_render_and_probe(self):
        params = RenderParams()
        session = prepare_session(
            input_path=ALIGN,
            params=params,
            reference_annotation_path=REF_ANNOT,
        )
        result = render_alignment(session, viewport=(0, 3000))
        self.assertIn("<svg", result.svg)

        probe = probe_alignment(session, 0.0)
        self.assertIn("query_base", probe)
        self.assertIn("reference_base", probe)


class ApiTests(unittest.TestCase):
    def setUp(self):
        self.client = app.test_client()

    def test_render_probe_export(self):
        payload = {
            "input_path": str(ALIGN),
            "reference_annotation_path": str(REF_ANNOT),
            "query_annotation_path": "",
            "viewport_start": 0,
            "viewport_end": 2000,
            "params": {"width": "auto", "height": 6.0, "dpi": 100},
        }

        render_resp = self.client.post("/api/render", json=payload)
        self.assertEqual(render_resp.status_code, 200)
        render_json = render_resp.get_json()
        token = render_json["token"]
        self.assertIn("svg", render_json)

        probe_resp = self.client.post("/api/probe", json={"token": token, "x_coord": 100.0})
        self.assertEqual(probe_resp.status_code, 200)
        probe_json = probe_resp.get_json()
        self.assertIn("column_index", probe_json)

        export_resp = self.client.post(
            "/api/export",
            json={"token": token, "format": "svg", "viewport_start": 0, "viewport_end": 1000},
        )
        self.assertEqual(export_resp.status_code, 200)
        self.assertIn("image/svg+xml", export_resp.headers["Content-Type"])


if __name__ == "__main__":
    unittest.main()
