from __future__ import annotations

import threading
import unittest
from pathlib import Path

import matplotlib

from core.mpl_backend import configure_headless_matplotlib
from core.inputs import parse_alignment_pair
from core.params import RenderParams
from core.service import prepare_session, probe_alignment, render_alignment
from webapp.app import app


ROOT = Path(__file__).resolve().parents[1]
ALIGN = ROOT / "test_alignment.fa"
REF_ANNOT = ROOT / "test_ref_annotation.txt"
BLAST_ALIGN = ROOT / "tests" / "fixtures" / "blast_pairwise_standard.txt"
BLAST_MULTI_HSP = ROOT / "tests" / "fixtures" / "blast_pairwise_multi_hsp.txt"
BLAST_INVALID = ROOT / "tests" / "fixtures" / "blast_pairwise_invalid.txt"


class CoreTests(unittest.TestCase):
    def test_backend_configuration_is_agg_and_idempotent(self):
        configure_headless_matplotlib()
        configure_headless_matplotlib()
        backend = str(matplotlib.get_backend()).lower()
        self.assertIn("agg", backend)

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

    def test_render_in_worker_thread(self):
        params = RenderParams()
        session = prepare_session(
            input_path=ALIGN,
            params=params,
            reference_annotation_path=REF_ANNOT,
        )

        outcome = {"error": None, "ok": False}

        def _run() -> None:
            try:
                result = render_alignment(session, viewport=(0, 2000))
                outcome["ok"] = "<svg" in result.svg
            except Exception as exc:  # pragma: no cover - test assertion path
                outcome["error"] = str(exc)

        worker = threading.Thread(target=_run)
        worker.start()
        worker.join(timeout=60)

        self.assertFalse(worker.is_alive(), "worker thread did not finish")
        self.assertIsNone(outcome["error"], f"thread render failed: {outcome['error']}")
        self.assertTrue(outcome["ok"])

    def test_parse_alignment_pair_fasta_still_supported(self):
        query, reference, query_name, reference_name = parse_alignment_pair(ALIGN)
        self.assertEqual(len(query), len(reference))
        self.assertTrue(query_name)
        self.assertTrue(reference_name)

    def test_parse_alignment_pair_blast_dot_identity_expansion(self):
        query, reference, query_name, reference_name = parse_alignment_pair(BLAST_ALIGN)
        self.assertEqual(len(query), len(reference))
        self.assertNotIn(".", reference)
        self.assertIn("-", query)
        self.assertIn("-", reference)
        self.assertEqual(query_name, "ExampleQuery")
        self.assertEqual(reference_name, "ExampleReference")

    def test_blast_multi_hsp_uses_top_score(self):
        query, reference, _, _ = parse_alignment_pair(BLAST_MULTI_HSP)
        self.assertEqual(query, "CCCC")
        self.assertEqual(reference, "TCCC")


class ApiTests(unittest.TestCase):
    def setUp(self):
        self.client = app.test_client()

    def test_render_probe_export(self):
        payload = {
            "input_path": str(ALIGN),
            "reference_annotation_path": str(REF_ANNOT),
            "query_annotation_path": "",
            "params": {"width": "auto", "height": 6.0, "dpi": 100},
        }

        render_resp = self.client.post("/api/render", json=payload)
        self.assertEqual(render_resp.status_code, 200)
        render_json = render_resp.get_json()
        token = render_json["token"]
        self.assertIn("svg", render_json)
        self.assertGreater(render_json["x_data_max"], render_json["x_data_min"])
        self.assertGreater(render_json["axes_right_px"], render_json["axes_left_px"])
        self.assertGreater(render_json["svg_width_px"], 0)
        self.assertGreater(render_json["svg_height_px"], 0)

        probe_resp = self.client.post("/api/probe", json={"token": token, "x_coord": 100.0})
        self.assertEqual(probe_resp.status_code, 200)
        probe_json = probe_resp.get_json()
        self.assertIn("column_index", probe_json)

        export_resp = self.client.post(
            "/api/export",
            json={"token": token, "format": "svg"},
        )
        self.assertEqual(export_resp.status_code, 200)
        self.assertIn("image/svg+xml", export_resp.headers["Content-Type"])

    def test_api_render_accepts_blast_input(self):
        payload = {
            "input_path": str(BLAST_ALIGN),
            "reference_annotation_path": "",
            "query_annotation_path": "",
            "params": {"width": "auto", "height": 6.0, "dpi": 100},
        }

        render_resp = self.client.post("/api/render", json=payload)
        self.assertEqual(render_resp.status_code, 200)
        render_json = render_resp.get_json()
        self.assertIn("<svg", render_json["svg"])
        self.assertEqual(render_json["query_name"], "ExampleQuery")
        self.assertEqual(render_json["reference_name"], "ExampleReference")

    def test_blast_invalid_format_errors_cleanly(self):
        payload = {
            "input_path": str(BLAST_INVALID),
            "reference_annotation_path": "",
            "query_annotation_path": "",
            "params": {"width": "auto", "height": 6.0, "dpi": 100},
        }
        render_resp = self.client.post("/api/render", json=payload)
        self.assertEqual(render_resp.status_code, 400)
        body = render_resp.get_json()
        self.assertIn("BLAST", body["error"])


if __name__ == "__main__":
    unittest.main()
