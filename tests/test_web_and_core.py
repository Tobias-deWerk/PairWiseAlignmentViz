from __future__ import annotations

import threading
import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile

import matplotlib

from core.mpl_backend import configure_headless_matplotlib
from core.inputs import parse_alignment_pair
from core.params import RenderParams
from core.service import extract_sequence_range, prepare_session, probe_alignment, render_alignment
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

    def test_extract_sequence_range_normalizes_reversed_bounds(self):
        session = prepare_session(
            input_path=ALIGN,
            params=RenderParams(),
            reference_annotation_path=REF_ANNOT,
        )
        left_idx = 200
        right_idx = 320
        result = extract_sequence_range(
            session,
            start_x=float(session.global_x[right_idx]),
            end_x=float(session.global_x[left_idx]),
            stream="query",
        )
        expected = session.data.query[left_idx : right_idx + 1].replace("-", "")
        self.assertEqual(result["start_column_index"], left_idx)
        self.assertEqual(result["end_column_index"], right_idx)
        self.assertEqual(result["sequence"], expected)
        self.assertEqual(result["sequence_length"], len(expected))
        self.assertIn(">SL4_qPV11 selected_region", result["fasta"])

    def test_extract_sequence_range_all_gap_returns_none_locals(self):
        session = prepare_session(
            input_path=ALIGN,
            params=RenderParams(),
            reference_annotation_path=REF_ANNOT,
        )
        first_non_gap_idx = next(
            idx for idx, is_gap in enumerate(session.data.is_reference_gap) if not is_gap
        )
        result = extract_sequence_range(
            session,
            start_x=float(session.global_x[0]),
            end_x=float(session.global_x[first_non_gap_idx - 1]),
            stream="reference",
        )
        self.assertEqual(result["sequence"], "")
        self.assertEqual(result["sequence_length"], 0)
        self.assertIsNone(result["start_local_position"])
        self.assertIsNone(result["end_local_position"])
        self.assertIn("local=NA-NA", result["fasta"])

    def test_prepare_session_swap_roles_swaps_stream_identity(self):
        default_session = prepare_session(
            input_path=ALIGN,
            params=RenderParams(),
            reference_annotation_path=REF_ANNOT,
        )
        swapped_session = prepare_session(
            input_path=ALIGN,
            params=RenderParams(),
            reference_annotation_path=REF_ANNOT,
            swap_roles=True,
        )

        self.assertEqual(swapped_session.data.query_name, default_session.data.reference_name)
        self.assertEqual(swapped_session.data.reference_name, default_session.data.query_name)
        self.assertEqual(swapped_session.data.query, default_session.data.reference)
        self.assertEqual(swapped_session.data.reference, default_session.data.query)
        self.assertEqual(len(swapped_session.query_annotations), len(default_session.query_annotations))
        self.assertEqual(len(swapped_session.reference_annotations), len(default_session.reference_annotations))

    def test_extract_sequence_range_swapped_query_uses_swapped_stream(self):
        default_session = prepare_session(
            input_path=ALIGN,
            params=RenderParams(),
            reference_annotation_path=REF_ANNOT,
        )
        swapped_session = prepare_session(
            input_path=ALIGN,
            params=RenderParams(),
            reference_annotation_path=REF_ANNOT,
            swap_roles=True,
        )
        left_idx = 200
        right_idx = 320
        swapped_query = extract_sequence_range(
            swapped_session,
            start_x=float(swapped_session.global_x[left_idx]),
            end_x=float(swapped_session.global_x[right_idx]),
            stream="query",
        )
        default_reference = extract_sequence_range(
            default_session,
            start_x=float(default_session.global_x[left_idx]),
            end_x=float(default_session.global_x[right_idx]),
            stream="reference",
        )

        self.assertEqual(swapped_query["sequence"], default_reference["sequence"])
        self.assertIn(f">{swapped_session.data.query_name} selected_region", swapped_query["fasta"])

    def test_probe_alignment_swapped_roles_report_swapped_stream_data(self):
        default_session = prepare_session(
            input_path=ALIGN,
            params=RenderParams(),
            reference_annotation_path=REF_ANNOT,
        )
        swapped_session = prepare_session(
            input_path=ALIGN,
            params=RenderParams(),
            reference_annotation_path=REF_ANNOT,
            swap_roles=True,
        )
        idx = 250
        default_probe = probe_alignment(default_session, float(default_session.global_x[idx]))
        swapped_probe = probe_alignment(swapped_session, float(swapped_session.global_x[idx]))

        self.assertEqual(swapped_probe["query_base"], default_probe["reference_base"])
        self.assertEqual(swapped_probe["reference_base"], default_probe["query_base"])
        self.assertEqual(
            swapped_probe["query_local_position"],
            default_probe["reference_local_position"],
        )
        self.assertEqual(
            swapped_probe["reference_local_position"],
            default_probe["query_local_position"],
        )


class ApiTests(unittest.TestCase):
    def setUp(self):
        self.client = app.test_client()

    def _write_temp_alignment(self, length: int) -> Path:
        with NamedTemporaryFile("w", suffix=".fa", delete=False) as handle:
            query = "A" * length
            reference = "A" * length
            handle.write(">q\n")
            handle.write(query)
            handle.write("\n>r\n")
            handle.write(reference)
            handle.write("\n")
            return Path(handle.name)

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
        self.assertIn("query_base", probe_json)
        self.assertIn("reference_base", probe_json)

        export_resp = self.client.post(
            "/api/export",
            json={"token": token, "format": "svg"},
        )
        self.assertEqual(export_resp.status_code, 200)
        self.assertIn("image/svg+xml", export_resp.headers["Content-Type"])

        extract_query = self.client.post(
            "/api/extract_sequence",
            json={"token": token, "stream": "query", "start_x": 100.0, "end_x": 300.0},
        )
        self.assertEqual(extract_query.status_code, 200)
        extract_query_json = extract_query.get_json()
        self.assertEqual(extract_query_json["stream"], "query")
        self.assertIn("fasta", extract_query_json)

        extract_reference = self.client.post(
            "/api/extract_sequence",
            json={"token": token, "stream": "reference", "start_x": 100.0, "end_x": 300.0},
        )
        self.assertEqual(extract_reference.status_code, 200)
        extract_reference_json = extract_reference.get_json()
        self.assertEqual(extract_reference_json["stream"], "reference")
        self.assertIn("sequence_length", extract_reference_json)

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

    def test_extract_sequence_validation_errors(self):
        missing_token_resp = self.client.post(
            "/api/extract_sequence",
            json={"stream": "query", "start_x": 1.0, "end_x": 2.0},
        )
        self.assertEqual(missing_token_resp.status_code, 400)
        self.assertIn("token is required", missing_token_resp.get_json()["error"])

        payload = {
            "input_path": str(ALIGN),
            "reference_annotation_path": str(REF_ANNOT),
            "query_annotation_path": "",
            "params": {"width": "auto", "height": 6.0, "dpi": 100},
        }
        render_resp = self.client.post("/api/render", json=payload)
        token = render_resp.get_json()["token"]

        invalid_stream_resp = self.client.post(
            "/api/extract_sequence",
            json={"token": token, "stream": "invalid", "start_x": 1.0, "end_x": 2.0},
        )
        self.assertEqual(invalid_stream_resp.status_code, 400)
        self.assertIn("stream must be", invalid_stream_resp.get_json()["error"])

        missing_start_resp = self.client.post(
            "/api/extract_sequence",
            json={"token": token, "stream": "query", "end_x": 2.0},
        )
        self.assertEqual(missing_start_resp.status_code, 400)
        self.assertIn("start_x is required", missing_start_resp.get_json()["error"])

        missing_end_resp = self.client.post(
            "/api/extract_sequence",
            json={"token": token, "stream": "query", "start_x": 1.0},
        )
        self.assertEqual(missing_end_resp.status_code, 400)
        self.assertIn("end_x is required", missing_end_resp.get_json()["error"])

    def test_api_render_and_extract_respect_swap_roles(self):
        base_payload = {
            "input_path": str(ALIGN),
            "reference_annotation_path": str(REF_ANNOT),
            "query_annotation_path": "",
            "params": {"width": "auto", "height": 6.0, "dpi": 100},
        }

        default_render = self.client.post("/api/render", json=base_payload)
        self.assertEqual(default_render.status_code, 200)
        default_json = default_render.get_json()
        default_token = default_json["token"]

        swapped_payload = dict(base_payload)
        swapped_payload["swap_roles"] = True
        swapped_render = self.client.post("/api/render", json=swapped_payload)
        self.assertEqual(swapped_render.status_code, 200)
        swapped_json = swapped_render.get_json()
        swapped_token = swapped_json["token"]

        self.assertEqual(swapped_json["query_name"], default_json["reference_name"])
        self.assertEqual(swapped_json["reference_name"], default_json["query_name"])

        default_reference_extract = self.client.post(
            "/api/extract_sequence",
            json={
                "token": default_token,
                "stream": "reference",
                "start_x": 100.0,
                "end_x": 300.0,
            },
        )
        self.assertEqual(default_reference_extract.status_code, 200)
        default_reference_json = default_reference_extract.get_json()

        swapped_query_extract = self.client.post(
            "/api/extract_sequence",
            json={
                "token": swapped_token,
                "stream": "query",
                "start_x": 100.0,
                "end_x": 300.0,
            },
        )
        self.assertEqual(swapped_query_extract.status_code, 200)
        swapped_query_json = swapped_query_extract.get_json()

        self.assertEqual(swapped_query_json["sequence"], default_reference_json["sequence"])
        self.assertEqual(swapped_query_json["stream"], "query")

    def test_preset_recommendation_short_for_small_input(self):
        alignment_path = self._write_temp_alignment(200)
        self.addCleanup(lambda: alignment_path.unlink(missing_ok=True))

        resp = self.client.post(
            "/api/preset_recommendation",
            json={"input_path": str(alignment_path)},
        )
        self.assertEqual(resp.status_code, 200)
        data = resp.get_json()
        self.assertEqual(data["recommended_preset"], "short_lt_5000")
        self.assertEqual(data["basis_length"], 200)
        self.assertEqual(data["threshold"], 5000)

    def test_preset_recommendation_long_for_large_input(self):
        alignment_path = self._write_temp_alignment(6000)
        self.addCleanup(lambda: alignment_path.unlink(missing_ok=True))

        resp = self.client.post(
            "/api/preset_recommendation",
            json={"input_path": str(alignment_path)},
        )
        self.assertEqual(resp.status_code, 200)
        data = resp.get_json()
        self.assertEqual(data["recommended_preset"], "long_ge_5000")
        self.assertEqual(data["basis_length"], 6000)
        self.assertEqual(data["threshold"], 5000)

    def test_preset_recommendation_requires_input_path(self):
        resp = self.client.post("/api/preset_recommendation", json={})
        self.assertEqual(resp.status_code, 400)
        self.assertIn("input_path is required", resp.get_json()["error"])

    def test_preset_recommendation_invalid_file_returns_400(self):
        resp = self.client.post(
            "/api/preset_recommendation",
            json={"input_path": str(ROOT / "tests" / "fixtures" / "does_not_exist.fa")},
        )
        self.assertEqual(resp.status_code, 400)
        self.assertIn("No such file", resp.get_json()["error"])


if __name__ == "__main__":
    unittest.main()
