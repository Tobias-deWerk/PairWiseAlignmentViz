from __future__ import annotations

import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from core import genome_workflow as gw


class GenomeWorkflowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)

        self._old_root = gw.WORK_ROOT
        self._old_uploads = gw.UPLOADS
        self._old_jobs = gw.JOBS

        gw.WORK_ROOT = Path(self.tempdir.name) / "work"
        gw.UPLOADS = {}
        gw.JOBS = {}

        self.addCleanup(self._restore_globals)

    def _restore_globals(self) -> None:
        gw.WORK_ROOT = self._old_root
        gw.UPLOADS = self._old_uploads
        gw.JOBS = self._old_jobs

    def test_create_upload_from_streams_accepts_fasta(self):
        payload = gw.create_upload_from_streams(
            "query.fa",
            io.BytesIO(b">query\nACGTACGT\n"),
            "reference.fa",
            io.BytesIO(b">reference\nACGTACGT\n"),
        )
        self.assertIn("upload_id", payload)
        self.assertEqual(payload["query_length"], 8)
        self.assertEqual(payload["reference_length"], 8)

    def test_build_selected_sequences_reorients_reverse_blocks(self):
        upload = gw.UploadRecord(
            upload_id="u1",
            dir_path=Path(self.tempdir.name),
            query_path=Path(self.tempdir.name) / "q.fa",
            reference_path=Path(self.tempdir.name) / "r.fa",
            query_name="q",
            reference_name="r",
            query_sequence="AAAACCCCGGGGTTTT",
            reference_sequence="AAAACCCCGGGGTTTT",
            query_length=16,
            reference_length=16,
            dotplot={
                "blocks": [
                    {
                        "block_id": "b0",
                        "q_start": 1,
                        "q_end": 4,
                        "r_start": 1,
                        "r_end": 4,
                        "orientation": "forward",
                        "score": 99.0,
                    },
                    {
                        "block_id": "b1",
                        "q_start": 5,
                        "q_end": 8,
                        "r_start": 5,
                        "r_end": 8,
                        "orientation": "reverse",
                        "score": 98.0,
                    },
                ]
            },
        )

        (
            query_concat,
            reference_concat,
            inversion_ranges,
            duplication_ranges,
            selected,
            inter_block_count,
            stitching_notes,
        ) = gw._build_selected_sequences(
            upload,
            [
                {"block_id": "b0", "include": True},
                {"block_id": "b1", "include": True},
            ],
        )

        self.assertEqual(query_concat, "AAAAGGGG")
        self.assertEqual(reference_concat, "AAAACCCC")
        self.assertEqual(inversion_ranges, [(5, 8)])
        self.assertEqual(duplication_ranges, [])
        self.assertEqual(len(selected), 2)
        self.assertEqual(inter_block_count, 0)
        self.assertEqual(stitching_notes, [])

    def test_build_selected_sequences_includes_intervals_when_enabled(self):
        upload = gw.UploadRecord(
            upload_id="u1b",
            dir_path=Path(self.tempdir.name),
            query_path=Path(self.tempdir.name) / "q.fa",
            reference_path=Path(self.tempdir.name) / "r.fa",
            query_name="q",
            reference_name="r",
            query_sequence="AAAACCCCGGGGTTTT",
            reference_sequence="AAAACCCCGGGGTTTT",
            query_length=16,
            reference_length=16,
            dotplot={
                "blocks": [
                    {
                        "block_id": "b0",
                        "q_start": 1,
                        "q_end": 4,
                        "r_start": 1,
                        "r_end": 4,
                        "orientation": "forward",
                        "score": 99.0,
                    },
                    {
                        "block_id": "b1",
                        "q_start": 9,
                        "q_end": 12,
                        "r_start": 9,
                        "r_end": 12,
                        "orientation": "forward",
                        "score": 99.0,
                    },
                ]
            },
        )

        (
            query_concat,
            reference_concat,
            _inversion_ranges,
            _duplication_ranges,
            _selected,
            inter_block_count,
            stitching_notes,
        ) = gw._build_selected_sequences(
            upload,
            [
                {"block_id": "b0", "include": True},
                {"block_id": "b1", "include": True},
            ],
            include_inter_block_intervals=True,
        )

        self.assertEqual(query_concat, "AAAACCCCGGGG")
        self.assertEqual(reference_concat, "AAAACCCCGGGG")
        self.assertEqual(inter_block_count, 1)
        self.assertEqual(stitching_notes, [])

    def test_build_selected_sequences_marks_duplicate_overlap(self):
        upload = gw.UploadRecord(
            upload_id="u1c",
            dir_path=Path(self.tempdir.name),
            query_path=Path(self.tempdir.name) / "q.fa",
            reference_path=Path(self.tempdir.name) / "r.fa",
            query_name="q",
            reference_name="r",
            query_sequence="AAAACCCCGGGGTTTT",
            reference_sequence="AAAACCCCGGGGTTTT",
            query_length=16,
            reference_length=16,
            dotplot={
                "blocks": [
                    {
                        "block_id": "b0",
                        "q_start": 1,
                        "q_end": 8,
                        "r_start": 1,
                        "r_end": 8,
                        "orientation": "forward",
                        "score": 99.0,
                    },
                    {
                        "block_id": "b1",
                        "q_start": 5,
                        "q_end": 12,
                        "r_start": 9,
                        "r_end": 16,
                        "orientation": "forward",
                        "score": 98.0,
                    },
                ]
            },
        )
        (
            _query_concat,
            _reference_concat,
            _inversion_ranges,
            duplication_ranges,
            selected,
            _inter_block_count,
            _stitching_notes,
        ) = gw._build_selected_sequences(
            upload,
            [
                {"block_id": "b0", "include": True},
                {"block_id": "b1", "include": True},
            ],
        )
        self.assertEqual(len(duplication_ranges), 2)
        self.assertTrue(all(item["is_duplicate"] for item in selected))
        self.assertIn("Duplicated on query 1-8; reference 1-8", duplication_ranges[0]["pointer_text"])

    def test_map_local_ranges_to_columns_handles_gaps(self):
        cols = gw._map_local_ranges_to_columns("AA--CCGG", [(3, 6)])
        self.assertEqual(cols, [(5, 8)])

    def test_start_alignment_job_requires_selected_blocks(self):
        upload = gw.UploadRecord(
            upload_id="u2",
            dir_path=Path(self.tempdir.name),
            query_path=Path(self.tempdir.name) / "q.fa",
            reference_path=Path(self.tempdir.name) / "r.fa",
            query_name="q",
            reference_name="r",
            query_sequence="AAAA",
            reference_sequence="AAAA",
            query_length=4,
            reference_length=4,
            dotplot={"blocks": [{"block_id": "b0", "q_start": 1, "q_end": 4, "r_start": 1, "r_end": 4, "orientation": "forward", "score": 100.0}]},
        )
        gw.UPLOADS[upload.upload_id] = upload

        with self.assertRaises(ValueError):
            gw.start_alignment_job(upload.upload_id, [{"block_id": "b0", "include": False}])

    def test_start_alignment_job_persists_align_options(self):
        upload = gw.UploadRecord(
            upload_id="u2b",
            dir_path=Path(self.tempdir.name),
            query_path=Path(self.tempdir.name) / "q.fa",
            reference_path=Path(self.tempdir.name) / "r.fa",
            query_name="q",
            reference_name="r",
            query_sequence="AAAA",
            reference_sequence="AAAA",
            query_length=4,
            reference_length=4,
            dotplot={"blocks": [{"block_id": "b0", "q_start": 1, "q_end": 4, "r_start": 1, "r_end": 4, "orientation": "forward", "score": 100.0}]},
        )
        gw.UPLOADS[upload.upload_id] = upload

        old_run = gw._run_alignment_job
        try:
            gw._run_alignment_job = lambda _job, _selected, _opts: None
            job_id = gw.start_alignment_job(
                upload.upload_id,
                [{"block_id": "b0", "include": True}],
                align_options={"include_inter_block_intervals": True},
            )
            self.assertIn(job_id, gw.JOBS)
            self.assertTrue(gw.JOBS[job_id].result["align_options_used"]["include_inter_block_intervals"])
        finally:
            gw._run_alignment_job = old_run

    def test_parse_dotplot_options_defaults(self):
        parsed = gw._parse_dotplot_options({})
        self.assertEqual(parsed["match_mode"], "maxmatch")
        self.assertNotIn("mincluster", parsed)
        self.assertNotIn("diagfactor", parsed)
        self.assertNotIn("breaklen", parsed)

    def test_parse_align_options_defaults_and_validation(self):
        parsed = gw._parse_align_options({})
        self.assertIn("include_inter_block_intervals", parsed)
        self.assertFalse(parsed["include_inter_block_intervals"])
        parsed_true = gw._parse_align_options({"include_inter_block_intervals": "true"})
        self.assertTrue(parsed_true["include_inter_block_intervals"])
        with self.assertRaises(ValueError):
            gw._parse_align_options({"include_inter_block_intervals": "invalid"})

    def test_parse_dotplot_options_valid_values(self):
        parsed = gw._parse_dotplot_options(
            {
                "match_mode": "mumreference",
                "mincluster": "42",
                "diagfactor": "0.3",
                "breaklen": 150,
            }
        )
        self.assertEqual(parsed["match_mode"], "mumreference")
        self.assertEqual(parsed["mincluster"], 42)
        self.assertAlmostEqual(parsed["diagfactor"], 0.3)
        self.assertEqual(parsed["breaklen"], 150)

    def test_parse_dotplot_options_invalid_values_raise(self):
        with self.assertRaises(ValueError):
            gw._parse_dotplot_options({"match_mode": "invalid"})
        with self.assertRaises(ValueError):
            gw._parse_dotplot_options({"mincluster": 0})
        with self.assertRaises(ValueError):
            gw._parse_dotplot_options({"diagfactor": -0.1})
        with self.assertRaises(ValueError):
            gw._parse_dotplot_options({"breaklen": "abc"})

    def test_build_nucmer_command_uses_single_mode_and_optional_flags(self):
        cmd = gw._build_nucmer_command(
            "/usr/bin/nucmer",
            options={"match_mode": "mum", "mincluster": 20, "diagfactor": 0.4, "breaklen": 99},
            prefix=Path("/tmp/prefix"),
            reference_path=Path("/tmp/ref.fa"),
            query_path=Path("/tmp/q.fa"),
        )
        self.assertIn("--mum", cmd)
        self.assertNotIn("--maxmatch", cmd)
        self.assertNotIn("--mumreference", cmd)
        self.assertIn("-c", cmd)
        self.assertIn("-D", cmd)
        self.assertIn("-b", cmd)

    def test_start_dotplot_job_persists_options_in_result(self):
        upload = gw.UploadRecord(
            upload_id="u3",
            dir_path=Path(self.tempdir.name),
            query_path=Path(self.tempdir.name) / "q.fa",
            reference_path=Path(self.tempdir.name) / "r.fa",
            query_name="q",
            reference_name="r",
            query_sequence="AAAA",
            reference_sequence="AAAA",
            query_length=4,
            reference_length=4,
            dotplot={},
        )
        gw.UPLOADS[upload.upload_id] = upload

        old_run = gw._run_dotplot_job
        try:
            gw._run_dotplot_job = lambda _job: None
            job_id = gw.start_dotplot_job(
                upload.upload_id,
                nucmer_options={"match_mode": "mumreference", "mincluster": 30},
            )
            self.assertIn(job_id, gw.JOBS)
            self.assertEqual(gw.JOBS[job_id].result["nucmer_options"]["match_mode"], "mumreference")
            self.assertEqual(gw.JOBS[job_id].result["nucmer_options"]["mincluster"], 30)
        finally:
            gw._run_dotplot_job = old_run

    def test_parse_show_coords_handles_separator_style_output(self):
        text = """
        /tmp/ref.fa /tmp/query.fa
        NUCMER
        [S1] [E1] | [S2] [E2] | [LEN 1] [LEN 2] | [% IDY]
        ===================================================
        2128 2198 | 74501 74432 | 71 70 | 95.83
        17507 34174 | 1292 18004 | 16668 16713 | 96.86
        """
        parsed = gw._parse_show_coords(text)
        self.assertGreaterEqual(parsed["parsed_rows_count"], 2)
        self.assertGreater(len(parsed["blocks"]), 0)
        self.assertGreater(len(parsed["points"]), 0)

    def test_parse_show_coords_handles_tabular_output(self):
        text = (
            "2128\t2198\t74501\t74432\t71\t70\t95.83\t142172\t74815\t0.05\t0.09\tref\tquery\n"
            "17507\t34174\t1292\t18004\t16668\t16713\t96.86\t142172\t74815\t11.72\t22.34\tref\tquery\n"
        )
        parsed = gw._parse_show_coords(text)
        self.assertEqual(parsed["parsed_rows_count"], 2)
        self.assertEqual(parsed["skipped_rows_count"], 0)
        self.assertEqual(len(parsed["blocks"]), 2)
        self.assertEqual(parsed["blocks"][0]["orientation"], "reverse")
        self.assertEqual(parsed["blocks"][1]["orientation"], "forward")

    def test_parse_show_coords_skips_malformed_lines_but_keeps_valid(self):
        text = (
            "bad line that should skip\n"
            "2128\t2198\t74501\t74432\t71\t70\t95.83\n"
            "too\tshort\trow\n"
            "17507\t34174\t1292\t18004\t16668\t16713\t96.86\n"
        )
        parsed = gw._parse_show_coords(text)
        self.assertEqual(parsed["parsed_rows_count"], 2)
        self.assertGreaterEqual(parsed["skipped_rows_count"], 1)
        self.assertEqual(len(parsed["blocks"]), 2)

    def test_run_dotplot_job_writes_metadata_counts(self):
        upload_dir = Path(self.tempdir.name) / "dotplot_job"
        upload_dir.mkdir(parents=True, exist_ok=True)
        upload = gw.UploadRecord(
            upload_id="u4",
            dir_path=upload_dir,
            query_path=upload_dir / "query.fa",
            reference_path=upload_dir / "reference.fa",
            query_name="q",
            reference_name="r",
            query_sequence="AAAA",
            reference_sequence="AAAA",
            query_length=4,
            reference_length=4,
            dotplot={},
        )
        upload.query_path.write_text(">q\nAAAA\n", encoding="utf-8")
        upload.reference_path.write_text(">r\nAAAA\n", encoding="utf-8")
        gw.UPLOADS[upload.upload_id] = upload

        job = gw.JobRecord(job_id="j1", job_type="dotplot", upload_id=upload.upload_id, result={"nucmer_options": {"match_mode": "maxmatch"}})

        coords_text = "2128\t2198\t74501\t74432\t71\t70\t95.83\t142172\t74815\t0.05\t0.09\tref\tquery\n"

        def fake_require_command(name: str) -> str:
            return f"/usr/bin/{name}"

        class DummyProcess:
            def __init__(self, stdout: str = "") -> None:
                self.stdout = stdout
                self.stderr = ""

        def fake_run_command(args, *, cwd, timeout_seconds=7200):
            if "show-coords" in str(args[0]):
                return DummyProcess(stdout=coords_text)
            (upload_dir / "dotplot.delta").write_text("delta", encoding="utf-8")
            return DummyProcess(stdout="")

        with patch.object(gw, "_require_command", side_effect=fake_require_command), patch.object(
            gw, "_run_command", side_effect=fake_run_command
        ):
            gw._run_dotplot_job(job)

        self.assertEqual(job.status, "done")
        self.assertIn("nucmer_command", job.result)
        self.assertIn("coords_command", job.result)
        self.assertIn("delta_path", job.result)
        self.assertGreaterEqual(job.result.get("parsed_rows_count", 0), 1)
        self.assertIn("skipped_rows_count", job.result)

    def test_map_labeled_local_ranges_to_columns_keeps_payload(self):
        mapped = gw._map_labeled_local_ranges_to_columns(
            "AA--CCGG",
            [{"start_local": 3, "end_local": 6, "label": "DUP", "pointer_text": "Duplicated on query 3-6; reference 3-6"}],
        )
        self.assertEqual(len(mapped), 1)
        self.assertEqual(mapped[0]["start_column"], 5)
        self.assertEqual(mapped[0]["end_column"], 8)
        self.assertEqual(mapped[0]["label"], "DUP")


if __name__ == "__main__":
    unittest.main()
