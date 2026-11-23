from io import BytesIO
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.parse import urlparse, parse_qs

import server


def valid_payload():
    return {
        "z": 2,
        "charge": 0,
        "mass": 4,
        "confirmed": True,
        "configuration": "CSF input\n",
        "isodata": "isotope data\n",
        "rangular_answers": "y\n",
        "rwfnestimate_answers": "y\n",
        "rmcdhf_answers": "y\n",
        "rci_answers": "",
    }


class BackendTests(unittest.TestCase):
    def test_fetch_ensdf_levels_converts_evaluated_metadata_without_inventing_wavefunctions(self):
        body = (
            "z,n,symbol,idx,energy,unc_e,jp,isospin,magnetic_dipole,electric_quadrupole,"
            "half_life,unit_hl,ENSDF_publication_cut-off,ENSDF_authors,Extraction_date\n"
            "6,7,C,0,0,,1/2-,0.5,0.702,,,STABLE,1-Aug-2024,Example Author,2026-11-23\n"
            "6,7,C,1,3089.451,0.019,(1/2+),,,,,fs,1-Aug-2024,Example Author,2026-11-23\n"
            "6,8,C,0,0,,0+,,,,,STABLE,1-Jan-2020,Other Author,2026-11-23\n"
        ).encode()
        with patch.object(server, "urlopen", return_value=BytesIO(body)) as open_url:
            result = server.fetch_ensdf_levels(6, 13)

        request = open_url.call_args.args[0]
        query = parse_qs(urlparse(request.full_url).query)
        self.assertEqual(query["nuclides"], ["13C"])
        self.assertEqual(len(result["levels"]), 2)
        self.assertAlmostEqual(result["levels"][1]["energy_mev"], 3.089451)
        self.assertEqual(result["levels"][1]["spin_parity"], "(1/2+)")
        self.assertFalse(result["levels"][1]["spin_parity_assignment_exact"])
        self.assertFalse(result["wavefunctions_included"])
        self.assertFalse(result["matrix_elements_included"])

    def test_fetch_ensdf_levels_rejects_invalid_species_and_remote_failures(self):
        with self.assertRaisesRegex(server.BackendError, "Mass number"):
            server.fetch_ensdf_levels(6, 5)
        with patch.object(server, "urlopen", side_effect=TimeoutError("timeout")):
            with self.assertRaisesRegex(server.BackendError, "IAEA LiveChart request failed"):
                server.fetch_ensdf_levels(6, 13)

    def test_validate_rejects_unconfirmed_or_inconsistent_species(self):
        payload = valid_payload()
        payload["confirmed"] = False
        with self.assertRaisesRegex(server.BackendError, "Confirm"):
            server.validate_job(payload)

        payload = valid_payload()
        payload["mass"] = 1
        with self.assertRaisesRegex(server.BackendError, "Mass number"):
            server.validate_job(payload)

    def test_validate_species_boundaries(self):
        payload = valid_payload()
        payload.update(z=118, charge=117, mass=294)
        result = server.validate_job(payload)
        self.assertEqual(result["z"], 118)
        self.assertEqual(result["charge"], 117)
        self.assertEqual(result["mass"], 294)

        payload = valid_payload()
        payload.update(z=118, charge=-119)
        with self.assertRaisesRegex(server.BackendError, "Charge must"):
            server.validate_job(payload)

    def test_runs_fixed_grasp_programs_and_returns_native_output(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            binary_dir = root / "bin"
            work_dir = root / "job"
            binary_dir.mkdir()
            work_dir.mkdir()
            original_bin = server.GRASP_BIN
            server.GRASP_BIN = binary_dir
            try:
                for name in server.CORE_PROGRAMS:
                    program = binary_dir / name
                    body = [
                        "#!/usr/bin/env python3",
                        "from pathlib import Path",
                        "import sys",
                        "sys.stdin.read()",
                    ]
                    if name == "rmcdhf":
                        body.append("Path('rwfn.out').write_bytes(b'GRASP orbital data')")
                    program.write_text("\n".join(body) + "\n", encoding="utf-8")
                    program.chmod(0o755)
                job_id = "test-job"
                with server.JOB_LOCK:
                    server.JOBS[job_id] = {}
                server.run_job(job_id, server.validate_job(valid_payload()), work_dir)
                with server.JOB_LOCK:
                    result = dict(server.JOBS[job_id])
                self.assertEqual(result["status"], "completed")
                self.assertEqual(result["progress"], 100)
                self.assertIn({"name": "rwfn.out", "size": len(b"GRASP orbital data")}, result["files"])
                self.assertIn("===== rmcdhf =====", result["log"])
            finally:
                server.GRASP_BIN = original_bin
                with server.JOB_LOCK:
                    server.JOBS.pop("test-job", None)

    def test_optional_rci_runs_only_when_transcript_is_supplied(self):
        payload = valid_payload()
        payload["rci_answers"] = "run rci\n"
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            binary_dir = root / "bin"
            work_dir = root / "job"
            binary_dir.mkdir()
            work_dir.mkdir()
            original_bin = server.GRASP_BIN
            server.GRASP_BIN = binary_dir
            try:
                for name in server.CORE_PROGRAMS + ("rci",):
                    program = binary_dir / name
                    program.write_text("#!/bin/sh\ncat >/dev/null\n", encoding="utf-8")
                    program.chmod(0o755)
                job_id = "test-rci"
                with server.JOB_LOCK:
                    server.JOBS[job_id] = {}
                server.run_job(job_id, server.validate_job(payload), work_dir)
                with server.JOB_LOCK:
                    result = dict(server.JOBS[job_id])
                self.assertIn("===== rci =====", result["log"])
            finally:
                server.GRASP_BIN = original_bin
                with server.JOB_LOCK:
                    server.JOBS.pop("test-rci", None)


if __name__ == "__main__":
    unittest.main()
