#!/usr/bin/env python3
"""Loopback-only static app server and GRASP MCDHF job runner."""

from __future__ import annotations

import csv
import io
import json
import mimetypes
import os
import re
import shutil
import subprocess
import tempfile
import threading
import time
import traceback
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.error import URLError
from urllib.parse import parse_qs, unquote, urlencode, urlparse
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parent
GRASP_BIN = Path(os.environ["GRASP_BIN"]).expanduser() if os.environ.get("GRASP_BIN") else None
JOB_ROOT = Path(os.environ.get("ORBITAL_JOB_DIR", Path.home() / ".cache" / "orbital-grasp")).expanduser()
PORT = int(os.environ.get("ORBITAL_PORT", "8765"))
JOB_TIMEOUT_SECONDS = 3600
MAX_REQUEST_BYTES = 1_500_000
MAX_INPUT_CHARS = 500_000
MAX_LOG_CHARS = 30_000
MAX_ENSDF_BYTES = 2_000_000
ENSDF_TIMEOUT_SECONDS = 20
CORE_PROGRAMS = ("rangular", "rwfnestimate", "rmcdhf")
OPTIONAL_PROGRAMS = ("rci",)
JOB_LOCK = threading.Lock()
JOBS = {}
FILENAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
ELEMENT_SYMBOLS = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)


class BackendError(Exception):
    """Expected input or GRASP execution error."""


def fetch_ensdf_levels(z: int, mass: int) -> dict:
    if not 1 <= z <= len(ELEMENT_SYMBOLS):
        raise BackendError("Atomic number must be between 1 and 118.")
    if mass < z or mass > 320:
        raise BackendError("Mass number must be at least Z and no greater than 320.")
    nuclide = f"{mass}{ELEMENT_SYMBOLS[z - 1]}"
    query = urlencode({"fields": "levels", "nuclides": nuclide})
    request = Request(
        f"https://nds.iaea.org/relnsd/v1/data?{query}",
        headers={"Accept": "text/csv", "User-Agent": "OrbitalAtomicStructure/1.0"},
    )
    try:
        with urlopen(request, timeout=ENSDF_TIMEOUT_SECONDS) as response:
            headers = getattr(response, "headers", None)
            charset = headers.get_content_charset() if headers else None
            data = response.read(MAX_ENSDF_BYTES + 1)
    except (OSError, URLError, TimeoutError) as error:
        raise BackendError(f"IAEA LiveChart request failed: {error}") from error
    if len(data) > MAX_ENSDF_BYTES:
        raise BackendError("IAEA LiveChart response exceeds the 2 MB limit.")
    try:
        text = data.decode("utf-8-sig" if charset in (None, "utf-8", "utf-8-sig") else charset)
        reader = csv.DictReader(io.StringIO(text))
        if not reader.fieldnames or not {"z", "n", "symbol", "idx", "energy", "jp"} <= set(reader.fieldnames):
            raise BackendError("IAEA LiveChart returned an unexpected levels dataset.")
        raw_rows = list(reader)
    except (LookupError, UnicodeDecodeError, csv.Error) as error:
        raise BackendError(f"Could not parse the IAEA LiveChart response: {error}") from error

    levels = []
    for row in raw_rows:
        try:
            row_z, neutron_count = int(row["z"]), int(row["n"])
            index = int(row["idx"])
            energy_kev = float(row["energy"] or 0)
        except (TypeError, ValueError):
            continue
        if row_z != z or neutron_count + row_z != mass or index < 0:
            continue
        raw_spin_parity = (row.get("jp") or "").strip()
        match = re.fullmatch(r"(\d+(?:/2)?)([+-])", raw_spin_parity)
        levels.append({
            "id": f"{nuclide}:idx{index}",
            "index": index,
            "k": index + 1,
            "energy_kev": energy_kev,
            "energy_mev": energy_kev / 1000,
            "energy_uncertainty_kev": row.get("unc_e") or None,
            "spin_parity": raw_spin_parity or None,
            "spin": match.group(1) if match else None,
            "parity": match.group(2) if match else None,
            "spin_parity_assignment_exact": match is not None,
            "isospin": row.get("isospin") or None,
            "magnetic_dipole_moment_mu_n": row.get("magnetic_dipole") or None,
            "electric_quadrupole_moment_b": row.get("electric_quadrupole") or None,
            "half_life": row.get("half_life") or None,
            "half_life_unit": row.get("unit_hl") or None,
            "ensdf_publication_cutoff": row.get("ENSDF_publication_cut-off") or None,
            "ensdf_authors": row.get("ENSDF_authors") or None,
        })
    if not levels:
        raise BackendError(f"IAEA LiveChart returned no evaluated levels for {nuclide}.")
    return {
        "source": "IAEA LiveChart / ENSDF",
        "dataset_id": f"IAEA-LiveChart:levels:{nuclide}",
        "nuclide": nuclide,
        "z": z,
        "mass": mass,
        "extraction_date": next((row.get("Extraction_date") for row in raw_rows if row.get("Extraction_date")), None),
        "levels": levels,
        "wavefunctions_included": False,
        "matrix_elements_included": False,
    }


def find_program(name: str) -> Path | None:
    if GRASP_BIN:
        candidate = (GRASP_BIN / name).resolve()
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
        return None
    located = shutil.which(name)
    return Path(located).resolve() if located else None


def engine_status() -> dict:
    programs = {name: find_program(name) is not None for name in CORE_PROGRAMS + OPTIONAL_PROGRAMS}
    missing = [name for name in CORE_PROGRAMS if not programs[name]]
    return {
        "engine": "GRASP",
        "available": not missing,
        "programs": programs,
        "missing": missing,
        "message": "Core GRASP programs are ready." if not missing else "Missing required GRASP programs: " + ", ".join(missing),
    }


def validate_job(payload: object) -> dict:
    if not isinstance(payload, dict):
        raise BackendError("Expected a JSON object.")

    def text_field(name: str, *, required: bool = True, max_chars: int = MAX_INPUT_CHARS) -> str:
        value = payload.get(name)
        if value is None and not required:
            return ""
        if not isinstance(value, str):
            raise BackendError(f"{name} must be text.")
        if len(value) > max_chars:
            raise BackendError(f"{name} exceeds the {max_chars:,}-character limit.")
        if "\0" in value:
            raise BackendError(f"{name} contains a null character.")
        if required and not value.strip():
            raise BackendError(f"{name} is required.")
        return value

    try:
        z = int(payload.get("z"))
        charge = int(payload.get("charge"))
        mass = int(payload.get("mass"))
    except (TypeError, ValueError):
        raise BackendError("Atomic number, charge, and mass number must be integers.") from None
    if not 1 <= z <= 118:
        raise BackendError("Atomic number must be between 1 and 118.")
    if not -118 <= charge < z:
        raise BackendError("Charge must be between -118 and Z−1 so the selected species has at least one electron.")
    if mass < z or mass > 320:
        raise BackendError("Mass number must be at least Z and no greater than 320.")
    if z - charge < 1:
        raise BackendError("A zero-electron species has no electronic wavefunction to calculate.")
    if payload.get("confirmed") is not True:
        raise BackendError("Confirm that the GRASP inputs match the selected species and state.")

    job = {
        "z": z,
        "charge": charge,
        "mass": mass,
        "configuration": text_field("configuration", max_chars=500_000),
        "isodata": text_field("isodata", max_chars=64_000),
        "rangular_answers": text_field("rangular_answers", max_chars=64_000),
        "rwfnestimate_answers": text_field("rwfnestimate_answers", max_chars=64_000),
        "rmcdhf_answers": text_field("rmcdhf_answers", max_chars=64_000),
        "rci_answers": text_field("rci_answers", required=False, max_chars=64_000),
    }
    return job


def _set_job(job_id: str, **values) -> None:
    with JOB_LOCK:
        JOBS[job_id].update(values)


def _run_program(job_id: str, work_dir: Path, program: str, answers: str, percent: int) -> None:
    executable = find_program(program)
    if executable is None:
        raise BackendError(f"GRASP program '{program}' is not executable. Check GRASP_BIN and the installation.")
    if answers and not answers.endswith("\n"):
        answers += "\n"
    log_path = work_dir / "run.log"
    with log_path.open("a", encoding="utf-8", errors="replace") as log_file:
        log_file.write(f"\n===== {program} =====\n")
        log_file.flush()
        try:
            result = subprocess.run(
                [str(executable)],
                input=answers,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                cwd=work_dir,
                text=True,
                timeout=JOB_TIMEOUT_SECONDS,
                check=False,
                shell=False,
            )
        except subprocess.TimeoutExpired as error:
            raise BackendError(f"{program} exceeded the {JOB_TIMEOUT_SECONDS // 60}-minute job limit.") from error
    log_tail = log_path.read_text(encoding="utf-8", errors="replace")[-MAX_LOG_CHARS:]
    if result.returncode != 0:
        raise BackendError(f"{program} exited with status {result.returncode}.\n{log_tail}")
    _set_job(job_id, progress=percent, step=f"Completed {program}", log=log_tail)


def _output_files(work_dir: Path) -> list[dict]:
    files = []
    for path in sorted(work_dir.iterdir()):
        if not path.is_file() or path.is_symlink() or not FILENAME_RE.fullmatch(path.name):
            continue
        if path.name in {"rcsf.inp", "isodata"}:
            continue
        if path.suffix.lower() not in {".w", ".c", ".out", ".log", ".rw"}:
            continue
        files.append({"name": path.name, "size": path.stat().st_size})
    return files


def run_job(job_id: str, payload: dict, work_dir: Path) -> None:
    (work_dir / "rcsf.inp").write_text(payload["configuration"], encoding="utf-8")
    (work_dir / "isodata").write_text(payload["isodata"], encoding="utf-8")
    stages = [
        ("rangular", "rangular_answers", 18, "Build angular data"),
        ("rwfnestimate", "rwfnestimate_answers", 36, "Estimate starting orbitals"),
        ("rmcdhf", "rmcdhf_answers", 84, "Run MCDHF self-consistent field"),
    ]
    if payload["rci_answers"].strip():
        stages.append(("rci", "rci_answers", 96, "Run requested RCI calculation"))
    _set_job(job_id, status="running", progress=5, step="Starting local GRASP workflow")
    for program, answer_key, percent, description in stages:
        _set_job(job_id, step=description)
        _run_program(job_id, work_dir, program, payload[answer_key], percent)
    files = _output_files(work_dir)
    log = (work_dir / "run.log").read_text(encoding="utf-8", errors="replace")[-MAX_LOG_CHARS:]
    _set_job(
        job_id,
        status="completed",
        progress=100,
        step="GRASP workflow completed",
        log=log,
        files=files,
        completed_at=time.time(),
    )


def _worker(job_id: str, payload: dict, work_dir: Path) -> None:
    try:
        run_job(job_id, payload, work_dir)
    except (BackendError, OSError) as error:
        _set_job(
            job_id,
            status="failed",
            step="GRASP workflow failed",
            error=str(error),
            log=((work_dir / "run.log").read_text(encoding="utf-8", errors="replace")[-MAX_LOG_CHARS:]
                 if (work_dir / "run.log").exists() else ""),
            completed_at=time.time(),
        )
    except Exception as error:
        _set_job(
            job_id,
            status="failed",
            step="Unexpected backend error",
            error=f"{type(error).__name__}: {error}",
            log=traceback.format_exc()[-MAX_LOG_CHARS:],
            completed_at=time.time(),
        )


class Handler(BaseHTTPRequestHandler):
    server_version = "OrbitalLocal/1.0"

    def log_message(self, format_string, *args):
        print(f"[orbital] {self.address_string()} {format_string % args}")

    def _json(self, status: int, value: dict) -> None:
        body = json.dumps(value).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _file(self, path: Path, download_name: str | None = None) -> None:
        body = path.read_bytes()
        content_type = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("X-Content-Type-Options", "nosniff")
        if download_name:
            self.send_header("Content-Disposition", f'attachment; filename="{download_name}"')
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        request = urlparse(self.path)
        parts = [unquote(part) for part in request.path.split("/") if part]
        if parts == ["api", "ensdf", "levels"]:
            query = parse_qs(request.query)
            try:
                z = int(query.get("z", [""])[0])
                mass = int(query.get("a", [""])[0])
                levels = fetch_ensdf_levels(z, mass)
            except (TypeError, ValueError):
                self._json(400, {"error": "Specify integer z and a query parameters."})
                return
            except BackendError as error:
                status = 400 if "must be" in str(error) or "Atomic number" in str(error) else 502
                self._json(status, {"error": str(error)})
                return
            self._json(200, levels)
            return
        if parts == ["api", "health"]:
            self._json(200, engine_status())
            return
        if len(parts) == 3 and parts[:2] == ["api", "jobs"]:
            with JOB_LOCK:
                job = JOBS.get(parts[2])
                public_job = dict(job) if job else None
            if not public_job:
                self._json(404, {"error": "Job not found."})
                return
            public_job.pop("work_dir", None)
            self._json(200, public_job)
            return
        if len(parts) == 5 and parts[:2] == ["api", "jobs"] and parts[3] == "files":
            # This route is retained for compatibility if the path has a trailing segment.
            self._file_for_job(parts[2], parts[4])
            return
        if len(parts) == 4 and parts[:2] == ["api", "jobs"] and parts[3] == "files":
            self._json(400, {"error": "Specify an output filename."})
            return
        if len(parts) == 5 and parts[:2] == ["api", "jobs"] and parts[3] == "file":
            self._file_for_job(parts[2], parts[4])
            return

        relative = request.path.lstrip("/") or "index.html"
        target = (ROOT / relative).resolve()
        if not target.is_relative_to(ROOT) or not target.is_file():
            self._json(404, {"error": "Not found."})
            return
        if target.name.startswith(".") or target.suffix == ".py":
            self._json(404, {"error": "Not found."})
            return
        self._file(target)

    def _file_for_job(self, job_id: str, filename: str) -> None:
        if not FILENAME_RE.fullmatch(filename):
            self._json(400, {"error": "Invalid output filename."})
            return
        with JOB_LOCK:
            job = JOBS.get(job_id)
            work_dir = Path(job["work_dir"]) if job and job.get("work_dir") else None
            allowed = {item["name"] for item in job.get("files", [])} if job else set()
        if work_dir is None or filename not in allowed:
            self._json(404, {"error": "Output file not found."})
            return
        candidate = work_dir / filename
        if candidate.is_symlink():
            self._json(404, {"error": "Output file not found."})
            return
        path = candidate.resolve()
        if path.parent != work_dir.resolve() or not path.is_file():
            self._json(404, {"error": "Output file not found."})
            return
        self._file(path, filename)

    def do_POST(self):
        if urlparse(self.path).path != "/api/jobs":
            self._json(404, {"error": "Not found."})
            return
        if self.headers.get("Content-Type", "").split(";", 1)[0].strip().lower() != "application/json":
            self._json(415, {"error": "Send application/json."})
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            self._json(400, {"error": "Invalid Content-Length."})
            return
        if length < 1 or length > MAX_REQUEST_BYTES:
            self._json(413, {"error": f"Request body must be between 1 and {MAX_REQUEST_BYTES:,} bytes."})
            return
        try:
            payload = validate_job(json.loads(self.rfile.read(length)))
        except (json.JSONDecodeError, UnicodeDecodeError):
            self._json(400, {"error": "Invalid JSON request body."})
            return
        except BackendError as error:
            self._json(400, {"error": str(error)})
            return

        engine = engine_status()
        if not engine["available"]:
            self._json(503, {"error": engine["message"], "engine": engine})
            return
        with JOB_LOCK:
            if any(job["status"] in {"queued", "running"} for job in JOBS.values()):
                self._json(429, {"error": "A GRASP job is already running. Wait for it to finish before starting another."})
                return
            job_id = uuid.uuid4().hex
            JOB_ROOT.mkdir(parents=True, exist_ok=True)
            work_dir = Path(tempfile.mkdtemp(prefix=f"{job_id}-", dir=JOB_ROOT)).resolve()
            JOBS[job_id] = {
                "id": job_id,
                "status": "queued",
                "step": "Queued",
                "progress": 0,
                "log": "",
                "error": None,
                "files": [],
                "work_dir": str(work_dir),
                "species": {"z": payload["z"], "charge": payload["charge"], "mass": payload["mass"]},
                "uses_rci": bool(payload["rci_answers"].strip()),
                "created_at": time.time(),
            }
        thread = threading.Thread(target=_worker, args=(job_id, payload, work_dir), daemon=True)
        thread.start()
        self._json(202, {"id": job_id, "status": "queued"})


def main() -> None:
    server = ThreadingHTTPServer(("127.0.0.1", PORT), Handler)
    print(f"Orbital local server: http://127.0.0.1:{PORT}")
    print("Bound to 127.0.0.1 only. GRASP status:", engine_status()["message"])
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping Orbital local server.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
