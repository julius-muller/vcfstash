from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict

import requests

import pytest

from vcfcache.utils.paths import get_vcfcache_root


def _run(cmd: list[str], env: dict[str, str], cwd: Path | None = None) -> str:
    res = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert res.returncode == 0, (
        f"Command failed: {' '.join(cmd)}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
    )
    return res.stdout


def _load_dotenv_if_present() -> Dict[str, str]:
    """Load simple KEY=VALUE pairs from repo .env if present.

    This avoids adding a dependency on python-dotenv.
    Only sets keys that aren't already in os.environ.
    """
    env_path = get_vcfcache_root() / ".env"
    out: Dict[str, str] = {}
    if not env_path.exists():
        return out
    for line in env_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, val = line.split("=", 1)
        key = key.strip()
        val = val.strip().strip("'").strip('"')
        if key and key not in os.environ:
            out[key] = val
    return out


def _delete_sandbox_deposition(dep_id: str, token: str) -> None:
    url = f"https://sandbox.zenodo.org/api/deposit/depositions/{dep_id}"
    try:
        resp = requests.delete(url, params={"access_token": token}, timeout=30)
        # Zenodo returns 204 on success; ignore failures since published records
        # may not be deletable.
        if resp.status_code not in (200, 202, 204, 404, 409):
            resp.raise_for_status()
    except Exception:
        # Best-effort teardown only.
        return


@pytest.mark.integration
@pytest.mark.skipif(
    os.environ.get("VCFCACHE_RUN_ZENODO_SANDBOX", "0") != "1",
    reason="Requires Zenodo sandbox network + VCFCACHE_RUN_ZENODO_SANDBOX=1",
)
def test_push_and_pull_roundtrip_via_sandbox(tmp_path: Path):
    os.environ.update(_load_dotenv_if_present())
    sandbox = os.environ.get("ZENODO_SANDBOX", "0") == "1"
    token = (
        os.environ.get("ZENODO_SANDBOX_TOKEN")
        if sandbox
        else os.environ.get("ZENODO_TOKEN")
    ) or os.environ.get("ZENODO_TOKEN")
    if not token:
        pytest.skip("Zenodo token not set; skipping sandbox roundtrip")

    vcfcache_root = get_vcfcache_root()
    input_bcf = vcfcache_root / "tests" / "data" / "nodata" / "gnomad_test.bcf"
    anno_config = vcfcache_root / "tests" / "config" / "test_annotation.yaml"
    params_template = vcfcache_root / "tests" / "config" / "test_params.yaml"

    work_dir = tmp_path / "work"
    work_dir.mkdir()
    cache_root = tmp_path / "db"
    pulled_dir = tmp_path / "pulled"
    home_dir = tmp_path / "home"
    home_dir.mkdir()

    params_file = tmp_path / "params.yaml"
    params_file.write_text(
        params_template.read_text().replace("${VCFCACHE_ROOT}", str(vcfcache_root))
    )

    env = os.environ.copy()
    env.update(
        {
            "ZENODO_SANDBOX": "1",
            "ZENODO_TOKEN": token,
            "ZENODO_SANDBOX_TOKEN": token,
            # Ensure downloads/extraction happen in this test's isolated HOME,
            # even if the repo .env sets VCFCACHE_DIR globally (e.g. VCFCACHE_DIR=/tmp).
            "VCFCACHE_DIR": str(home_dir / ".cache" / "vcfcache"),
            "HOME": str(home_dir),
        }
    )

    # Build a tiny cache locally (mock annotation)
    _run(
        [
            "vcfcache",
            "blueprint-init",
            "--vcf",
            str(input_bcf),
            "--output",
            str(cache_root),
            "-y",
            str(params_file),
            "--force",
        ],
        env=env,
        cwd=work_dir,
    )
    _run(
        [
            "vcfcache",
            "cache-build",
            "--name",
            "test_anno",
            "--db",
            str(cache_root),
            "-a",
            str(anno_config),
            "-y",
            str(params_file),
            "--force",
        ],
        env=env,
        cwd=work_dir,
    )

    dep_id = None
    doi = None
    try:
        # Push to sandbox Zenodo (upload each run)
        push_out = _run(
            [
                "vcfcache",
                "push",
                "--cache-dir",
                str(cache_root),
                "--test",
                "--publish",
            ],
            env=env,
            cwd=work_dir,
        )

        m = re.search(r"Deposition ID:\s*(\S+)", push_out)
        if m:
            dep_id = m.group(1)
        m = re.search(r"DOI:\s*(\S+)", push_out)
        assert m, f"Could not parse DOI from push output:\n{push_out}"
        doi = m.group(1)

        # Pull back from sandbox and validate structure
        # cache-build --doi downloads to ~/.cache/vcfcache/caches/ or blueprints/
        # --debug flag uses sandbox mode for Zenodo operations
        _run(
            [
                "vcfcache",
                "cache-build",
                "--debug",
                "--doi",
                doi,
            ],
            env=env,
            cwd=work_dir,
        )

        # Check if downloaded to caches or blueprints directory
        caches_dir = home_dir / ".cache" / "vcfcache" / "caches"
        blueprints_dir = home_dir / ".cache" / "vcfcache" / "blueprints"

        # Should be in caches (since we pushed a cache)
        if caches_dir.exists():
            extracted_roots = list(caches_dir.iterdir())
            assert len(extracted_roots) > 0, f"No cache found in {caches_dir}"
            extracted_root = extracted_roots[0]
        elif blueprints_dir.exists():
            extracted_roots = list(blueprints_dir.iterdir())
            assert len(extracted_roots) > 0, f"No blueprint found in {blueprints_dir}"
            extracted_root = extracted_roots[0]
        else:
            raise AssertionError("Downloaded cache not found in expected locations")

        assert extracted_root.exists()
        assert (extracted_root / "blueprint" / "vcfcache.bcf").exists()
        annotation_dir = extracted_root / "cache" / "test_anno"
        assert annotation_dir.is_dir()

        # Annotate a small sample using the pulled cache
        out_anno = tmp_path / "annotated"
        _run(
            [
                "vcfcache",
                "annotate",
                "-a",
                str(annotation_dir),
                "--vcf",
                str(vcfcache_root / "tests" / "data" / "nodata" / "sample4.bcf"),
                "--output",
                str(out_anno),
                "-y",
                str(params_file),
                "--force",
            ],
            env=env,
            cwd=work_dir,
        )
        produced = out_anno / "sample4_vc.bcf"
        assert produced.exists()

    finally:
        if dep_id and token:
            _delete_sandbox_deposition(dep_id, token)
        shutil.rmtree(work_dir, ignore_errors=True)
