import os
import shutil
import subprocess
import tempfile
from pathlib import Path
import pytest


def run_cmd(cmd, cwd=None, env=None):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd, env=env)
    assert result.returncode == 0, f"Command failed: {cmd}\nSTDOUT:{result.stdout}\nSTDERR:{result.stderr}"
    return result.stdout


@pytest.mark.skipif(
    os.environ.get("VCFCACHE_RUN_INTEGRATION", "0") != "1",
    reason="Integration test requires network and VCFCACHE_RUN_INTEGRATION=1",
)
def test_full_integration_annotation():
    """
    Integration test: download cache via alias + manifest, run annotate, and validate output exists.

    This test is enabled by default and relies on network access to Zenodo.
    It uses the public cache alias cache-hg38-gnomad-4.1wgs-AF0100-vep-115.2-basic.
    """

    alias = "cache-hg38-gnomad-4.1wgs-AF0100-vep-115.2-basic"
    sample_vcf = Path(__file__).resolve().parent / "data" / "nodata" / "sample4.bcf"
    params = Path(__file__).resolve().parent / "config" / "test_params.yaml"

    outdir = Path(tempfile.mkdtemp(prefix="vcfcache_integration_"))
    home = outdir / "home"
    output_bcf = outdir / "sample4_vc.bcf"
    stats_dir = outdir / "stats"
    home.mkdir(parents=True, exist_ok=True)

    try:
        cmd = (
            f"vcfcache annotate -a {alias} "
            f"--vcf {sample_vcf} "
            f"--output {output_bcf} "
            f"--stats-dir {stats_dir} "
            f"-y {params} "
            f"--force "
        )
        env = os.environ.copy()
        env["HOME"] = str(home)
        run_cmd(cmd, env=env)

        # Validate output exists
        assert output_bcf.exists(), f"Annotated BCF missing: {output_bcf}"

        # bcftools sanity check
        run_cmd(f"bcftools view -h {output_bcf}")

    finally:
        shutil.rmtree(outdir, ignore_errors=True)
