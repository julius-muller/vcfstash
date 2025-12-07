import os
from pathlib import Path

import pytest

from vcfcache.integrations.zenodo import download_doi


@pytest.mark.integration
def test_download_doi_if_configured(tmp_path):
    doi = os.environ.get("ZENODO_DOI_DOWNLOAD")
    if not doi:
        pytest.skip("ZENODO_DOI_DOWNLOAD not set; skipping real download")

    dest = tmp_path / "cache.tar.gz"
    download_doi(doi, dest)
    assert dest.exists()
    assert dest.stat().st_size > 0

