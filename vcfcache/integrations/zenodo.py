"""Minimal Zenodo REST client helpers for cache upload/download."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import requests

ZENODO_API = "https://zenodo.org/api"
ZENODO_SANDBOX_API = "https://sandbox.zenodo.org/api"


class ZenodoError(RuntimeError):
    pass


def _api_base(sandbox: bool) -> str:
    return ZENODO_SANDBOX_API if sandbox else ZENODO_API


def _auth_headers(token: Optional[str]) -> dict:
    return {"Authorization": f"Bearer {token}"} if token else {}


def download_doi(doi: str, dest: Path) -> Path:
    """Download the first file of a Zenodo record given a DOI.

    Note: For simplicity we pick the first attached file. Records intended for
    vcfcache caches should contain a single tarball.
    """

    rec_id = doi.split(".")[-1] if "zenodo" in doi else doi
    url = f"{ZENODO_API}/records/{rec_id}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    record = resp.json()
    files = record.get("files", [])
    if not files:
        raise ZenodoError(f"No files found in record {doi}")
    file_url = files[0]["links"]["self"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(file_url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    f.write(chunk)
    return dest


def create_deposit(token: str, sandbox: bool = False) -> dict:
    url = f"{_api_base(sandbox)}/deposit/depositions"
    resp = requests.post(url, params={"access_token": token}, json={}, timeout=30)
    resp.raise_for_status()
    return resp.json()


def upload_file(
    deposition: dict, path: Path, token: str, sandbox: bool = False
) -> dict:
    bucket = deposition["links"]["bucket"]
    filename = path.name
    with open(path, "rb") as fp:
        resp = requests.put(
            f"{bucket}/{filename}",
            data=fp,
            params={"access_token": token},
            timeout=120,
        )
    resp.raise_for_status()
    return resp.json()


def publish_deposit(deposition: dict, token: str, sandbox: bool = False) -> dict:
    url = deposition["links"]["publish"]
    resp = requests.post(url, params={"access_token": token}, timeout=30)
    resp.raise_for_status()
    return resp.json()
