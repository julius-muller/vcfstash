"""Minimal Zenodo REST client helpers for cache upload/download."""

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


def download_doi(doi: str, dest: Path, sandbox: bool = False) -> Path:
    """Download the first file of a Zenodo record given a DOI or record ID.

    Args:
        doi: Zenodo DOI (e.g. 10.5281/zenodo.12345) or record id.
        dest: Local path to write the downloaded file.
        sandbox: If True, target the Zenodo sandbox API.

    Note: For simplicity we pick the first attached file. Records intended for
    vcfcache caches should contain a single tarball.
    """

    rec_id = doi.split(".")[-1] if "zenodo" in doi else doi
    url = f"{_api_base(sandbox)}/records/{rec_id}"
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


def search_zenodo_records(
    item_type: str = "blueprints",
    genome: Optional[str] = None,
    source: Optional[str] = None,
    sandbox: bool = False
) -> list:
    """Search Zenodo for vcfcache blueprints or caches.

    Args:
        item_type: Type of item to search for ("blueprints" or "caches")
        genome: Optional genome build filter (e.g., "GRCh38", "GRCh37")
        source: Optional data source filter (e.g., "gnomad")
        sandbox: If True, search sandbox Zenodo

    Returns:
        List of record dictionaries with metadata
    """
    # Build search query using Elasticsearch query string syntax
    # Search in keywords field using field-specific syntax
    query_parts = ["keywords:vcfcache"]

    if item_type == "blueprints":
        query_parts.append("keywords:blueprint")
    else:
        query_parts.append("keywords:cache")

    if genome:
        query_parts.append(f"keywords:{genome}")
    if source:
        query_parts.append(f"keywords:{source}")

    # Combine with AND to require all terms
    query = " AND ".join(query_parts)

    # Zenodo search API
    search_url = f"{_api_base(sandbox)}/records/"

    try:
        # Note: unauthenticated requests limited to 25 results
        # Could be increased to 100 with authentication if needed
        resp = requests.get(
            search_url,
            params={"q": query, "size": 25},
            timeout=30
        )

        # Better error handling - show actual response
        if not resp.ok:
            error_detail = ""
            try:
                error_detail = f" - {resp.json()}"
            except Exception:
                error_detail = f" - {resp.text[:200]}"
            raise ZenodoError(
                f"Zenodo API error ({resp.status_code}): {resp.reason}{error_detail}"
            )

        data = resp.json()

        records = []
        for hit in data.get("hits", {}).get("hits", []):
            metadata = hit.get("metadata", {})

            # Calculate total size
            files = hit.get("files", [])
            total_size = sum(f.get("size", 0) for f in files)

            records.append({
                "title": metadata.get("title", "Unknown"),
                "doi": hit.get("doi", "Unknown"),
                "created": metadata.get("publication_date", "Unknown"),
                "description": metadata.get("description", ""),
                "keywords": metadata.get("keywords", []),
                "size_mb": total_size / (1024 * 1024),
                "creators": metadata.get("creators", []),
            })

        return records

    except requests.exceptions.RequestException as e:
        raise ZenodoError(f"Failed to search Zenodo: {e}") from e


def resolve_zenodo_alias(
    alias_or_doi: str,
    item_type: str = "caches",
    sandbox: bool = False,
) -> tuple[str, str]:
    """Resolve a cache/blueprint alias to a Zenodo DOI via keywords search.

    Args:
        alias_or_doi: Alias (e.g. cache-hg38-...) or a DOI/record id.
        item_type: "caches" or "blueprints".
        sandbox: If True, query the Zenodo sandbox API.

    Returns:
        (doi, alias)
    """
    # If the user passed a DOI or record id directly, just use it.
    if alias_or_doi.startswith("10.") or "zenodo." in alias_or_doi:
        rec_id = alias_or_doi.split(".")[-1] if "zenodo" in alias_or_doi else alias_or_doi
        return alias_or_doi, f"zenodo-{rec_id}"

    alias = alias_or_doi
    keyword_type = "cache" if item_type == "caches" else "blueprint"
    escaped = alias.replace('"', '\\"')
    query = f'keywords:vcfcache AND keywords:{keyword_type} AND keywords:"{escaped}"'

    search_url = f"{_api_base(sandbox)}/records/"
    try:
        resp = requests.get(search_url, params={"q": query, "size": 1}, timeout=30)
        if not resp.ok:
            error_detail = ""
            try:
                error_detail = f" - {resp.json()}"
            except Exception:
                error_detail = f" - {resp.text[:200]}"
            raise ZenodoError(
                f"Zenodo API error ({resp.status_code}): {resp.reason}{error_detail}"
            )
        data = resp.json()
        hits = data.get("hits", {}).get("hits", [])
        if not hits:
            raise ZenodoError(
                f"Could not resolve alias '{alias}' on Zenodo ({'sandbox' if sandbox else 'production'})."
            )
        doi = hits[0].get("doi")
        if not doi:
            raise ZenodoError(f"Resolved record for '{alias}' has no DOI (not published?).")
        return doi, alias
    except requests.exceptions.RequestException as e:
        raise ZenodoError(f"Failed to resolve alias on Zenodo: {e}") from e
