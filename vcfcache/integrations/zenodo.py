# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Minimal Zenodo REST client helpers for cache upload/download."""

import os
import threading
from contextlib import contextmanager
from pathlib import Path
from typing import Optional

import requests

DEFAULT_HTTP_TIMEOUT = (10, 30)  # (connect, read) seconds
DEFAULT_ZENODO_HARD_TIMEOUT_S = 60  # wall-clock seconds (best-effort; prevents DNS stalls)

ZENODO_API = "https://zenodo.org/api"
ZENODO_SANDBOX_API = "https://sandbox.zenodo.org/api"


class ZenodoError(RuntimeError):
    pass


def _zenodo_hard_timeout_seconds() -> int:
    raw = os.environ.get("VCFCACHE_ZENODO_HARD_TIMEOUT", "").strip()
    if not raw:
        return DEFAULT_ZENODO_HARD_TIMEOUT_S
    try:
        seconds = int(raw)
    except ValueError:
        raise ZenodoError(
            f"Invalid VCFCACHE_ZENODO_HARD_TIMEOUT={raw!r}; must be an integer seconds value."
        )
    return max(0, seconds)


@contextmanager
def _hard_timeout(seconds: int, *, message: str):
    """Best-effort wall-clock timeout for Zenodo calls.

    Requests' connect/read timeouts do not reliably bound DNS resolution time in
    some environments. On POSIX and the main thread, we use SIGALRM to interrupt
    long waits; elsewhere, this becomes a no-op.
    """
    if (
        seconds <= 0
        or os.name != "posix"
        or threading.current_thread() is not threading.main_thread()
    ):
        yield
        return

    import signal

    old_handler = signal.getsignal(signal.SIGALRM)

    def _handler(_signum, _frame):
        raise ZenodoError(message)

    signal.signal(signal.SIGALRM, _handler)
    signal.setitimer(signal.ITIMER_REAL, float(seconds))
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, old_handler)


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
    sandbox: bool = False,
    min_size_mb: float = 1.0,
) -> list:
    """Search Zenodo for vcfcache blueprints or caches.

    Args:
        item_type: Type of item to search for ("blueprints" or "caches")
        genome: Optional genome build filter (e.g., "GRCh38", "GRCh37")
        source: Optional data source filter (e.g., "gnomad")
        sandbox: If True, search sandbox Zenodo

    Returns:
        List of record dictionaries with blueprints
    """
    def _search(query: str) -> list[dict]:
        search_url = f"{_api_base(sandbox)}/records/"
        resp = None
        last_err: Exception | None = None
        for attempt in range(3):
            try:
                hard = _zenodo_hard_timeout_seconds()
                with _hard_timeout(
                    hard,
                    message=(
                        f"Zenodo request exceeded {hard}s (DNS/network stall?). "
                        "Try again, or use `vcfcache list --local`."
                    ),
                ):
                    resp = requests.get(
                        search_url,
                        params={"q": query, "size": 25},
                        timeout=DEFAULT_HTTP_TIMEOUT,
                    )
                last_err = None
                break
            except requests.exceptions.RequestException as e:
                last_err = e
                if attempt < 2:
                    import time
                    time.sleep(1 + attempt)
        if resp is None:
            raise ZenodoError(f"Failed to search Zenodo after retries: {last_err}")

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

            files = hit.get("files", [])
            total_size = sum(f.get("size", 0) for f in files)
            size_mb = total_size / (1024 * 1024)
            if size_mb < min_size_mb:
                continue

            records.append(
                {
                    "title": metadata.get("title", "Unknown"),
                    "doi": hit.get("doi", "Unknown"),
                    "created": metadata.get("publication_date", hit.get("created", "Unknown")),
                    "description": metadata.get("description", ""),
                    "keywords": metadata.get("keywords", []),
                    "size_mb": size_mb,
                    "creators": metadata.get("creators", []),
                }
            )
        return records

    def _matches_item_type(record: dict) -> bool:
        keywords_raw = record.get("keywords") or []
        keywords = [k.lower() for k in keywords_raw if isinstance(k, str)]
        title = str(record.get("title", "")).lower()
        desc = str(record.get("description", "")).lower()

        has_bp_kw = any(k.startswith("bp-") or k == "blueprint" for k in keywords)
        has_cache_kw = any(k.startswith("cache-") or k == "cache" for k in keywords)

        title_has_bp = "blueprint" in title
        title_has_cache = "cache" in title

        if item_type == "blueprints":
            if has_cache_kw:
                return False
            if has_bp_kw:
                return True
            return title_has_bp and not title_has_cache

        # caches
        if has_bp_kw:
            return False
        if has_cache_kw:
            return True
        return title_has_cache and not title_has_bp

    # Keyword-based search (keywords are required for vcfcache records).
    query_parts = ["keywords:vcfcache"]
    if item_type == "blueprints":
        query_parts.append("keywords:blueprint")
    else:
        query_parts.append("keywords:cache")
    if genome:
        query_parts.append(f"keywords:{genome}")
    if source:
        query_parts.append(f"keywords:{source}")
    primary_query = " AND ".join(query_parts)

    try:
        primary = _search(primary_query)
        records = [r for r in primary if _matches_item_type(r)]
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
        hard = _zenodo_hard_timeout_seconds()
        with _hard_timeout(
            hard,
            message=(
                f"Zenodo request exceeded {hard}s (DNS/network stall?). "
                "Try again, or pass a DOI directly."
            ),
        ):
            resp = requests.get(
                search_url,
                params={"q": query, "size": 1},
                timeout=DEFAULT_HTTP_TIMEOUT,
            )
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
