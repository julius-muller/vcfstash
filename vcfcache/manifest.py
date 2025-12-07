"""Public cache manifest handling.

Defines a small schema for mapping cache aliases to Zenodo DOIs and metadata.
"""

from __future__ import annotations

import dataclasses
import yaml
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclasses.dataclass
class CacheEntry:
    alias: str
    doi: str
    genome: str
    af: str
    tool: str
    image_tag: str
    updated_at: str
    md5: Optional[str] = None
    size: Optional[int] = None
    annotation_yaml_md5: Optional[str] = None


def load_manifest(path_or_url: str) -> List[CacheEntry]:
    """Load manifest from local file or URL.

    For URLs we expect a simple HTTP GET; callers can supply a local path for
    offline use.
    """

    import requests

    if path_or_url.startswith("http://") or path_or_url.startswith("https://"):
        resp = requests.get(path_or_url, timeout=30)
        resp.raise_for_status()
        content = resp.text
    else:
        content = Path(path_or_url).read_text()

    data = yaml.safe_load(content) or []
    entries: List[CacheEntry] = []
    for item in data:
        entries.append(CacheEntry(**item))
    return entries


def find_alias(entries: List[CacheEntry], alias: str) -> Optional[CacheEntry]:
    for e in entries:
        if e.alias == alias:
            return e
    return None


def format_manifest(entries: List[CacheEntry]) -> str:
    lines = ["alias\tgenome\taf\ttool\timage\tdoi\tupdated"]
    for e in entries:
        lines.append(
            f"{e.alias}\t{e.genome}\t{e.af}\t{e.tool}\t{e.image_tag}\t{e.doi}\t{e.updated_at}"
        )
    return "\n".join(lines)
