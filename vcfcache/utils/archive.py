# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Archiving helpers for VCFcache caches.

This module bundles cache directories into tarballs and extracts them back.
It's intentionally light-weight (tarfile + hashlib) so it works in both
local and container environments without extra deps.
"""

import hashlib
import tarfile
from pathlib import Path
from typing import Iterable


def tar_cache(cache_dir: Path, tar_path: Path, compression: str = "gz") -> Path:
    """Create a tar archive of a cache directory.

    Args:
        cache_dir: Path to the cache directory containing blueprint/cache.
        tar_path: Output tar path (".tar.gz" recommended).
        compression: tarfile mode suffix ("gz", "bz2", "xz", or "").

    Returns:
        Path to the created tarball.
    """

    cache_dir = cache_dir.expanduser().resolve()
    tar_path = tar_path.expanduser().resolve()
    mode = f"w:{compression}" if compression else "w"

    tar_path.parent.mkdir(parents=True, exist_ok=True)

    with tarfile.open(tar_path, mode) as tf:
        tf.add(cache_dir, arcname=cache_dir.name)

    return tar_path


def tar_cache_subset(
    cache_dir: Path,
    tar_path: Path,
    *,
    include_blueprint: bool = True,
    include_cache_name: str | None = None,
    include_empty_cache_dir: bool = True,
    compression: str = "gz",
) -> Path:
    """Create a tar archive of a cache directory with selective contents.

    Args:
        cache_dir: Base directory containing blueprint/ and cache/.
        tar_path: Output tar path (".tar.gz" recommended).
        include_blueprint: Whether to include blueprint/ in the archive.
        include_cache_name: If set, include only cache/<name>/ contents.
        include_empty_cache_dir: Include an empty cache/ dir when excluding caches.
        compression: tarfile mode suffix ("gz", "bz2", "xz", or "").

    Returns:
        Path to the created tarball.
    """

    cache_dir = cache_dir.expanduser().resolve()
    tar_path = tar_path.expanduser().resolve()
    mode = f"w:{compression}" if compression else "w"

    tar_path.parent.mkdir(parents=True, exist_ok=True)

    base_name = cache_dir.name

    def _filter(tarinfo: tarfile.TarInfo) -> tarfile.TarInfo | None:
        name = tarinfo.name
        if name == base_name:
            return tarinfo
        parts = Path(name).parts
        if not parts or parts[0] != base_name:
            return None
        rel = parts[1:]
        if not rel:
            return tarinfo
        if rel[0] == "blueprint":
            return tarinfo if include_blueprint else None
        if rel[0] == "cache":
            if len(rel) == 1:
                return tarinfo if include_empty_cache_dir else None
            if include_cache_name and rel[1] == include_cache_name:
                return tarinfo
            return None
        # Keep root-level files (e.g., .vcfcache_complete, README/info files).
        return tarinfo

    with tarfile.open(tar_path, mode) as tf:
        tf.add(cache_dir, arcname=cache_dir.name, filter=_filter)

    return tar_path


def extract_cache(tar_path: Path, dest_dir: Path) -> Path:
    """Extract a cache tarball to destination directory.

    Args:
        tar_path: Path to tarball.
        dest_dir: Directory to extract into.

    Returns:
        Path to the extracted cache root (dest_dir/<cache_name>).
    """

    tar_path = tar_path.expanduser().resolve()
    dest_dir = dest_dir.expanduser().resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)

    with tarfile.open(tar_path, "r:*") as tf:
        tf.extractall(dest_dir, filter="data")

    with tarfile.open(tar_path, "r:*") as tf:
        members = [m for m in tf.getmembers() if m.name and not m.name.endswith("/")]
        top = Path(members[0].name).parts[0] if members else ""
    return dest_dir / top


def file_md5(path: Path, chunk: int = 1 << 20) -> str:
    """Compute md5 hash of a file (Zenodo supplies md5)."""

    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def dir_md5(paths: Iterable[Path]) -> str:
    """Compute md5 over multiple files deterministically (path + data)."""

    h = hashlib.md5()
    for p in sorted(paths, key=lambda p: str(p)):
        h.update(str(p).encode())
        with open(p, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                h.update(block)
    return h.hexdigest()
