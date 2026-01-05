# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Naming helpers for blueprint/cache aliases.

This centralises the alias formats we use across manifests, Docker tags, and
Zenodo deposits so tests and CLI stay consistent.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class CacheName:
    genome: str
    source: str
    release: str
    filt: str  # e.g. AF0100
    tool: str | None = None
    tool_version: str | None = None
    preset: str | None = None

    @property
    def is_cache(self) -> bool:
        return all([self.tool, self.tool_version, self.preset])

    @property
    def cache_alias(self) -> str:
        if not self.is_cache:
            raise ValueError("tool, tool_version, preset required for cache alias")
        return (
            f"cache-{self.genome}-{self.source}-{self.release}-"
            f"{self.filt}-{self.tool}-{self.tool_version}-{self.preset}"
        )

    @property
    def blueprint_alias(self) -> str:
        return f"bp-{self.genome}-{self.source}-{self.release}-{self.filt}"

    @classmethod
    def parse(cls, alias: str) -> "CacheName":
        """Parse either blueprint or cache alias into components."""
        if alias.startswith("cache-"):
            parts = alias.split("-")
            if len(parts) < 8:
                raise ValueError(f"Invalid cache alias: {alias}")
            _, genome, source, release, filt, tool, tool_version, *rest = parts
            preset = "-".join(rest) if rest else None
            return cls(
                genome=genome,
                source=source,
                release=release,
                filt=filt,
                tool=tool,
                tool_version=tool_version,
                preset=preset,
            )
        if alias.startswith("bp-"):
            _, genome, source, release, filt = alias.split("-", 4)
            return cls(genome=genome, source=source, release=release, filt=filt)
        raise ValueError(f"Unknown alias format: {alias}")
