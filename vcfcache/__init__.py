# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""VCF Annotation Cache package.

VCFcache is a tool to accelerate VCF annotations of large VCF files by maintaining
a cache of frequently shared variants across human WGS samples. It manages a variant
cache database and runs VCF annotations only on novel variants not present in the cache,
significantly reducing annotation time.
"""

__version__ = "0.4.0b2"
EXPECTED_BCFTOOLS_VERSION = "1.20"