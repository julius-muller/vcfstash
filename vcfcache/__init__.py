"""VCF Annotation Cache package.

VCFcache is a tool to accelerate VCF annotations of large VCF files by maintaining
a cache of frequently shared variants across human WGS samples. It manages a variant
cache database and runs VCF annotations only on novel variants not present in the cache,
significantly reducing annotation time.
"""

# Package-wide constants
EXPECTED_BCFTOOLS_VERSION = "1.22+htslib-1.22"