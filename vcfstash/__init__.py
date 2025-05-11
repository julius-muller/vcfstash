"""VCF Annotation Cache package.

VCFstash is a tool to accelerate VCF annotations of large VCF files by maintaining
a cache of frequently shared variants across human WGS samples. It manages a variant
cache database and runs VCF annotations only on novel variants not present in the cache,
significantly reducing annotation time.
"""
