# VCFcache Testing (quick ref)

This short file stays for legacy links. The canonical, detailed guide is `tests/README.md`.

## Quick commands

### Host (vanilla)
```bash
python -m pytest tests/ \
  --color=yes --disable-warnings -rA --durations=10 -vv --maxfail=1
```

### Blueprint image
```bash
docker run --rm -t \
  --entrypoint /bin/bash \
  ghcr.io/julius-muller/vcfcache-blueprint:TAG \
  -lc 'cd /app && export PYTHONPATH=/app/venv/lib/python3.13/site-packages && PYTEST_ADDOPTS="--color=yes --disable-warnings -rA --durations=10 -vv --maxfail=1" \
       python3 -m pytest tests'
```

### Annotated image (mount your VEP cache)
```bash
docker run --rm -t \
  --entrypoint /bin/bash \
  -v /path/to/vep/cache:/opt/vep/.vep:ro \
  ghcr.io/julius-muller/vcfcache-annotated:TAG \
  -lc 'cd /app && export PYTHONPATH=/app/venv/lib/python3.13/site-packages && PYTEST_ADDOPTS="--color=yes --disable-warnings -rA --durations=10 -vv --maxfail=1" \
       VCFCACHE_LOGLEVEL=INFO VCFCACHE_FILE_LOGLEVEL=ERROR \
       python3 -m pytest tests'
```

Add `-s` inside `PYTEST_ADDOPTS` if you want live stdout.

## Scenarios (auto-detected)
- **vanilla:** no `/cache`, no VEP
- **blueprint:** `/cache` present, no VEP
- **annotated:** `/cache` present, VEP available (mount required at runtime)

For fixtures, test categories, and behavior per scenario, read `tests/README.md`.
