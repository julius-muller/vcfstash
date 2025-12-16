# Local Docker build helpers (maintainers)

This directory contains helper scripts for building container images locally.

## Supported: build the runtime image

The maintained runtime image is built from `docker/Dockerfile.vcfcache`.

Build + test + (optionally) push `ghcr.io/julius-muller/vcfcache:latest`:
```bash
./scripts/local-build/build-and-push-final.sh
```

Options:
- `--skip-tests`: skip the container pytest run
- `--skip-push`: don’t push to GHCR (build locally only)
- `--force`: rebuild even if the image already exists

## Publishing versioned images

For versioned tags (`ghcr.io/julius-muller/vcfcache:vX.Y.Z`), use `./scripts/release.sh X.Y.Z`, which will build/tag/push consistently with the Python package release.

## Legacy scripts

The remaining scripts in this folder were used for earlier “blueprint/annotated image” experiments and for generating population BCF subsets. They may still be useful for ad-hoc work, but they are not part of the current supported release pipeline.
