# GRCh37 + VEP 115.2 (Apptainer) Minimal Recipe

This is a minimal recipe for annotating a VCFcache blueprint with Ensembl VEP 115.2 when VEP is provided via Apptainer/Singularity.

## Assumptions

- VEP runs via:
  - `apptainer exec -B /mnt/data:/mnt/data /mnt/data/apps/ensembl-vep/115/vep.sif vep`
- A GRCh37 VEP cache is available on the host and reachable inside the container (configured via `vep_cache` in `params.yaml`).

## Use

Download your blueprint:
```bash
vcfcache blueprint-init --doi 10.5072/zenodo.414604 -o new --debug
```

Annotate the blueprint:
```bash
vcfcache cache-build \
  --db new/<your_blueprint_dir> \
  --name vep115.2_grch37 \
  -a vcfcache/recipes/grch37_vep115.2_apptainer/annotation.yaml \
  -y vcfcache/recipes/grch37_vep115.2_apptainer/params.yaml \
  --force
```

Notes:
- Update `vcfcache/recipes/grch37_vep115.2_apptainer/params.yaml` to point `vep_cache` to your GRCh37 VEP cache path.
- If your Apptainer `.sif` lives elsewhere, edit `annotation_tool_cmd` / `tool_version_command` accordingly.
