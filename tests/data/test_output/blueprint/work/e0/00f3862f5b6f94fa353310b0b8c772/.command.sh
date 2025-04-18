#!/bin/bash -ue
# Check if files exist
[ -f "crayz_db.bcf" ] || { echo "VCF file not found: crayz_db.bcf"; exit 1; }
[ -f "crayz_db.bcf.csi" ] || { echo "VCF index not found: crayz_db.bcf.csi"; exit 1; }
[ -f "reference.fasta" ] || { echo "Reference file not found: reference.fasta"; exit 1; }
[ -f "reference.fasta.fai" ] || { echo "Reference index not found: reference.fasta.fai"; exit 1; }
[ -f "chr_add.txt" ] || { echo "Chr add file not found: chr_add.txt"; exit 1; }

# Basic file format checks
/home/j380r/projects/vcfstash/tools/bcftools view -h "crayz_db.bcf" >/dev/null || { echo "Invalid VCF/BCF format"; exit 1; }
