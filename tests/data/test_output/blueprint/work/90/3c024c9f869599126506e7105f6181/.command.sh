#!/bin/bash -ue
# First, convert to BCF and index to ensure proper format
/home/j380r/projects/vcfstash/tools/bcftools view -G -Ou crayz_db.bcf --threads 0 |
/home/j380r/projects/vcfstash/tools/bcftools annotate -x INFO --rename-chrs chr_add.txt --threads 0 -Ou |
/home/j380r/projects/vcfstash/tools/bcftools filter -i 'CHROM ~ "^chr[1-9,X,Y,M]" && CHROM ~ "[0-9,X,Y,M]$"' --threads 0 -Ou | \
 \
/home/j380r/projects/vcfstash/tools/bcftools norm -m- -c x -f reference.fasta -o crayz_db_norm.bcf -Ob --threads 0 --write-index
