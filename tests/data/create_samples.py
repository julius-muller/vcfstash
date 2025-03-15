import subprocess
from pathlib import Path

# Create a temporary working directory
TEST_DATA = Path("tests/data/nodata")
GNOMAD_TEST = TEST_DATA / "gnomad_test.vcf"
DBSNP_TEST = TEST_DATA / "dbsnp_test.vcf"
ref_fa = TEST_DATA / "reference.fasta"

subprocess.run([
    "bcftools", "norm",
    "-m-", "-c", "x",  # split multiallelic sites
    "-f", str(TEST_DATA / "reference.fasta"),  # reference FASTA
    "-o", str(GNOMAD_TEST).replace('.vcf','.bcf'), "-Ob", '--write-index',
    str(GNOMAD_TEST)
], check=True)

subprocess.run([
    "bcftools", "norm",
    "-m-", "-c", "x",  # split multiallelic sites
    "-f", str(TEST_DATA / "reference.fasta"),  # reference FASTA
    "-o", str(DBSNP_TEST ).replace('.vcf','.bcf'), "-Ob", '--write-index',
    str(DBSNP_TEST )
], check=True)

