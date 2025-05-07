#!/bin/bash
# Simple VCF/BCF Reference Validator
# Checks chromosome compatibility and reference sequence integrity

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Print functions
print_step() { echo -e "${BLUE}→ ${1}${NC}"; }
print_success() { echo -e "${GREEN}✓ ${1}${NC}"; }
print_warning() { echo -e "${YELLOW}⚠ ${1}${NC}"; }
print_error() { echo -e "${RED}${BOLD}✗ ${1}${NC}"; exit 1; }

# Timer functions
timer_start() { start_time=$(date +%s.%N); }
timer_end() { duration=$(echo "$(date +%s.%N) - $start_time" | bc); echo -e "  ${YELLOW}Time: ${duration}s${NC}"; }

# Usage info
if [ $# -lt 3 ]; then
    print_error "Missing required arguments\nUsage: $0 <VCF_FILE> <REFERENCE_FASTA> <CHR_MAPPING> [VARIANTS_TO_CHECK]"
fi

VCF_FILE=$1
REF_FASTA=$2
CHR_MAP=$3
VARIANTS_TO_CHECK=${4:-5}

# Verify files exist
print_step "Checking input files..."
timer_start
[ -f "$VCF_FILE" ] || print_error "VCF/BCF file not found: $VCF_FILE"
[ -f "$REF_FASTA" ] || print_error "Reference genome not found: $REF_FASTA"
[ -f "$CHR_MAP" ] || print_error "Chromosome mapping file not found: $CHR_MAP"

# Check for index files
[ -f "${VCF_FILE}.csi" ] || [ -f "${VCF_FILE}.tbi" ] || print_error "VCF/BCF index not found. Please index your file."
[ -f "${REF_FASTA}.fai" ] || print_error "Reference index not found. Please run 'samtools faidx' on your reference."
print_success "All input files exist"
print_success "Found both VCF and reference indexes"
timer_end

# Get reference chromosomes from index
print_step "Extracting reference genome information..."
timer_start
REF_CHROMS=$(cut -f1 "${REF_FASTA}.fai" | sort)
REF_CHROM_COUNT=$(echo "$REF_CHROMS" | wc -l)
echo "  Found $REF_CHROM_COUNT chromosomes in reference genome"
timer_end

# Load reference index into memory
print_step "Loading reference index..."
timer_start
declare -A REF_INDEX
while IFS=$'\t' read -r chrom length offset linelen linewidth; do
    REF_INDEX["$chrom"]=$length
done < "${REF_FASTA}.fai"
echo "  Loaded ${#REF_INDEX[@]} chromosomes from reference index"
timer_end

# Get VCF/BCF chromosomes
print_step "Extracting VCF chromosome information..."
timer_start
VCF_INDEX="${VCF_FILE}.csi"
[ -f "$VCF_INDEX" ] || VCF_INDEX="${VCF_FILE}.tbi"
VCF_CHROMS=$(bcftools index -s "$VCF_FILE" 2>/dev/null | cut -f1 | sort)
VCF_CHROM_COUNT=$(echo "$VCF_CHROMS" | wc -l)
echo "  Retrieved chromosomes from VCF index"
echo "  Found $VCF_CHROM_COUNT unique chromosomes in VCF"
timer_end

# Load chromosome mapping
print_step "Processing chromosome mappings..."
timer_start
declare -A CHR_MAPPING_FORWARD
declare -A CHR_MAPPING_REVERSE
while IFS=$'\t' read -r from to; do
    CHR_MAPPING_FORWARD["$from"]="$to"
    CHR_MAPPING_REVERSE["$to"]="$from"
done < "$CHR_MAP"
echo "  Loaded ${#CHR_MAPPING_FORWARD[@]} chromosome mappings from $CHR_MAP"
timer_end

# Check chromosome compatibility
print_step "Checking chromosome compatibility..."
timer_start
INCOMPATIBLE=0
FIRST_VCF_CHROM=$(echo "$VCF_CHROMS" | head -1)
FIRST_REF_CHROM=$(echo "$REF_CHROMS" | head -1)

# Determine naming conventions
if [[ "$FIRST_VCF_CHROM" == chr* && "$FIRST_REF_CHROM" != chr* ]]; then
    echo "  VCF uses 'chr' prefix, reference doesn't"
    MAPPING_DIRECTION="vcf_to_ref"
elif [[ "$FIRST_VCF_CHROM" != chr* && "$FIRST_REF_CHROM" == chr* ]]; then
    echo "  Reference uses 'chr' prefix, VCF doesn't"
    MAPPING_DIRECTION="ref_to_vcf"
else
    echo "  Both VCF and reference use the same chromosome naming convention"
    MAPPING_DIRECTION="none"
fi

# Check each VCF chromosome exists in reference (with mapping if needed)
for chrom in $VCF_CHROMS; do
    REF_CHROM="$chrom"  # Default to no mapping

    # Apply mapping based on direction
    if [[ "$MAPPING_DIRECTION" == "vcf_to_ref" ]]; then
        if [[ -n "${CHR_MAPPING_REVERSE[$chrom]}" ]]; then
            REF_CHROM="${CHR_MAPPING_REVERSE[$chrom]}"
        else
            # Try removing chr prefix
            REF_CHROM="${chrom#chr}"
        fi
    elif [[ "$MAPPING_DIRECTION" == "ref_to_vcf" ]]; then
        if [[ -n "${CHR_MAPPING_FORWARD[$chrom]}" ]]; then
            REF_CHROM="${CHR_MAPPING_FORWARD[$chrom]}"
        else
            # Try adding chr prefix
            REF_CHROM="chr${chrom}"
        fi
    fi

    # Check if mapped chromosome exists in reference
    if [[ -z "${REF_INDEX[$REF_CHROM]}" ]]; then
        echo "  ❌ Chromosome $chrom (mapped to $REF_CHROM) not found in reference"
        ((INCOMPATIBLE++))
    fi
done

if [[ $INCOMPATIBLE -eq 0 ]]; then
    print_success "All VCF chromosomes can be mapped to reference"
else
    print_warning "$INCOMPATIBLE VCF chromosomes could not be mapped to reference"
fi
timer_end

# Function to map chromosome names
map_chromosome() {
    local chrom="$1"
    local mapped="$chrom"  # Default to same name

    if [[ "$MAPPING_DIRECTION" == "vcf_to_ref" ]]; then
        if [[ -n "${CHR_MAPPING_REVERSE[$chrom]}" ]]; then
            mapped="${CHR_MAPPING_REVERSE[$chrom]}"
        else
            mapped="${chrom#chr}"
        fi
    elif [[ "$MAPPING_DIRECTION" == "ref_to_vcf" ]]; then
        if [[ -n "${CHR_MAPPING_FORWARD[$chrom]}" ]]; then
            mapped="${CHR_MAPPING_FORWARD[$chrom]}"
        else
            mapped="chr${chrom}"
        fi
    fi

    echo "$mapped"
}

# Check variant reference alleles
print_step "Checking reference alleles in variants..."
timer_start

# Get chromosomes to check (limited to what we can process)
CHROMS_TO_CHECK=$(echo "$VCF_CHROMS" | head -$VARIANTS_TO_CHECK)
VARIANTS_CHECKED=0
VARIANTS_MATCHED=0

for chrom in $CHROMS_TO_CHECK; do
    echo "  Checking chromosome: $chrom"

    # Get first variant from this chromosome
    VARIANT=$(bcftools view -H -r "$chrom" "$VCF_FILE" 2>/dev/null | head -1)

    if [[ -z "$VARIANT" ]]; then
        echo "  No variants found on chromosome $chrom"
        continue
    fi

    # Parse variant
    VCF_CHROM=$(echo "$VARIANT" | cut -f1)
    POS=$(echo "$VARIANT" | cut -f2)
    ID=$(echo "$VARIANT" | cut -f3)
    REF_ALLELE=$(echo "$VARIANT" | cut -f4)

    # Skip if reference allele is too large
    if [[ ${#REF_ALLELE} -gt 1000 ]]; then
        echo "  Skipping large reference allele (${#REF_ALLELE} bp) at $VCF_CHROM:$POS"
        continue
    fi

    # Map chromosome to reference
    REF_CHROM=$(map_chromosome "$VCF_CHROM")

    echo "  Checking variant at $VCF_CHROM:$POS ($ID) with REF=$REF_ALLELE"

    # Check if chromosome exists in reference
    if [[ -z "${REF_INDEX[$REF_CHROM]}" ]]; then
        echo "  ❌ Chromosome $VCF_CHROM (mapped to $REF_CHROM) not found in reference"
        continue
    fi

    # Extract reference sequence
    REF_SEQ=$(samtools faidx "$REF_FASTA" "${REF_CHROM}:${POS}-$((POS + ${#REF_ALLELE} - 1))" 2>/dev/null | grep -v "^>" | tr -d '\n')

    ((VARIANTS_CHECKED++))

    if [[ "${REF_ALLELE,,}" == "${REF_SEQ,,}" ]]; then
        echo -e "  ${GREEN}✓ Reference allele matches: $REF_SEQ${NC}"
        ((VARIANTS_MATCHED++))
    else
        echo -e "  ${RED}✗ Reference mismatch at $VCF_CHROM:$POS (VCF: $REF_ALLELE, FASTA: $REF_SEQ)${NC}"
    fi
done

if [[ $VARIANTS_CHECKED -eq 0 ]]; then
    print_warning "Could not extract any variants for reference validation"
elif [[ $VARIANTS_MATCHED -eq $VARIANTS_CHECKED ]]; then
    print_success "All $VARIANTS_CHECKED checked variants match the reference genome"
else
    print_warning "$VARIANTS_MATCHED out of $VARIANTS_CHECKED variants match the reference genome"
fi
timer_end

# Final summary
if [[ $INCOMPATIBLE -eq 0 && ($VARIANTS_MATCHED -eq $VARIANTS_CHECKED || $VARIANTS_CHECKED -eq 0) ]]; then
    print_success "Validation completed successfully"
    exit 0
else
    print_warning "Validation completed with warnings"
    exit 0
fi
