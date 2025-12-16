#!/bin/bash
# Debug script to check variant counts at each step

CACHED_DIR="$1"
UNCACHED_DIR="$2"

if [ -z "$CACHED_DIR" ] || [ -z "$UNCACHED_DIR" ]; then
    echo "Usage: $0 <cached_dir> <uncached_dir>"
    echo "Example: $0 ft1/cached ft1/uncached"
    exit 1
fi

echo "=========================================="
echo "Variant Count Analysis"
echo "=========================================="
echo ""

echo "Input file (should be same for both):"
echo "  Original: $(bcftools index -n ft1.bcf 2>/dev/null || echo 'N/A')"
echo ""

echo "Uncached workflow:"
UNCACHED_WORK="$UNCACHED_DIR/work/annotate-nocache"
if [ ! -d "$UNCACHED_WORK" ]; then
    UNCACHED_WORK="$UNCACHED_DIR/work"
fi
if [ -f "$UNCACHED_WORK/"*"_normalized.bcf" ]; then
    NORM=$(ls "$UNCACHED_WORK/"*"_normalized.bcf" 2>/dev/null | head -1)
    echo "  After norm:  $(bcftools index -n "$NORM" 2>/dev/null || echo 'N/A')"
fi
# Try both _vc.bcf and _vst.bcf naming
UNCACHED_OUT=$(ls "$UNCACHED_DIR/"*"_vc.bcf" "$UNCACHED_DIR/"*"_vst.bcf" 2>/dev/null | head -1)
if [ -n "$UNCACHED_OUT" ]; then
    echo "  Final output: $(bcftools index -n "$UNCACHED_OUT" 2>/dev/null || echo 'N/A')"
else
    echo "  Final output: N/A"
fi
echo ""

echo "Cached workflow:"
CACHED_WORK="$CACHED_DIR/work/annotate"
if [ ! -d "$CACHED_WORK" ]; then
    CACHED_WORK="$CACHED_DIR/work"
fi
if [ -f "$CACHED_WORK/"*"_normalized.bcf" ]; then
    NORM=$(ls "$CACHED_WORK/"*"_normalized.bcf" 2>/dev/null | head -1)
    echo "  After norm:  $(bcftools index -n "$NORM" 2>/dev/null || echo 'N/A')"
fi
# Try both naming conventions for intermediate files
for pattern in "_isecvst.bcf" "_isec_vc.bcf"; do
    if [ -f "$CACHED_WORK/"*"$pattern" ]; then
        ISEC=$(ls "$CACHED_WORK/"*"$pattern" 2>/dev/null | head -1)
        echo "  Step 1 (with cache): $(bcftools index -n "$ISEC" 2>/dev/null || echo 'N/A')"
        break
    fi
done
for pattern in "_isecvst_miss.bcf" "_isec_vc_miss.bcf"; do
    if [ -f "$CACHED_WORK/"*"$pattern" ]; then
        MISS=$(ls "$CACHED_WORK/"*"$pattern" 2>/dev/null | head -1)
        echo "  Step 2 (missing):    $(bcftools index -n "$MISS" 2>/dev/null || echo 'N/A')"
        break
    fi
done
if [ -f "$CACHED_WORK/"*"_missing_annotated.bcf" ]; then
    ANNO=$(ls "$CACHED_WORK/"*"_missing_annotated.bcf" 2>/dev/null | head -1)
    echo "  Step 3 (annotated):  $(bcftools index -n "$ANNO" 2>/dev/null || echo 'N/A')"
fi
# Try both _vc.bcf and _vst.bcf naming
CACHED_OUT=$(ls "$CACHED_DIR/"*"_vc.bcf" "$CACHED_DIR/"*"_vst.bcf" 2>/dev/null | head -1)
if [ -n "$CACHED_OUT" ]; then
    echo "  Final output: $(bcftools index -n "$CACHED_OUT" 2>/dev/null || echo 'N/A')"
else
    echo "  Final output: N/A"
fi
echo ""

if [ -n "$CACHED_OUT" ] && [ -n "$UNCACHED_OUT" ]; then
    CACHED_COUNT=$(bcftools index -n "$CACHED_OUT" 2>/dev/null || echo 0)
    UNCACHED_COUNT=$(bcftools index -n "$UNCACHED_OUT" 2>/dev/null || echo 0)
    echo "Difference: $((CACHED_COUNT - UNCACHED_COUNT)) variants"
else
    echo "Difference: Cannot calculate (output files not found)"
fi
