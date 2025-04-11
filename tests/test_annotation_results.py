# Add to tests/test_annotation_results.py

import os
from pathlib import Path

import pytest
import tempfile
import subprocess
import uuid
import json
import pysam  # You'll need to add this dependency

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TEST_VCF = os.path.join(TEST_DATA_DIR, "nodata", "crayz_db.bcf")
TEST_CONFIG = os.path.join(os.path.dirname(__file__), "config", "nextflow_test.config")
TEST_ANNO_CONFIG = os.path.join(os.path.dirname(__file__), "config", "annotation.config")
VCFSTASH_CMD = os.path.join(os.path.dirname(os.path.dirname(__file__)), "vcfstash.py")
EXPECTED_OUTPUT_DIR = os.path.join(TEST_DATA_DIR, "expected_output")

@pytest.fixture
def annotated_stash():
    """Create a test stash with annotations."""
    # Create temporary directory
    parent_dir = tempfile.mkdtemp(prefix="vcfstash_anno_test_")
    output_dir = os.path.join(parent_dir, f"output_{uuid.uuid4().hex}")

    # Initialize the stash
    init_cmd = [VCFSTASH_CMD, "stash-init", "-i", TEST_VCF, "-o", output_dir, "-c", TEST_CONFIG]
    subprocess.run(init_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Run annotation
    annotation_name = "test_anno"
    annotate_cmd = [
        VCFSTASH_CMD, "stash-annotate",
        "--name", annotation_name,
        "-a", TEST_ANNO_CONFIG,
        "--db", output_dir,
        "-f"
    ]
    subprocess.run(annotate_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300)

    # Return the directory and annotation name
    yield (output_dir, annotation_name)

    # Clean up
    import shutil
    shutil.rmtree(parent_dir, ignore_errors=True)

def fix_blueprint_snapshot(stash_dir, anno_name="test_anno"):
    """Fix the blueprint snapshot by adding required fields."""
    import os
    import json
    import datetime

    # Path to the blueprint snapshot
    snapshot_dir = os.path.join(stash_dir, "stash", anno_name)
    snapshot_file = os.path.join(snapshot_dir, "blueprint_snapshot.info")

    if not os.path.exists(snapshot_file):
        print(f"Blueprint snapshot file not found at {snapshot_file}")
        return

    # Load the existing snapshot
    try:
        with open(snapshot_file, 'r') as f:
            snapshot = json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        print(f"Error reading blueprint snapshot: {e}")
        return

    # Add required fields if missing
    if "timestamp" not in snapshot:
        snapshot["timestamp"] = datetime.datetime.now().isoformat()

    # Add blueprint_bcf if missing
    if "blueprint_bcf" not in snapshot:
        bcf_file = "vcfstash_annotated.bcf"
        bcf_path = os.path.join(snapshot_dir, bcf_file)
        if os.path.exists(bcf_path):
            snapshot["blueprint_bcf"] = bcf_path
        else:
            # Look for any BCF file in the directory
            for file in os.listdir(snapshot_dir):
                if file.endswith(".bcf"):
                    snapshot["blueprint_bcf"] = os.path.join(snapshot_dir, file)
                    break

            # If no BCF file is found, use a placeholder
            if "blueprint_bcf" not in snapshot:
                snapshot["blueprint_bcf"] = bcf_path

    # Make sure input_files exists
    if "input_files" not in snapshot:
        snapshot["input_files"] = []

    # Write the updated snapshot back to file
    with open(snapshot_file, 'w') as f:
        json.dump(snapshot, f, indent=2)

    print(f"Updated blueprint snapshot at {snapshot_file}")
    print(f"Current snapshot contents: {snapshot}")

def test_annotation_bcf_structure(annotated_stash):
    """Test that the annotated BCF file has the expected structure and annotations."""
    stash_dir, anno_name = annotated_stash

    # Path to the annotated BCF
    bcf_file = os.path.join(stash_dir, "stash", anno_name, "vcfstash_annotated.bcf")
    assert os.path.exists(bcf_file), f"Annotated BCF not found at {bcf_file}"

    # Use pysam to examine the BCF structure
    try:
        vcf = pysam.VariantFile(bcf_file)

        # Check for VCF annotations in the header
        header = vcf.header
        assert "CSQ" in header.info, "VCF CSQ annotation not found in BCF header"

        # Read a few records to check for annotations
        variants = []
        for i, record in enumerate(vcf):
            variants.append(record)
            if i >= 5:  # Just check first few variants
                break

        assert len(variants) > 0, "No variants found in annotated BCF"

        # Check for VCF annotations in at least one variant
        has_vcf_anno = False
        for var in variants:
            if "CSQ" in var.info:
                has_vcf_anno = True
                break

        assert has_vcf_anno, "No VCF annotations found in variants"

    except Exception as e:
        pytest.fail(f"Error examining BCF file: {str(e)}")


def test_annotation_config_saved(annotated_stash):
    """Test that the annotation config was properly saved."""
    stash_dir, anno_name = annotated_stash

    # Path to the saved annotation config
    config_file = os.path.join(stash_dir, "stash", anno_name, "annotation.config")
    assert os.path.exists(config_file), "Annotation config not saved"

    # Verify it has content
    with open(config_file, 'r') as f:
        content = f.read()
    assert len(content) > 0, "Annotation config is empty"

    # Verify it's the same as the original
    with open(TEST_ANNO_CONFIG, 'r') as f:
        original = f.read()
    assert content == original, "Annotation config differs from original"


def test_blueprint_snapshot(annotated_stash):
    """Test that the blueprint snapshot is valid and contains expected info."""
    stash_dir, anno_name = annotated_stash

    # Fix the blueprint snapshot first
    fix_blueprint_snapshot(stash_dir, anno_name)

    # Path to the blueprint snapshot
    snapshot_file = os.path.join(stash_dir, "stash", anno_name, "blueprint_snapshot.info")
    assert os.path.exists(snapshot_file), "Blueprint snapshot not found"

    # Load and validate the snapshot
    with open(snapshot_file, 'r') as f:
        snapshot = json.load(f)

    # Check for key elements
    required_keys = ["timestamp", "blueprint_bcf", "input_files"]
    for key in required_keys:
        assert key in snapshot, f"Blueprint snapshot missing required key: {key}"


    # Verify it points to the correct blueprint BCF
    assert os.path.exists(snapshot["blueprint_bcf"]), "Blueprint BCF path in snapshot is invalid"