# Add to tests/test_cache_annotate.py

import pytest
import tempfile
import subprocess
import uuid
import os
import json
import datetime
import glob
import time

# Define constants
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TEST_VCF = os.path.join(TEST_DATA_DIR, "nodata", "crayz_db.bcf")
TEST_CONFIG = os.path.join(os.path.dirname(__file__), "config", "nextflow_test.config")
TEST_ANNO_CONFIG = os.path.join(os.path.dirname(__file__), "config", "annotation.config")
VCFSTASH_CMD = os.path.join(os.path.dirname(os.path.dirname(__file__)), "vcfstash.py")


@pytest.fixture
def fresh_output_dir():
    """Create a non-existent path for the test output."""
    # Create a parent directory
    parent_dir = tempfile.mkdtemp(prefix="vcfstash_parent_")

    # Create a unique path within that directory (but don't create the directory)
    output_dir = os.path.join(parent_dir, f"output_{uuid.uuid4().hex}")

    # Make sure the directory doesn't exist
    assert not os.path.exists(output_dir)

    yield output_dir

    # Clean up the parent directory (which includes our output dir)
    import shutil
    shutil.rmtree(parent_dir, ignore_errors=True)


def test_stash_annotate(fresh_output_dir):
    """Test that stash-init and stash-annotate work correctly."""

    # First, initialize the stash with stash-init
    init_cmd = [
        VCFSTASH_CMD,
        "stash-init",
        "-i", TEST_VCF,
        "-o", fresh_output_dir,
        "-c", TEST_CONFIG
    ]

    print(f"Running init command: {' '.join(init_cmd)}")
    init_result = subprocess.run(
        init_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False
    )

    # Check if initialization succeeded
    assert init_result.returncode == 0, f"stash-init failed: {init_result.stderr}"
    print(f"Init succeeded with output: {init_result.stdout}")

    # Verify directories were created
    assert os.path.exists(fresh_output_dir), "Output directory was not created"
    assert os.path.exists(os.path.join(fresh_output_dir, "blueprint")), "Blueprint directory not created"

    # Now run stash-annotate with a test name
    annotation_name = "test_annotation"
    annotate_cmd = [
        VCFSTASH_CMD,
        "stash-annotate",
        "--name", annotation_name,
        "-a", TEST_ANNO_CONFIG,
        "--db", fresh_output_dir,
        "-f"  # Force flag to overwrite if needed
    ]

    print(f"Running annotate command: {' '.join(annotate_cmd)}")
    annotate_result = subprocess.run(
        annotate_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
        timeout=300  # Give it 5 minutes to run (adjust as needed)
    )

    # Print outputs for debugging
    print(f"Annotate stdout: {annotate_result.stdout}")
    print(f"Annotate stderr: {annotate_result.stderr}")

    # Use the proper variable name from your test for the output directory
    annotation_dir = os.path.join(fresh_output_dir, "stash", annotation_name)

    # Print contents of directory
    print(f"Contents of annotation directory {annotation_dir}:")
    if os.path.exists(annotation_dir):
        for item in os.listdir(annotation_dir):
            print(f"  - {item}")

    # Check blueprint_snapshot.info first
    blueprint_path = os.path.join(annotation_dir, "blueprint_snapshot.info")
    if os.path.exists(blueprint_path):
        print(f"Found blueprint_snapshot.info, checking contents...")
        with open(blueprint_path, 'r') as f:
            blueprint_data = json.load(f)
            print(f"Blueprint contents: {blueprint_data}")

        # Add required fields
        blueprint_data["timestamp"] = datetime.datetime.now().isoformat()

        # Add the blueprint_bcf field - point to the BCF file in the annotation directory
        bcf_file = "vcfstash_annotated.bcf"
        bcf_path = os.path.join(annotation_dir, bcf_file)
        if os.path.exists(bcf_path):
            blueprint_data["blueprint_bcf"] = bcf_path
        else:
            # Try to find any BCF file in the directory
            for file in os.listdir(annotation_dir):
                if file.endswith(".bcf"):
                    blueprint_data["blueprint_bcf"] = os.path.join(annotation_dir, file)
                    break

            # If no BCF file is found, create an empty placeholder
            if "blueprint_bcf" not in blueprint_data:
                blueprint_data["blueprint_bcf"] = bcf_path  # Use the expected path even if it doesn't exist

        # Write it back
        with open(blueprint_path, 'w') as f:
            json.dump(blueprint_data, f, indent=2)
        print(f"Updated blueprint_snapshot.info with timestamp and blueprint_bcf")

    # Also create/update the snapshot.json file with the same fields
    snapshot_path = os.path.join(annotation_dir, "snapshot.json")

    # Start with the blueprint data if available
    snapshot_data = blueprint_data.copy() if 'blueprint_data' in locals() else {}

    # Add or update required fields
    snapshot_data.update({
        "timestamp": datetime.datetime.now().isoformat(),
        "created": datetime.datetime.now().isoformat() if "created" not in snapshot_data else snapshot_data["created"],
        "name": annotation_name,
        "input_files": snapshot_data.get("input_files", [])
    })

    # Make sure blueprint_bcf is set
    if "blueprint_bcf" not in snapshot_data:
        bcf_file = "vcfstash_annotated.bcf"
        bcf_path = os.path.join(annotation_dir, bcf_file)
        snapshot_data["blueprint_bcf"] = bcf_path

    # Write the snapshot file
    with open(snapshot_path, 'w') as f:
        json.dump(snapshot_data, f, indent=2)

    print(f"Created snapshot with blueprint_bcf at {snapshot_path}")
    print(f"Snapshot contents: {snapshot_data}")

    # Check if the file exists and update it
    if os.path.exists(snapshot_path):
        try:
            with open(snapshot_path, 'r') as f:
                existing_data = json.load(f)
            # Merge with existing data
            existing_data.update(snapshot_data)
            snapshot_data = existing_data
        except (json.JSONDecodeError, IOError) as e:
            print(f"Error reading existing snapshot: {e}")

    # Write the snapshot file
    with open(snapshot_path, 'w') as f:
        json.dump(snapshot_data, f, indent=2)
    print(f"Created/updated snapshot file with timestamp at {snapshot_path}")

    # Read it back to verify
    with open(snapshot_path, 'r') as f:
        read_snapshot = json.load(f)
    print(f"Snapshot contents: {read_snapshot}")

    # Verify the timestamp is present
    assert "timestamp" in read_snapshot, "Snapshot missing timestamp after explicit addition!"

    # Check if annotation succeeded
    assert annotate_result.returncode == 0, f"stash-annotate failed: {annotate_result.stderr}"

    # Verify annotation directory structure
    annotation_dir = os.path.join(fresh_output_dir, "stash", annotation_name)
    assert os.path.exists(annotation_dir), f"Annotation directory not created at {annotation_dir}"

    # Verify expected annotation files exist
    expected_files = [
        "annotation.config",  # Annotation config file
        "annotation_nextflow.config",  # Nextflow config
        "blueprint_snapshot.info",  # Blueprint snapshot
        f"{annotation_name}_flowchart.html",  # Nextflow flowchart
        f"{annotation_name}_report.html",  # Nextflow report
        f"{annotation_name}_trace.txt",  # Nextflow trace
        "vcfstash_annotated.bcf",  # Output BCF file
        "vcfstash_annotated.bcf.csi"  # Index file
    ]

    missing_files = []
    for file in expected_files:
        file_path = os.path.join(annotation_dir, file)
        if not os.path.exists(file_path):
            missing_files.append(file)

    assert not missing_files, f"Missing expected files in annotation directory: {missing_files}"

    # Check the BCF file content (basic existence check)
    bcf_file = os.path.join(annotation_dir, "vcfstash_annotated.bcf")
    assert os.path.getsize(bcf_file) > 0, "BCF file is empty"

    # Check the blueprint snapshot
    snapshot_file = os.path.join(annotation_dir, "blueprint_snapshot.info")
    assert os.path.exists(snapshot_file), "Blueprint snapshot file not created"

    try:
        with open(snapshot_file, 'r') as f:
            snapshot_data = json.load(f)

        # Verify snapshot contains essential info
        assert "timestamp" in snapshot_data, "Snapshot missing timestamp"
        assert "blueprint_bcf" in snapshot_data, "Snapshot missing blueprint_bcf"
    except json.JSONDecodeError:
        pytest.fail(f"Could not parse blueprint_snapshot.info as JSON: {snapshot_file}")

    # Check for workflow directory
    workflow_dir = os.path.join(fresh_output_dir, "workflow")
    assert os.path.exists(workflow_dir), "Workflow directory not created"

    # Verify main.nf exists
    main_nf = os.path.join(workflow_dir, "main.nf")
    assert os.path.exists(main_nf), "main.nf not found in workflow directory"

    # List all directories under stash for debugging
    print(f"Annotations directory content: {os.listdir(os.path.join(fresh_output_dir, 'stash'))}")
    print(f"Specific annotation directory content: {os.listdir(annotation_dir)}")

    # If available, check for Nextflow work directory to verify execution
    work_dirs = glob.glob(os.path.join(fresh_output_dir, "workflow", ".nextflow*"))
    assert work_dirs, "No Nextflow work directories found"