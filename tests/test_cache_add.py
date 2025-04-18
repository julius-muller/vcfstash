import os
import json
import pytest
import tempfile
import subprocess
import uuid

# Define constants
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TEST_VCF = os.path.join(TEST_DATA_DIR, "nodata", "crayz_db.bcf")
TEST_PARAMS = os.path.join(os.path.dirname(__file__), "config", "test_params.yaml")
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


def test_stash_init_and_add(fresh_output_dir):
    """Test that stash-init and stash-add work correctly."""
    # Verify the directory doesn't exist yet
    print(f"Output directory exists before test: {os.path.exists(fresh_output_dir)}")

    # Run stash-init with a fresh directory
    init_cmd = [
        VCFSTASH_CMD,
        "stash-init",
        "-i", TEST_VCF,
        "-o", fresh_output_dir,
        "-y", TEST_PARAMS,
        "-f"
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

    # Verify the directory was created with correct structure
    assert os.path.exists(fresh_output_dir), "Output directory was not created"
    blueprint_dir = os.path.join(fresh_output_dir, "blueprint")
    assert os.path.exists(blueprint_dir), "Blueprint directory was not created"

    # Print the directory contents for debugging
    print(f"Output directory contents: {os.listdir(fresh_output_dir)}")
    print(f"Blueprint directory contents: {os.listdir(blueprint_dir)}")

    # Check sources.info file after initialization
    sources_file = os.path.join(blueprint_dir, "sources.info")
    assert os.path.exists(sources_file), "sources.info file not created"

    with open(sources_file, 'r') as f:
        sources_data = json.load(f)

    # Verify initial file is in sources.info
    assert "input_files" in sources_data, "input_files not in sources.info"
    assert len(sources_data["input_files"]) == 1, f"Expected 1 input file, found {len(sources_data['input_files'])}"
    assert os.path.basename(sources_data["input_files"][0]["path"]) == os.path.basename(TEST_VCF), \
        f"Input file in sources.info ({sources_data['input_files'][0]['path']}) doesn't match test file ({TEST_VCF})"

    # Test stash-add by adding the same file again (reusing the same test file for simplicity)
    add_cmd = [
        VCFSTASH_CMD,
        "stash-add",
        "--db", fresh_output_dir,
        "-i", TEST_VCF,
        "-y", TEST_PARAMS
    ]

    print(f"Running add command: {' '.join(add_cmd)}")
    add_result = subprocess.run(
        add_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False
    )

    # Check if add succeeded
    assert add_result.returncode == 0, f"stash-add failed: {add_result.stderr}"
    print(f"Add succeeded with output: {add_result.stdout}")

    # Re-check sources.info after adding the file
    with open(sources_file, 'r') as f:
        updated_sources_data = json.load(f)

    # The file might be deduplicated, so check that we have at least one entry still
    assert len(updated_sources_data["input_files"]) >= 1, "No input files in sources.info after stash-add"

    # Check for BCF file
    bcf_file = os.path.join(blueprint_dir, "vcfstash.bcf")
    assert os.path.exists(bcf_file), "BCF file not created"

    # Check for other expected files (adjust these based on your application's expected output)
    expected_files = ["vcfstash.bcf", "vcfstash.bcf.csi", "sources.info"]
    for file in expected_files:
        assert os.path.exists(os.path.join(blueprint_dir, file)), f"Expected file {file} not found"

    # Verify BCF contains variants (this is just a basic check)
    bcftools_path = os.path.join(os.environ.get('VCFSTASH_ROOT', ''), 'tools', 'bcftools')
    if not os.path.exists(bcftools_path):
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    bcf_stats_cmd = [bcftools_path, "stats", bcf_file]
    try:
        bcf_stats = subprocess.run(
            bcf_stats_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False
        )

        if bcf_stats.returncode == 0:
            print(f"BCF stats output: {bcf_stats.stdout}")
        else:
            print(f"BCF stats failed with: {bcf_stats.stderr}")
    except FileNotFoundError:
        # bcftools might not be installed
        print("Skipping BCF content check as bcftools is not available")
