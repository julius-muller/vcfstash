import pytest
import subprocess
import shutil
import tarfile
from pathlib import Path
# Import required functions from cache.py
from cache import ensure_indexed
import os


# Test data paths - define relative to test file location
TEST_ROOT = Path(__file__).parent / "data"
NODATA_DIR = TEST_ROOT / "nodata"  # Path to nodata directory
TEST_DATA = TEST_ROOT / "test_data"
TEST_DB = TEST_ROOT / "test_db"
TEST_OUT = TEST_ROOT / "test_out"
WORKFLOW_DIR = TEST_ROOT / "workflow"
TEST_PARAMS = TEST_ROOT / "test_params.yml"
REF_FASTA = NODATA_DIR / "reference.fasta"



# Expected MD5 hashes for test cases
EXPECTED_HASHES = {
    "init": "4015425fc735e3ebae192d77334c7ee9",
    "add": "663863e6157bcc14cc6a893e89538e54",
    "annotated": "HASH_TO_BE_FILLED"
}

@pytest.fixture(scope="function")
def setup_test_env():
    """Set up test environment and clean up afterwards"""

    # Clean up existing test directories before test
    shutil.rmtree(TEST_OUT, ignore_errors=True)
    shutil.rmtree(WORKFLOW_DIR, ignore_errors=True)
    if TEST_PARAMS.exists():
        TEST_PARAMS.unlink()

    # Create complete output directory structure including blueprint
    TEST_OUT.mkdir(parents=True, exist_ok=True)
    (TEST_OUT / "test_db" / "blueprint").mkdir(parents=True)
    (TEST_OUT / "test_db" / "annotations").mkdir(parents=True)

    # Create workflow files
    WORKFLOW_DIR.mkdir(parents=True)
    workflow_content = '''
    workflow {
        input = Channel.fromPath(params.input)
        output = params.output
        db_mode = params.db_mode

        ANNOTATE(input, output, db_mode)
    }

    process ANNOTATE {
        input:
            path vcf
            val outdir
            val db_mode
        output:
            path "tool_version.log"

        publishDir "${outdir}", mode: 'copy'

        script:
            """
            echo "test_tool v1.0" > tool_version.log
            """
    }
    '''
    (WORKFLOW_DIR / "main.nf").write_text(workflow_content)

    config_content = '''
    params {
        input = null
        output = null
        db_mode = false
    }

    process {
        publishDir = null
    }
    '''
    (WORKFLOW_DIR / "nextflow.config").write_text(config_content)

    # Create test params file
    TEST_PARAMS.parent.mkdir(parents=True, exist_ok=True)
    TEST_PARAMS.write_text("custom_param: test_value")

    return

def compute_md5(file_path):
    """Compute MD5 hash of a BCF file, filtering out Date headers"""
    cleaned_output = Path(str(file_path) + ".clean")
    try:
        # Filter out Date lines from bcftools headers
        cmd = f"bcftools view {file_path} | grep -v 'Date=' > {cleaned_output}"
        subprocess.run(cmd, shell=True, check=True)
        result = subprocess.run(["md5sum", cleaned_output], check=True, capture_output=True, text=True)
        return result.stdout.split()[0]
    finally:
        cleaned_output.unlink(missing_ok=True)

TEST_CASES = [
    ("init", [
        "python3", "src/cache.py", "init",
        "-n", "test_db",
        "-o", str(TEST_OUT),
        "-i", str(NODATA_DIR / "gnomad_test.bcf"),
        "-f", str(REF_FASTA),
        "-t", "2"
    ]),
    ("add", [
        "python3", "src/cache.py", "add",
        "-d", str(TEST_OUT / "test_db"),
        "-i", str(NODATA_DIR / "dbsnp_test.bcf"),
        "-f", str(REF_FASTA),
        "-t", "2"
    ])
]

def test_file_processing(setup_test_env):
    """Test file indexing"""
    test_files = [
        NODATA_DIR / "gnomad_test.bcf",
        NODATA_DIR / "dbsnp_test.bcf"
    ]
    for file_path in test_files:
        ensure_indexed(file_path)

# Update test_cache_operations to use NODATA_DIR files
@pytest.mark.parametrize("test_case", TEST_CASES)
def test_cache_operations(setup_test_env, test_case):
    """Test cache operations with result validation"""
    case_name, cmd = test_case

    # For add operations, initialize the database first
    if case_name == "add":
        init_cmd = [
            "python3", "src/cache.py", "init",
            "-n", "test_db",
            "-o", str(TEST_OUT),
            "-i", str(NODATA_DIR / "gnomad_test.bcf"),
            "-f", str(REF_FASTA),
            "-t", "2"
        ]
        init_result = subprocess.run(init_cmd, capture_output=True, text=True)
        assert init_result.returncode == 0, f"Database initialization failed: {init_result.stderr}"

    # Run test command
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Command failed: {result.stderr}"

    # Verify file content
    actual_hash = compute_md5(TEST_OUT / "test_db" / "blueprint" / "variants.bcf")
    if EXPECTED_HASHES[case_name] == "HASH_TO_BE_FILLED":
        print(f"\nMD5 hash for {case_name}: {actual_hash}")
    else:
        assert actual_hash == EXPECTED_HASHES[case_name], \
            f"MD5 mismatch for {case_name}: expected {EXPECTED_HASHES[case_name]}, got {actual_hash}"

def test_annotation(setup_test_env):
    """Test annotation workflow on the database"""
    # Initialize database first
    init_cmd = [
        "python3", "src/cache.py", "init",
        "-n", "test_db",
        "-o", str(TEST_OUT),
        "-i", str(NODATA_DIR / "gnomad_test.bcf"),
        "-f", str(REF_FASTA),
        "-t", "2"
    ]
    result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert result.returncode == 0, "Database initialization failed"

    db_path = TEST_OUT / "test_db"
    annotations_dir = db_path / "annotations"

    # Copy workflow files
    real_workflow = Path(__file__).parent.parent / "workflow"
    shutil.copytree(real_workflow, WORKFLOW_DIR, dirs_exist_ok=True)

    # Run annotation with Nextflow args
    annotate_cmd = [
        "python3", "src/cache.py", "annotate",
        "-d", str(db_path),
        "-w", str(WORKFLOW_DIR),
        "-p", str(TEST_PARAMS),
        "--",  # separator for Nextflow args
        "-profile", "test",
        "--max_cpus", "4"
    ]
    result = subprocess.run(annotate_cmd, capture_output=True, text=True)

    # Print full output for debugging
    print("\nSTDOUT:")
    print(result.stdout)
    print("\nSTDERR:")
    print(result.stderr)

    assert result.returncode == 0, f"Annotation failed. Exit code: {result.returncode}\nOutput: {result.stdout}\nError: {result.stderr}"

    # Get the archive
    archives = sorted(annotations_dir.glob("*.vepstash"))
    assert archives, "No annotation archive found"
    archive = archives[0]

    # Extract the BCF file and compute its hash
    with tarfile.open(archive, "r:gz") as tar:
        print("\nArchive contents:")
        for member in tar.getmembers():
            print(f"  {member.name}")
        bcf_member = next(m for m in tar.getmembers() if m.name.endswith(".bcf"))
        tar.extract(bcf_member, path=annotations_dir)

        bcf_file = annotations_dir / bcf_member.name
        actual_hash = compute_md5(bcf_file)
        print(f"\nAnnotated BCF MD5: {actual_hash}")

        if EXPECTED_HASHES["annotated"] == "HASH_TO_BE_FILLED":
            print(f"Fill in EXPECTED_HASHES['annotated'] with: {actual_hash}")
        else:
            assert actual_hash == EXPECTED_HASHES["annotated"], \
                f"Annotation MD5 mismatch: expected {EXPECTED_HASHES['annotated']}, got {actual_hash}"

        # Clean up
        bcf_file.unlink()


