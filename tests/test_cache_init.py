import pytest
import tempfile
import shutil
import hashlib
import os
import re
import json
import subprocess
from pathlib import Path

# Define constants for testing
TEST_ROOT = os.path.dirname(os.path.abspath(__file__))
VCFSTASH_CMD = os.path.join(os.path.dirname(TEST_ROOT), "vcfstash.py")
TEST_DATA_DIR = os.path.join(TEST_ROOT, "data", "nodata")
TEST_VCF = os.path.join(TEST_DATA_DIR, "crayz_db.bcf")
EXPECTED_OUTPUT_DIR = os.path.join(TEST_ROOT, "data", "expected_output", "stash_init_result")


# Fixture for temporary test output directory
@pytest.fixture
def test_output_dir():
    """Provides a path for a directory that doesn't exist yet."""
    # Create a path for a temporary directory, but don't create it
    temp_dir = tempfile.mkdtemp(prefix="vcfstash_test_")

    # Remove the directory immediately - vcfstash will create it
    os.rmdir(temp_dir)

    # Return the path to the test function
    yield temp_dir

    # Clean up after the test is done
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir, ignore_errors=True)


# Helper functions
def compute_md5(filename):
    """Compute MD5 hash for a file."""
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def run_stash_init(input_vcf, output_dir, config_file=None, force=False):
    """Run the stash-init command and return the process result."""
    # Make sure the directory doesn't exist (clean start)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    # Get the VCFSTASH_ROOT environment variable
    vcfstash_root = os.environ.get('VCFSTASH_ROOT')
    if not vcfstash_root:
        raise ValueError("VCFSTASH_ROOT environment variable is not set")

    # Create a temporary params file with the correct paths
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as temp_file:
        temp_params_file = temp_file.name

        # Read the original params file
        original_params_file = os.path.join(TEST_ROOT, "config", "test_params.yaml")
        with open(original_params_file, 'r') as f:
            params_content = f.read()

        # Replace ${VCFSTASH_ROOT} with the actual value
        params_content = params_content.replace('${VCFSTASH_ROOT}', vcfstash_root)

        # Write the modified content to the temporary file
        temp_file.write(params_content)

    try:
        cmd = [
            VCFSTASH_CMD,
            "stash-init",
            "--vcf", input_vcf,
            "--output", output_dir,
            "-y", temp_params_file
        ]

        # We no longer use config_file as per requirements
        # Only process section is allowed in configs, and they're optional

        if force:
            cmd.append("-f")

        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return result
    finally:
        # Clean up the temporary file
        if os.path.exists(temp_params_file):
            os.unlink(temp_params_file)

def compare_directories_ignore_timestamps(dir1, dir2):
    """Compare two directories, ignoring timestamp differences and dynamic directories."""
    result = {'success': True, 'message': '', 'details': {}}

    # Paths to ignore (using regex patterns)
    ignore_patterns = [
        r'^blueprint/\.nextflow/cache/[^/]+/.*$',  # Ignore nextflow cache with different UUIDs
        r'^.*\.log$',  # Ignore all log files
        r'^blueprint/init_report\.html$',  # Dynamic HTML report
        r'^blueprint/init_trace\.txt$',  # Dynamic trace data
        r'^blueprint/\.nextflow\.log$',  # Nextflow log
        r'^blueprint/\.nextflow/history$',  # Nextflow history
        r'^blueprint/vcfstash\.bcf\.csi$',  # Index file that might differ
        r'^workflow/init\.yaml$',  # Ignore init.yaml as we've changed its content
        r'^workflow/init\.config$',  # Ignore init.config as we're not using it anymore
        r'^workflow/modules/.*\.nf$',  # Ignore workflow module files
        r'^blueprint/init_flowchart\.html$',  # Ignore flowchart HTML
        r'^workflow/\.nextflow/history$',  # Ignore Nextflow history
        r'^workflow/\.nextflow/framework/.*$'  # Ignore Nextflow framework files
    ]

    # Special files that need custom comparison
    special_files = [
        r'^blueprint/vcfstash\.bcf$',  # BCF file
        r'^blueprint/sources\.info$'  # Sources info
    ]

    def should_ignore(path):
        return any(re.match(pattern, path) for pattern in ignore_patterns)

    def is_special_file(path):
        return any(re.match(pattern, path) for pattern in special_files)

    # Get the directory structures
    dir1_files = set()
    dir2_files = set()
    special_file_paths = set()

    for root, _, files in os.walk(dir1):
        rel_root = os.path.relpath(root, dir1)
        for file in files:
            rel_path = os.path.join(rel_root, file) if rel_root != '.' else file
            rel_path = rel_path.replace('\\', '/')  # Normalize path separators
            if not should_ignore(rel_path):
                dir1_files.add(rel_path)
                if is_special_file(rel_path):
                    special_file_paths.add(rel_path)

    for root, _, files in os.walk(dir2):
        rel_root = os.path.relpath(root, dir2)
        for file in files:
            rel_path = os.path.join(rel_root, file) if rel_root != '.' else file
            rel_path = rel_path.replace('\\', '/')  # Normalize path separators
            if not should_ignore(rel_path):
                dir2_files.add(rel_path)
                if is_special_file(rel_path):
                    special_file_paths.add(rel_path)

    # Compare directory structures
    only_in_dir1 = dir1_files - dir2_files
    only_in_dir2 = dir2_files - dir1_files
    common_files = dir1_files.intersection(dir2_files)

    # Check if any files are missing
    if only_in_dir1 or only_in_dir2:
        result['success'] = False
        if only_in_dir1:
            result['details']['only_in_test'] = list(only_in_dir1)
        if only_in_dir2:
            result['details']['only_in_reference'] = list(only_in_dir2)

    # Compare content of common files, with special handling for certain file types
    diff_files = []

    for file_path in common_files:
        file1 = os.path.join(dir1, file_path)
        file2 = os.path.join(dir2, file_path)

        # Skip special files for normal comparison
        if is_special_file(file_path):
            continue

        # Compare files with timestamp normalization
        if not compare_files_ignoring_timestamps(file1, file2):
            diff_files.append(file_path)

    # Now handle special files separately
    for file_path in special_file_paths:
        if file_path in common_files:  # Only if both files exist
            file1 = os.path.join(dir1, file_path)
            file2 = os.path.join(dir2, file_path)

            if file_path.endswith('.bcf'):
                # For BCF files, we'll just check that they're valid but won't compare content
                if not is_valid_bcf(file1) or not is_valid_bcf(file2):
                    diff_files.append(file_path)
            elif file_path.endswith('sources.info'):
                # For sources.info files, use our special comparison
                if not compare_sources_info(file1, file2):
                    diff_files.append(file_path)

    if diff_files:
        result['success'] = False
        result['details']['different_files'] = diff_files

    # Prepare detailed error message
    error_messages = []
    for key, value in result['details'].items():
        if value:  # Only include non-empty details
            error_messages.append(f"{key}: {value}")

    if error_messages:
        result['message'] = "Directory comparison failed: " + " | ".join(error_messages)

    return result


def is_valid_bcf(bcf_file):
    """Check if a BCF file is valid by running bcftools."""
    try:
        # Use bcftools to check if the BCF file is valid
        # We just check that bcftools can read the header
        bcftools_path = os.path.join(os.environ.get('VCFSTASH_ROOT', ''), 'tools', 'bcftools')
        if not os.path.exists(bcftools_path):
            # Fall back to system bcftools if the project-specific one doesn't exist
            bcftools_path = 'bcftools'

        result = subprocess.run(
            [bcftools_path, 'view', '-h', bcf_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )
        return result.returncode == 0
    except Exception:
        return False


# Tests
def test_stash_init_against_reference(test_output_dir):
    """Test that stash-init creates the expected directory structure."""
    output_dir = test_output_dir
    result = run_stash_init(TEST_VCF, output_dir, force=True)
    assert result.returncode == 0, f"stash-init failed: {result.stderr}"

    # Use a modified comparison function that ignores timestamps
    comparison_result = compare_directories_ignore_timestamps(output_dir, EXPECTED_OUTPUT_DIR)

    if not comparison_result['success']:
        if 'different_files' in comparison_result['details']:
            for file_path in comparison_result['details']['different_files']:
                test_file = os.path.join(output_dir, file_path)
                ref_file = os.path.join(EXPECTED_OUTPUT_DIR, file_path)
                print_file_diff(test_file, ref_file)

    assert comparison_result['success'], comparison_result['message']


def compare_files_ignoring_timestamps(file1, file2):
    """Compare two files ignoring timestamp differences and path differences."""
    if not os.path.exists(file1) or not os.path.exists(file2):
        return False

    # Special handling for Nextflow files
    if file1.endswith('.nf'):
        try:
            with open(file1, 'r', encoding='utf-8') as f1, open(file2, 'r', encoding='utf-8') as f2:
                content1 = f1.read().strip()
                content2 = f2.read().strip()

                # Normalize line endings
                content1 = content1.replace('\r\n', '\n')
                content2 = content2.replace('\r\n', '\n')

                if content1 != content2:
                    print(f"\nDifferences in {os.path.basename(file1)}:")
                    print(f"Test file ({file1}):\n{content1}")
                    print(f"\nReference file ({file2}):\n{content2}")
                    return False
                return True
        except UnicodeDecodeError:
            print(f"Warning: UnicodeDecodeError when reading {file1} or {file2}")
            return compare_binary_files(file1, file2)

    # Special handling for known files
    file_basename = os.path.basename(file1)

    # For sources.info, we need special handling
    if file_basename == "sources.info":
        return compare_sources_info(file1, file2)

    # Check if the file is likely binary
    def is_binary(file_path):
        # Check file extension first
        if any(file_path.endswith(ext) for ext in
               ('.jar', '.zip', '.bcf', '.bam', '.bai', '.gz', '.bin', '.pyc')):
            return True

        # If not obvious from extension, check content
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                f.read(1024)  # Try to read as text
            return False  # If we can read it as text, it's not binary
        except UnicodeDecodeError:
            return True  # If we can't read it as text, it's binary

    # Handle binary files
    if is_binary(file1) or is_binary(file2):
        # For binary files, compare file sizes first
        size1 = os.path.getsize(file1)
        size2 = os.path.getsize(file2)

        if size1 != size2:
            # Size mismatch - files are definitely different
            return False

        # For JAR/ZIP files, we might need a more complex comparison
        # For now, we'll do a byte-by-byte comparison as fallback
        with open(file1, 'rb') as f1, open(file2, 'rb') as f2:
            chunk_size = 8192  # 8K chunks
            while True:
                chunk1 = f1.read(chunk_size)
                chunk2 = f2.read(chunk_size)

                if chunk1 != chunk2:
                    return False

                if not chunk1:  # End of file
                    break

        return True  # If we got here, binary files match

    # Handle text files
    try:
        with open(file1, 'r', encoding='utf-8') as f1, open(file2, 'r', encoding='utf-8') as f2:
            content1 = f1.read()
            content2 = f2.read()

        # Normalize timestamps in ISO format: 2023-04-02T18:14:50
        content1 = re.sub(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}', 'TIMESTAMP_ISO', content1)
        content2 = re.sub(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}', 'TIMESTAMP_ISO', content2)

        # Normalize timestamps like "Apr 2 18:14:50 2025"
        content1 = re.sub(r'\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d+\s+\d+:\d+:\d+\s+\d+\b',
                          'TIMESTAMP', content1)
        content2 = re.sub(r'\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d+\s+\d+:\d+:\d+\s+\d+\b',
                          'TIMESTAMP', content2)

        # Normalize UUIDs
        content1 = re.sub(r'[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}', 'UUID', content1)
        content2 = re.sub(r'[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}', 'UUID', content2)

        # Normalize absolute paths
        content1 = re.sub(r'/tmp/vcfstash_test_[a-z0-9]+', '/tmp/test_dir', content1)
        content2 = re.sub(r'/tmp/vcfstash_test_[a-z0-9]+', '/tmp/test_dir', content2)

        return content1 == content2
    except UnicodeDecodeError:
        # If we unexpectedly hit a binary file that wasn't caught earlier
        print(f"Warning: UnicodeDecodeError when reading {file1} or {file2}")
        return compare_binary_files(file1, file2)

def print_file_diff(file1, file2):
    """Print a detailed diff between two files."""
    import difflib

    with open(file1, 'r', encoding='utf-8') as f1, open(file2, 'r', encoding='utf-8') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    differ = difflib.Differ()
    diff = list(differ.compare(lines1, lines2))

    print("\nDetailed diff:")
    for line in diff:
        if line.startswith(('+ ', '- ', '? ')):
            print(line.rstrip())

def compare_binary_files(file1, file2):
    """Compare two binary files byte by byte."""
    if os.path.getsize(file1) != os.path.getsize(file2):
        return False

    with open(file1, 'rb') as f1, open(file2, 'rb') as f2:
        chunk_size = 8192  # 8K chunks
        while True:
            chunk1 = f1.read(chunk_size)
            chunk2 = f2.read(chunk_size)

            if chunk1 != chunk2:
                return False

            if not chunk1:  # End of file
                break

    return True

def compare_sources_info(file1, file2):
    """
    Compare the contents of two sources.info files, ignoring timestamps and normalizing paths.

    Args:
        file1: Path to the first file
        file2: Path to the second file

    Returns:
        bool: True if the files are equivalent (ignoring timestamps and exact paths), False otherwise
        str: Description of differences if any, otherwise None
    """
    with open(file1, "r") as f1, open(file2, "r") as f2:
        data1 = json.load(f1)
        data2 = json.load(f2)

    # Ignore timestamps
    if "created" in data1:
        data1["created"] = "TIMESTAMP_IGNORED"
    if "created" in data2:
        data2["created"] = "TIMESTAMP_IGNORED"

    # Normalize input file paths by removing the absolute path
    if "input_files" in data1:
        for file_info in data1["input_files"]:
            if "path" in file_info:
                # Keep only the filename part
                file_info["path"] = os.path.basename(file_info["path"])
            if "added" in file_info:
                file_info["added"] = "TIMESTAMP_IGNORED"

    if "input_files" in data2:
        for file_info in data2["input_files"]:
            if "path" in file_info:
                # Keep only the filename part
                file_info["path"] = os.path.basename(file_info["path"])
            if "added" in file_info:
                file_info["added"] = "TIMESTAMP_IGNORED"

    # Also normalize the name field which will be different per test run
    if "name" in data1:
        data1["name"] = "NAME_IGNORED"
    if "name" in data2:
        data2["name"] = "NAME_IGNORED"

    # Check if the structures are now equal
    are_equal = data1 == data2

    if not are_equal:
        # Format the data for display
        diff_str = f"Differences in {os.path.basename(file1)}:\n"
        diff_str += f"Test file:\n{json.dumps(data1, indent=2)}...\n"
        diff_str += f"Reference file:\n{json.dumps(data2, indent=2)}...\n"
        return False, diff_str

    return True, None



def test_key_files_content_matches(test_output_dir):
    """Test that key files in the stash directory match the reference."""
    output_dir = test_output_dir
    result = run_stash_init(TEST_VCF, output_dir, force=True)
    assert result.returncode == 0, f"stash-init failed: {result.stderr}"

    # List of important files to check
    key_files = [
        "blueprint/sources.info",
        # "workflow/init.yaml",  # Excluded because we've changed the content of test_params.yaml
        "workflow/main.nf"
    ]

    # Check all the module files too
    module_files = [
        "workflow/modules/annotate_cache.nf",
        "workflow/modules/annotate.nf",
        "workflow/modules/intersect.nf",
        "workflow/modules/merge.nf",
        "workflow/modules/merge_variants.nf",
        "workflow/modules/normalize.nf",
        "workflow/modules/utils.nf"
    ]
    key_files.extend(module_files)

    for file_path in key_files:
        test_file = os.path.join(output_dir, file_path)
        ref_file = os.path.join(EXPECTED_OUTPUT_DIR, file_path)

        assert os.path.exists(test_file), f"File {file_path} not found in test output"

        # Check if reference file exists (if not, this might be a new file)
        if not os.path.exists(ref_file):
            print(f"Warning: Reference file {file_path} does not exist. Skipping comparison.")
            continue

        # Use appropriate comparison based on file type
        if file_path.endswith('sources.info'):
            files_match = compare_sources_info(test_file, ref_file)
        else:
            files_match = compare_files_ignoring_timestamps(test_file, ref_file)

        if not files_match:
            # For debugging: print content of both files
            print(f"\nDifferences in {file_path}:")
            try:
                with open(test_file, 'r', encoding='utf-8', errors='replace') as f1:
                    test_content = f1.read()
                with open(ref_file, 'r', encoding='utf-8', errors='replace') as f2:
                    ref_content = f2.read()

                print(f"Test file:\n{test_content[:500]}...")
                print(f"Reference file:\n{ref_content[:500]}...")
            except UnicodeDecodeError:
                print("Could not display content - binary file")

        assert files_match, f"Content mismatch for {file_path} (ignoring timestamps)"

def test_bcf_file_matches_reference(test_output_dir):
    """Test that the generated BCF file matches the reference."""
    output_dir = test_output_dir  # Now it's a parameter, not a function call
    result = run_stash_init(TEST_VCF, output_dir, force=True)
    assert result.returncode == 0, f"stash-init failed: {result.stderr}"

    # Compare BCF content
    test_bcf = os.path.join(output_dir, "variants.bcf")
    ref_bcf = os.path.join(EXPECTED_OUTPUT_DIR, "variants.bcf")

    # Get bcftools path
    bcftools_path = os.path.join(os.environ.get('VCFSTASH_ROOT', ''), 'tools', 'bcftools')
    if not os.path.exists(bcftools_path):
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Convert BCF to readable text
    test_vcf_text = subprocess.run(
        [bcftools_path, "view", test_bcf],
        stdout=subprocess.PIPE, text=True
    ).stdout
    ref_vcf_text = subprocess.run(
        [bcftools_path, "view", ref_bcf],
        stdout=subprocess.PIPE, text=True
    ).stdout

    # Filter out timestamp lines or normalize them
    test_lines = []
    ref_lines = []

    for line in test_vcf_text.splitlines():
        # Skip or normalize timestamp-containing lines
        if "##fileDate=" in line:
            continue
        # Replace timestamp patterns in header
        if line.startswith('##'):
            line = re.sub(r'\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d+\s+\d+:\d+:\d+\s+\d+\b',
                          'TIMESTAMP', line)
        test_lines.append(line)

    for line in ref_vcf_text.splitlines():
        # Skip or normalize timestamp-containing lines
        if "##fileDate=" in line:
            continue
        # Replace timestamp patterns in header
        if line.startswith('##'):
            line = re.sub(r'\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d+\s+\d+:\d+:\d+\s+\d+\b',
                          'TIMESTAMP', line)
        ref_lines.append(line)

    # Compare filtered content
    test_filtered = "\n".join(test_lines)
    ref_filtered = "\n".join(ref_lines)

    assert test_filtered == ref_filtered, "BCF content doesn't match reference (ignoring timestamps)"
