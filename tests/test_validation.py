import os
import tempfile
import pytest
from pathlib import Path
from src.utils.validation import compute_md5

def test_compute_md5():
    """Test that compute_md5 correctly calculates MD5 hash of a file."""
    # Create a temporary file with known content
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write(b"test content for MD5 calculation")
        temp_file_path = temp_file.name

    try:
        # Calculate MD5 hash
        calculated_md5 = compute_md5(Path(temp_file_path))

        # Expected MD5 hash for "test content for MD5 calculation"
        # This was calculated from the actual file created by the test
        expected_md5 = "b4e714277ba501ef9b2ed937048dc9a4"

        # Assert that the calculated hash matches the expected hash
        assert calculated_md5 == expected_md5, f"Expected {expected_md5}, got {calculated_md5}"
    finally:
        # Clean up the temporary file
        os.unlink(temp_file_path)
