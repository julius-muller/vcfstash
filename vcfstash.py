#!/usr/bin/env python3
"""
VCFStash - A VCF Annotation Cache Manager

Author: Julius Müller, PhD
Organization: GHGA - German Human Genome-Phenome Archive
Date: 16-03-2025
"""

import os
import sys
from pathlib import Path
from src.cli import main

if __name__ == "__main__":
    # Get the absolute path to the package root directory
    package_root = Path(sys.argv[0]).parent.resolve()

    # Set VCFSTASH_ROOT environment variable
    os.environ['VCFSTASH_ROOT'] = str(package_root)

    # Add package root to Python path if needed
    if package_root != Path.cwd():
        sys.path.insert(0, str(package_root))

    main()
