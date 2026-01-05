# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Allow vcfcache to be run as a module: python -m vcfcache"""

from vcfcache.cli import main

if __name__ == "__main__":
    main()
