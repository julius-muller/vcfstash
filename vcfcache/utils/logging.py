# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Logging utilities for the vcfcache package.

This module provides functions for setting up logging with different verbosity levels
and for logging command execution.
"""

import logging
import sys
from pathlib import Path
from typing import Optional


def setup_logging(verbosity: int, log_file: Optional[Path] = None) -> logging.Logger:
    """Setup logging with verbosity levels.

    Args:
        verbosity: Verbosity level (0-2)
        log_file: Optional path to log file. If None, file logging is disabled

    Returns:
        Configured logger
    """
    # Console is one level less verbose than file
    console_levels = {
        0: logging.WARNING,  # -q (quiet)
        1: logging.INFO,  # default
        2: logging.DEBUG,  # -v
    }
    file_levels = {
        0: logging.INFO,  # -q (quiet)
        1: logging.DEBUG,  # default
        2: logging.DEBUG,  # -v
    }

    console_level = console_levels.get(min(verbosity, 2), logging.WARNING)
    file_level = file_levels.get(min(verbosity, 2), logging.INFO)

    logger = logging.getLogger("vcfcache")
    logger.setLevel(logging.DEBUG)  # Allow all levels to handlers

    # Create formatters - different for console and file
    # Console: Clean format, no timestamps (unless in debug mode)
    if verbosity >= 2:
        # Debug mode: show detailed info
        console_formatter = logging.Formatter(
            "[%(asctime)s] %(levelname)s [%(filename)s.%(funcName)s:%(lineno)d] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    elif verbosity == 1:
        # Info mode: show timestamp but cleaner
        console_formatter = logging.Formatter(
            "[%(asctime)s] %(message)s",
            datefmt="%H:%M:%S",
        )
    else:
        # Default: minimal, clean output
        console_formatter = logging.Formatter("%(message)s")

    # File: Always detailed
    file_formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s [%(name)s.%(filename)s.%(funcName)s:%(lineno)d] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Prevent adding handlers multiple times
    if logger.hasHandlers():
        logger.handlers.clear()

    # Console handler - less verbose
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(console_level)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # File handler - more verbose (only if log_file is provided)
    if log_file is not None:
        # Create parent directory if it doesn't exist
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, mode="a")
        file_handler.setLevel(file_level)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
        logger.debug(f"File logging enabled at: {log_file}")
        if not log_file.exists():
            raise FileNotFoundError(
                f"File logging enabled but log file does not exist: {log_file}"
            )

    logger.debug(
        f"Log levels - Console: {logging.getLevelName(console_level)}"
        + (f", File: {logging.getLevelName(file_level)}" if log_file else "")
    )
    return logger


def log_command(logger: logging.Logger, info: bool = False) -> None:
    """Log the command used to execute the script.

    Args:
        logger: Logger instance to use for logging
        info: If True, log at INFO level; otherwise, log at DEBUG level
    """
    command = f"python3 {sys.argv[0]} {' '.join(sys.argv[1:])}"
    if info:
        logger.info("Script command: %s", command)
    else:
        logger.debug("Script command: %s", command)
