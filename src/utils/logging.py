import logging
import sys
from pathlib import Path
from typing import Optional


def setup_logger(
    log_file: Path,
    level: int = logging.INFO,
    name: str = "vepstash",
    format_string: Optional[str] = None,
) -> logging.Logger:
    """
    Configure and return a logger with both file and console handlers.

    Args:
        log_file: Path to the log file
        level: Logging level (default: INFO)
        name: Logger name (default: vepstash)
        format_string: Custom format string for log messages

    Returns:
        logging.Logger: Configured logger instance
    """
    if format_string is None:
        format_string = "[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s"

    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Prevent adding handlers multiple times
    if logger.hasHandlers():
        logger.handlers.clear()

    # Create formatters and handlers
    formatter = logging.Formatter(format_string, datefmt="%Y-%m-%d %H:%M:%S")

    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


def log_command(logger: logging.Logger) -> None:
    """
    Log the command used to execute the script.

    Args:
        logger: Logger instance to use
    """
    command = f"python3 {sys.argv[0]} {' '.join(sys.argv[1:])}"
    logger.info("Script command: %s", command)


# Example usage
if __name__ == "__main__":
    logger = setup_logger(
        log_file=Path("vepstash.log"),
        level=logging.DEBUG,
    )
    logger.debug("Debug message")
    logger.info("Info message")
    logger.warning("Warning message")
    logger.error("Error message")
    log_command(logger)