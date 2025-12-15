from __future__ import annotations

import argparse
import re
from pathlib import Path


def extract_section(changelog: str, version: str) -> str:
    header_re = re.compile(r"^##\s+(.+?)\s*$", re.M)
    headers = list(header_re.finditer(changelog))

    target_idx = None
    for i, m in enumerate(headers):
        if m.group(1).startswith(version):
            target_idx = i
            break
        if m.group(1).startswith(f"[{version}]"):
            target_idx = i
            break

    if target_idx is None:
        raise SystemExit(f"Could not find a '## {version}' section in CHANGELOG.md")

    start = headers[target_idx].start()
    end = headers[target_idx + 1].start() if target_idx + 1 < len(headers) else len(changelog)
    section = changelog[start:end].strip()
    return section + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract a version section from CHANGELOG.md for GitHub Release notes."
    )
    parser.add_argument("version", help="Version to extract (e.g. 0.4.0b1)")
    parser.add_argument(
        "--changelog",
        default="CHANGELOG.md",
        help="Path to CHANGELOG.md (default: CHANGELOG.md)",
    )
    args = parser.parse_args()

    changelog_path = Path(args.changelog)
    text = changelog_path.read_text(encoding="utf-8")
    print(extract_section(text, args.version))


if __name__ == "__main__":
    main()

