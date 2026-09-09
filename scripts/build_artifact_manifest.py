"""Create a SHA-256 inventory of all files under a release root except itself."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    rows = []
    for path in sorted(item for item in root.rglob("*") if item.is_file() and item.resolve() != output):
        if any(part in {"__pycache__", ".git"} for part in path.relative_to(root).parts) or path.suffix == ".pyc":
            continue
        rows.append({"file": path.relative_to(root).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)})
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False, encoding="utf-8-sig", lineterminator="\n")
    print(f"Wrote {len(rows)} entries to {output}")


if __name__ == "__main__":
    main()
