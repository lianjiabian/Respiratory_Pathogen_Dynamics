"""Validate RDA artifacts with Python/pyreadr and, when supplied, R itself."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import pyreadr


R_CODE = r'''
args <- commandArgs(trailingOnly=TRUE)
e <- new.env(parent=emptyenv())
objects <- load(args[[1]], envir=e)
expected <- args[[2]]
if (!identical(objects, expected)) stop(paste("unexpected objects", paste(objects, collapse=",")))
x <- e[[expected]]
if (!is.data.frame(x)) stop("object is not a data.frame")
cat(sprintf("RDA_OK|%s|%s|%d|%d\n", expected, paste(class(x), collapse="/"), nrow(x), ncol(x)))
'''


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rscript", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = pd.read_csv(args.manifest.resolve(), encoding="utf-8-sig")
    rscript = args.rscript.resolve() if args.rscript else None
    rows: list[dict[str, object]] = []
    for i, item in enumerate(manifest.itertuples(index=False), start=1):
        path = root / item.rda_file
        expected_object = str(item.object_name)
        expected_rows = int(item.rows)
        expected_columns = int(item.columns)
        loaded = pyreadr.read_r(str(path))
        python_ok = list(loaded) == [expected_object]
        if python_ok:
            frame = loaded[expected_object]
            python_ok = frame.shape == (expected_rows, expected_columns)
        r_status = "NOT_RUN"
        r_detail = "Rscript not supplied"
        if rscript:
            with tempfile.TemporaryDirectory(prefix="rda_validation_") as temp_dir:
                temp_path = Path(temp_dir) / f"artifact_{i:02d}.rda"
                shutil.copy2(path, temp_path)
                result = subprocess.run(
                    [str(rscript), "-e", R_CODE, str(temp_path), expected_object],
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=False,
                )
                marker = next((line for line in result.stdout.splitlines() if line.startswith("RDA_OK|")), "")
                expected_marker = f"RDA_OK|{expected_object}|data.frame|{expected_rows}|{expected_columns}"
                r_status = "PASS" if result.returncode == 0 and marker == expected_marker else "FAIL"
                r_detail = marker or (result.stderr.strip()[-1000:] if result.stderr else "no R marker")
        rows.append({
            "object_name": expected_object,
            "rda_file": item.rda_file,
            "expected_rows": expected_rows,
            "expected_columns": expected_columns,
            "python_pyreadr_status": "PASS" if python_ok else "FAIL",
            "r_load_status": r_status,
            "r_detail": r_detail,
        })
    report = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.output, index=False, encoding="utf-8-sig", lineterminator="\n")
    if not report["python_pyreadr_status"].eq("PASS").all():
        raise SystemExit("Python RDA validation failed")
    if rscript and not report["r_load_status"].eq("PASS").all():
        raise SystemExit("R load validation failed")
    print(report.to_string(index=False))


if __name__ == "__main__":
    main()
