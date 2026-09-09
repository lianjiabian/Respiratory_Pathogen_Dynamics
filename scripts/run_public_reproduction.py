"""Validate the 11 public tables and recreate compact reported summaries."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from importlib.metadata import version
from pathlib import Path

import pandas as pd
import pyreadr


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().eq("true")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(root / "release_manifest.csv", encoding="utf-8-sig")
    expected_tables = [f"S{i}" for i in range(1, 12)]
    if manifest["table"].tolist() != expected_tables:
        raise ValueError("release_manifest.csv must contain S1-S11 in order")
    privacy = pd.read_csv(root / "qa" / "privacy_scan.csv", encoding="utf-8-sig")
    if privacy["table"].tolist() != expected_tables or not privacy["status"].eq("PASS").all():
        raise ValueError("All 11 public tables must have privacy status PASS")

    tables: dict[str, pd.DataFrame] = {}
    checks: list[dict[str, object]] = []
    for item in manifest.itertuples(index=False):
        csv_path = root / item.csv_file
        rda_path = root / item.rda_file
        csv_frame = pd.read_csv(csv_path, encoding="utf-8-sig")
        objects = pyreadr.read_r(str(rda_path))
        object_names = list(objects)
        rda_frame = objects[str(item.object_name)] if object_names == [str(item.object_name)] else None
        values_match = False
        if rda_frame is not None:
            try:
                pd.testing.assert_frame_equal(
                    csv_frame.reset_index(drop=True),
                    rda_frame.reset_index(drop=True),
                    check_dtype=False,
                    check_exact=False,
                    rtol=1e-12,
                    atol=1e-12,
                )
                values_match = True
            except AssertionError:
                values_match = False
        passed = (
            sha256(csv_path) == str(item.csv_sha256)
            and sha256(rda_path) == str(item.rda_sha256)
            and csv_frame.shape == (int(item.rows), int(item.columns))
            and rda_frame is not None
            and rda_frame.shape == (int(item.rows), int(item.columns))
            and list(csv_frame.columns) == list(rda_frame.columns)
            and values_match
        )
        checks.append(
            {
                "table": item.table,
                "object_name": item.object_name,
                "rows": int(item.rows),
                "columns": int(item.columns),
                "csv_sha256_match": sha256(csv_path) == str(item.csv_sha256),
                "rda_sha256_match": sha256(rda_path) == str(item.rda_sha256),
                "python_rda_load": object_names == [str(item.object_name)],
                "shape_and_columns_match": (
                    rda_frame is not None
                    and csv_frame.shape == rda_frame.shape
                    and list(csv_frame.columns) == list(rda_frame.columns)
                ),
                "csv_rda_values_match": values_match,
                "status": "PASS" if passed else "FAIL",
            }
        )
        if not passed:
            raise ValueError(f"Public table validation failed: {item.table}")
        tables[str(item.table)] = csv_frame

    check_frame = pd.DataFrame(checks)
    check_frame.to_csv(
        output_dir / "public_table_validation.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )

    cohort = tables["S1"].sort_values(["phase", "pathogen"])
    demographics = tables["S2"].sort_values("phase")
    age = tables["S6"].copy()
    age = age.loc[~as_bool(age["suppressed"])].copy()
    age["rate_pct"] = pd.to_numeric(age["rate_pct"], errors="coerce")
    age = age.loc[age["rate_pct"].notna()].copy()
    age_peaks = age.loc[
        age.groupby(["phase", "pathogen"], observed=True)["rate_pct"].idxmax()
    ].sort_values(["phase", "pathogen"])
    codetection_age = tables["S7"].copy()
    pairs = tables["S8"].loc[~as_bool(tables["S8"]["suppressed"])].copy()
    pairs["observed_expected_ratio"] = pd.to_numeric(
        pairs["observed_expected_ratio"], errors="coerce"
    )
    pairs = pairs.sort_values(
        ["observed_expected_ratio", "pathogen_a", "pathogen_b"],
        ascending=[False, True, True],
        na_position="last",
    )
    environment = (
        tables["S10"]
        .groupby("period", observed=True)
        .agg(
            cells=("q_period", "size"),
            estimable=("q_period", "count"),
            fdr05=("q_period", lambda values: int(pd.to_numeric(values, errors="coerce").le(0.05).sum())),
            low_n=("low_n", lambda values: int(as_bool(values).sum())),
        )
        .reset_index()
    )
    selected_rf = tables["S11"].loc[
        as_bool(tables["S11"]["selected_by_training_cv"])
    ].copy()
    selected_rf = selected_rf.sort_values("pathogen")

    outputs = {
        "table_s1_cohort_pathogen_reproduced.csv": cohort,
        "table_s2_demographics_reproduced.csv": demographics,
        "age_peak_pathogen_rates_reproduced.csv": age_peaks,
        "codetection_by_age_reproduced.csv": codetection_age,
        "joint_tested_codetection_ranked_reproduced.csv": pairs,
        "environment_period_summary_reproduced.csv": environment,
        "rf_training_cv_selected_rows_reproduced.csv": selected_rf,
    }
    for filename, frame in outputs.items():
        frame.to_csv(output_dir / filename, index=False, encoding="utf-8-sig", lineterminator="\n")

    summary = {
        "status": "PASS",
        "mode": "privacy_safe_aggregate_tables",
        "python": platform.python_version(),
        "packages": {name: version(name) for name in ("pandas", "pyreadr")},
        "verified_tables": len(check_frame),
        "derived_summaries": list(outputs),
    }
    (output_dir / "public_table_reproduction_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
