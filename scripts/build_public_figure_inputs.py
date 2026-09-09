"""Create the minimum aggregate input set needed to redraw the nine figures."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import pandas as pd


BRANCH = "minimal_corrected"

# Only columns consumed by input validation or drawing are retained.  In
# particular, numerator/denominator columns that are not drawn are omitted.
PUBLIC_COLUMNS: dict[str, tuple[str, ...]] = {
    "Fig1_monthly_phase2.csv": (
        "Year", "Month", "rate_display_pct", "virus", "branch", "phase"
    ),
    "Fig2_age_rates.csv": (
        "Age_Group", "rate_display_pct", "virus", "branch", "phase"
    ),
    "Fig2_monthly_age_phase2.csv": (
        "Year", "Month", "Age_Group", "rate_display_pct", "virus", "branch", "phase"
    ),
    "Fig3_monthly_age_codetection.csv": (
        "Year", "Month", "Age_Group", "rate_display_pct", "branch"
    ),
    "Fig3_age_codetection.csv": (
        "Age_Group", "codetection_rate_pct", "branch"
    ),
    # Exact pair counts are the values displayed by the author-designated Fig. 3.
    # They are aggregate, but cells below 20 require explicit author/institutional
    # approval before public upload.
    "Fig3_codetection_matrix.csv": (
        "virus_a", "virus_b", "joint_tested", "both_positive", "expected_count",
        "observed_expected_ratio", "branch"
    ),
    "Fig4_predictions.csv": (
        "branch", "virus", "lag", "date", "actual", "predicted", "persistence",
        "seasonal_naive", "selected_for_plot"
    ),
    "Fig4_selected_feature_importance.csv": (
        "branch", "virus", "lag", "feature", "importance", "rank"
    ),
    "Fig4_selected_performance.csv": (
        "branch", "virus", "selected_lag", "best_baseline",
        "delta_rmse_rf_minus_baseline", "rf_rmse_common", "rf_r2_common",
        "baseline_rmse_common", "baseline_r2_common", "ci95_low", "ci95_high",
        "rf_win_ci95"
    ),
    "FigS1_monthly_rates.csv": (
        "Year", "Month", "rate_display_pct", "virus", "branch", "phase"
    ),
    "FigS2_environment_cells.csv": (
        "branch", "period", "age_group", "virus", "environment", "n_months",
        "rho_display", "q_period", "fdr05_period", "status"
    ),
    "FigS3_environment_cells.csv": (
        "branch", "period", "age_group", "virus", "environment", "n_months",
        "rho_display", "q_period", "fdr05_period", "status"
    ),
    "FigS4_environment_cells.csv": (
        "branch", "period", "age_group", "virus", "environment", "n_months",
        "rho_display", "q_period", "fdr05_period", "status"
    ),
    "FigS5_rf_lag_metrics.csv": (
        "branch", "virus", "lag", "rf_rmse_common", "rf_r2_common",
        "best_baseline", "best_baseline_rmse_common", "best_baseline_r2_common",
        "delta_r2_rf_minus_best", "selected_by_training_cv", "legacy_test_selected"
    ),
}


FORBIDDEN_COLUMN = re.compile(
    r"(^|_)(?:id|patient|sample|record|mrn|name|姓名|病案|样本)(?:_|$)", re.IGNORECASE
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    source_dir = args.source_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for filename, columns in PUBLIC_COLUMNS.items():
        source = source_dir / filename
        frame = pd.read_csv(source, encoding="utf-8-sig")
        if "branch" not in frame:
            raise ValueError(f"Missing branch column: {filename}")
        frame = frame.loc[frame["branch"].astype(str).eq(BRANCH)].copy()
        missing = [column for column in columns if column not in frame]
        if missing:
            raise ValueError(f"{filename}: missing columns {missing}")
        frame = frame.loc[:, list(columns)]
        forbidden = [column for column in frame if FORBIDDEN_COLUMN.search(str(column))]
        if forbidden:
            raise ValueError(f"{filename}: identifier-like columns {forbidden}")
        if frame.empty:
            raise ValueError(f"{filename}: no {BRANCH} rows")

        temporal_resolution = "none"
        if {"Year", "Month"}.issubset(frame.columns):
            temporal_resolution = "calendar_month"
            if "Day" in frame.columns:
                raise ValueError(f"{filename}: unexpected day-level column")
        if "date" in frame:
            dates = pd.to_datetime(frame["date"], errors="raise")
            if not dates.dt.dayofweek.eq(6).all():
                raise ValueError(f"{filename}: dates are not Sunday week anchors")
            temporal_resolution = "weekly_aggregate_sunday_anchor"

        small_cells = 0
        if "both_positive" in frame:
            counts = pd.to_numeric(frame["both_positive"], errors="coerce")
            small_cells = int(counts.between(1, 19, inclusive="both").sum())

        destination = output_dir / filename
        frame.to_csv(destination, index=False, encoding="utf-8-sig", lineterminator="\n")
        rows.append(
            {
                "file": filename,
                "rows": len(frame),
                "columns": len(frame.columns),
                "sha256": sha256(destination),
                "branch": BRANCH,
                "temporal_resolution": temporal_resolution,
                "direct_identifier_columns": "",
                "individual_exact_dates": False,
                "displayed_small_pair_count_cells": small_cells,
                "privacy_status": (
                    "AUTHOR_INSTITUTION_REVIEW_REQUIRED_FIG3"
                    if small_cells or filename == "Fig3_monthly_age_codetection.csv"
                    else "PASS_AGGREGATE_MINIMIZED"
                ),
            }
        )

    manifest = pd.DataFrame(rows)
    manifest_path = args.manifest.resolve()
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_path, index=False, encoding="utf-8-sig", lineterminator="\n")
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
