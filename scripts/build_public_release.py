"""Build the public, aggregate-only supplementary data release.

The script is deliberately portable: all paths are supplied on the command
line, all RDA files are written by Python/pyreadr, and source files are never
modified.  Each logical supplementary table is written as one CSV and one
single-data-frame RDA file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import pandas as pd
import pyreadr


@dataclass(frozen=True)
class TableSpec:
    number: str
    object_name: str
    title: str
    source_file: str


TABLES = [
    TableSpec("S1", "cohort_pathogen_summary", "Cohort and pathogen-specific summary by diagnostic phase", "Supplementary_Table_S1_cohort_pathogen_summary.csv"),
    TableSpec("S2", "demographics_summary", "Aggregate demographic characteristics by diagnostic phase", "Supplementary_Table_S1_demographics.csv"),
    TableSpec("S3", "panel_history", "Observed diagnostic-panel signatures and target-specific denominators", "Supplementary_Table_S1_panel_history.csv"),
    TableSpec("S4", "complete_year_sensitivity", "Complete-calendar-year and boundary sensitivity analysis", "Supplementary_Table_S1_complete_year_sensitivity.csv"),
    TableSpec("S5", "age_standardized_sensitivity", "Age-standardized cross-phase sensitivity analysis", "Supplementary_Table_S1_age_standardized_sensitivity.csv"),
    TableSpec("S6", "age_pathogen_rates", "Age-stratified pathogen positivity", "Supplementary_Table_S2_age_pathogen_rates.csv"),
    TableSpec("S7", "codetection_by_age", "Co-detection by age group", "Supplementary_Table_S2_codetection_by_age.csv"),
    TableSpec("S8", "joint_tested_codetection", "Pairwise co-detection among jointly tested events", "Supplementary_Table_S2_joint_tested_codetection.csv"),
    TableSpec("S9", "panel_size_denominators", "Co-detection by number of targets tested", "Supplementary_Table_S2_panel_denominators.csv"),
    TableSpec("S10", "environment_spearman", "Exploratory monthly ecological Spearman correlations", "Supplementary_Table_S3_environment_spearman.csv"),
    TableSpec("S11", "rf_lag_metrics", "Random-forest performance versus prespecified baselines", "Supplementary_Table_S3_rf_lag_metrics.csv"),
]


EXPECTED_SCHEMAS = {
    "cohort_pathogen_summary": ["phase", "pathogen", "n_events", "n_test", "n_pos", "rate_pct", "ci95_low_pct", "ci95_high_pct", "suppressed", "note"],
    "demographics_summary": ["phase", "n_events", "age_median_years", "age_q1_years", "age_q3_years", "male_n", "male_pct", "female_n", "female_pct", "unknown_sex_n"],
    "panel_history": ["phase", "panel_signature", "panel_size", "targets", "start_month", "end_month", "n_events", "pathogen", "n_test", "rare_signature_merged", "suppressed", "note"],
    "complete_year_sensitivity": ["phase", "window", "pathogen", "n_test", "n_pos", "rate_pct", "ci95_low_pct", "ci95_high_pct"],
    "age_standardized_sensitivity": ["standard_population", "pathogen", "phase", "standardized_rate_pct", "ci95_low_pct", "ci95_high_pct", "method", "n_bootstrap", "seed", "note"],
    "age_pathogen_rates": ["phase", "age_group", "pathogen", "n_test", "n_pos", "rate_pct", "ci95_low_pct", "ci95_high_pct", "low_n_lt20", "suppressed"],
    "codetection_by_age": ["age_group", "valid_events", "codetected", "codetection_rate_pct", "ci95_low_pct", "ci95_high_pct", "suppressed"],
    "joint_tested_codetection": ["pathogen_a", "pathogen_b", "joint_tested", "both_positive", "rate_pct", "expected_count", "observed_expected_ratio", "suppressed"],
    "panel_size_denominators": ["n_targets_tested", "events", "codetected", "codetection_rate_pct", "suppressed"],
    "environment_spearman": ["period", "age_group", "pathogen", "environment", "n_months", "rho", "p_value", "q_period", "q_age_panel", "low_n", "status"],
    "rf_lag_metrics": ["pathogen", "lag", "test_n", "rf_rmse", "rf_r2", "persistence_rmse", "persistence_r2", "seasonal_naive_rmse", "seasonal_naive_r2", "cv_mean_rmse", "cv_mean_r2", "selected_by_training_cv"],
}


COUNT_COLUMNS = {
    "cohort_pathogen_summary": ["n_events", "n_test", "n_pos"],
    "demographics_summary": ["n_events", "male_n", "female_n", "unknown_sex_n"],
    "panel_history": ["n_events", "n_test"],
    "complete_year_sensitivity": ["n_test", "n_pos"],
    "age_pathogen_rates": ["n_test", "n_pos"],
    "codetection_by_age": ["valid_events", "codetected"],
    "joint_tested_codetection": ["joint_tested", "both_positive"],
    "panel_size_denominators": ["events", "codetected"],
}


VARIABLE_INFO = {
    "phase": ("Diagnostic phase", "category"),
    "pathogen": ("Pathogen or pooled pathogen target", "category"),
    "pathogen_a": ("First pathogen in an unordered pair", "category"),
    "pathogen_b": ("Second pathogen in an unordered pair", "category"),
    "n_events": ("Number of analytic same-day laboratory testing events", "events"),
    "n_test": ("Number of events with a reported result for the target", "events"),
    "n_pos": ("Number of target-positive events", "events"),
    "rate_pct": ("Positivity or pairwise co-detection rate", "percentage points"),
    "ci95_low_pct": ("Lower bound of the Wilson 95% confidence interval", "percentage points"),
    "ci95_high_pct": ("Upper bound of the Wilson 95% confidence interval", "percentage points"),
    "suppressed": ("Primary or complementary disclosure-control flag", "boolean"),
    "note": ("Methodological note", "text"),
    "age_median_years": ("Median age", "years"),
    "age_q1_years": ("First quartile of age", "years"),
    "age_q3_years": ("Third quartile of age", "years"),
    "male_n": ("Events recorded as male", "events"),
    "male_pct": ("Events recorded as male", "percentage points"),
    "female_n": ("Events recorded as female", "events"),
    "female_pct": ("Events recorded as female", "percentage points"),
    "unknown_sex_n": ("Events with unknown or missing sex", "events"),
    "panel_signature": ("Observed tested-target panel signature; rare signatures are pooled", "category"),
    "panel_size": ("Number of targets in the observed panel signature", "targets"),
    "targets": ("Targets represented by the panel signature", "text"),
    "start_month": ("First aggregate calendar month in which the signature was observed", "YYYY-MM"),
    "end_month": ("Last aggregate calendar month in which the signature was observed", "YYYY-MM"),
    "rare_signature_merged": ("Whether multiple signatures with fewer than 20 events were pooled", "boolean"),
    "window": ("Calendar-year sensitivity-analysis window", "text"),
    "standard_population": ("Reference age distribution for direct standardization", "text"),
    "standardized_rate_pct": ("Directly age-standardized positivity", "percentage points"),
    "method": ("Estimation method", "text"),
    "n_bootstrap": ("Number of bootstrap replicates", "replicates"),
    "seed": ("Random seed", "integer"),
    "age_group": ("Prespecified age category", "category"),
    "low_n_lt20": ("Target-specific tested denominator below 20", "boolean"),
    "valid_events": ("Events with at least one reported target", "events"),
    "codetected": ("Events with at least two positive targets", "events"),
    "codetection_rate_pct": ("Co-detection rate", "percentage points"),
    "joint_tested": ("Events jointly tested for both pathogens", "events"),
    "both_positive": ("Jointly tested events positive for both pathogens", "events"),
    "expected_count": ("Expected both-positive count under marginal independence", "events"),
    "observed_expected_ratio": ("Observed divided by expected both-positive count", "ratio"),
    "n_targets_tested": ("Number of reported pathogen targets in an event", "targets"),
    "events": ("Eligible events", "events"),
    "period": ("Prespecified aggregate calendar period", "category"),
    "environment": ("Environmental variable", "category"),
    "n_months": ("Pairwise-complete calendar months", "months"),
    "rho": ("Spearman rank correlation coefficient", "dimensionless"),
    "p_value": ("Two-sided nominal Spearman P value", "probability"),
    "q_period": ("Benjamini-Hochberg q value within period", "probability"),
    "q_age_panel": ("Benjamini-Hochberg q value within age panel", "probability"),
    "low_n": ("Correlation did not meet the minimum paired-month rule", "boolean"),
    "status": ("Computation status", "category"),
    "lag": ("Model lag", "weeks"),
    "test_n": ("Held-out weeks available for the comparison", "weeks"),
    "rf_rmse": ("Random-forest root mean squared error", "percentage points"),
    "rf_r2": ("Random-forest coefficient of determination", "dimensionless"),
    "persistence_rmse": ("Persistence-baseline root mean squared error", "percentage points"),
    "persistence_r2": ("Persistence-baseline coefficient of determination", "dimensionless"),
    "seasonal_naive_rmse": ("Seasonal-naive baseline root mean squared error", "percentage points"),
    "seasonal_naive_r2": ("Seasonal-naive baseline coefficient of determination", "dimensionless"),
    "cv_mean_rmse": ("Training-only cross-validation mean RMSE", "percentage points"),
    "cv_mean_r2": ("Training-only cross-validation mean R-squared", "dimensionless"),
    "selected_by_training_cv": ("Lag selected using training-only cross-validation", "boolean"),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().eq("true")


def suppress_values(frame: pd.DataFrame, index: int, columns: list[str]) -> None:
    for column in columns:
        frame.at[index, column] = pd.NA
    frame.at[index, "suppressed"] = True


def apply_complementary_suppression(spec: TableSpec, frame: pd.DataFrame) -> list[dict[str, object]]:
    """Add a second hidden cell where a published subtotal isolates one hidden cell."""

    log: list[dict[str, object]] = []
    if spec.object_name == "age_pathogen_rates":
        for (phase, pathogen), group in frame.groupby(["phase", "pathogen"], sort=False, dropna=False):
            suppressed = bool_series(group["suppressed"])
            if int(suppressed.sum()) != 1:
                continue
            candidates = group.loc[~suppressed].copy()
            candidates["_n_pos"] = pd.to_numeric(candidates["n_pos"], errors="coerce")
            candidates["_n_test"] = pd.to_numeric(candidates["n_test"], errors="coerce")
            candidates = candidates.sort_values(["_n_pos", "_n_test", "age_group"], na_position="last")
            chosen = int(candidates.index[0])
            age_group = str(frame.at[chosen, "age_group"])
            suppress_values(frame, chosen, ["n_test", "n_pos", "rate_pct", "ci95_low_pct", "ci95_high_pct"])
            log.append({"table": spec.number, "group": f"{phase}|{pathogen}", "secondary_cell": age_group, "reason": "prevent subtotal subtraction from recovering the sole primary-suppressed cell"})
    elif spec.object_name == "codetection_by_age":
        suppressed = bool_series(frame["suppressed"])
        if int(suppressed.sum()) == 1:
            candidates = frame.loc[~suppressed].copy()
            candidates["_count"] = pd.to_numeric(candidates["codetected"], errors="coerce")
            candidates = candidates.sort_values(["_count", "valid_events", "age_group"], na_position="last")
            chosen = int(candidates.index[0])
            age_group = str(frame.at[chosen, "age_group"])
            suppress_values(frame, chosen, ["valid_events", "codetected", "codetection_rate_pct", "ci95_low_pct", "ci95_high_pct"])
            log.append({"table": spec.number, "group": "all age groups", "secondary_cell": age_group, "reason": "prevent overall-total subtraction from recovering the sole primary-suppressed cell"})
    return log


def write_and_roundtrip_rda(frame: pd.DataFrame, path: Path, object_name: str) -> str:
    pyreadr.write_rdata(str(path), frame, df_name=object_name, compress="gzip")
    loaded = pyreadr.read_r(str(path))
    if list(loaded) != [object_name]:
        raise AssertionError(f"{path.name}: unexpected R objects {list(loaded)}")
    actual = loaded[object_name].reset_index(drop=True)
    expected = frame.reset_index(drop=True)
    pd.testing.assert_frame_equal(actual, expected, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    return "PASS"


def privacy_findings(spec: TableSpec, frame: pd.DataFrame) -> list[str]:
    findings: list[str] = []
    forbidden_header = re.compile(r"(?:^id$|_id$|sample|patient|subject|token|mrn|medical_record|病案|age_days|date_of|event_date|birth)", re.I)
    for column in frame.columns:
        if forbidden_header.search(column):
            findings.append(f"forbidden_header:{column}")

    patterns = {
        "exact_date": re.compile(r"\b(?:19|20)\d{2}[-/]\d{1,2}[-/]\d{1,2}\b"),
        "mrn_date_composite": re.compile(r"\b[A-Za-z0-9.]+_(?:19|20)\d{2}_\d{1,2}_\d{1,2}\b"),
        "absolute_path": re.compile(r"(?:[A-Za-z]:[\\/]|\\\\[A-Za-z0-9_.-]+[\\/])"),
        "secret": re.compile(r"(?:sk-[A-Za-z0-9_-]{16,}|api[_ -]?key|bearer\s+[A-Za-z0-9._-]+)", re.I),
    }
    text = "\n".join(frame.select_dtypes(include=["object", "string"]).fillna("").astype(str).stack().tolist())
    for name, pattern in patterns.items():
        if pattern.search(text):
            findings.append(name)

    for column in COUNT_COLUMNS.get(spec.object_name, []):
        values = pd.to_numeric(frame[column], errors="coerce")
        if bool(((values > 0) & (values < 20)).any()):
            findings.append(f"visible_count_below_20:{column}")

    if spec.object_name == "age_pathogen_rates":
        for key, group in frame.groupby(["phase", "pathogen"], dropna=False):
            n_suppressed = int(bool_series(group["suppressed"]).sum())
            if n_suppressed == 1:
                findings.append("complementary_disclosure:" + "|".join(map(str, key)))
    if spec.object_name == "codetection_by_age" and int(bool_series(frame["suppressed"]).sum()) == 1:
        findings.append("complementary_disclosure:overall_age_total")
    return findings


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    source_dir = args.source_dir.resolve()
    output_dir = args.output_dir.resolve()
    data_dir = output_dir / "data"
    docs_dir = output_dir / "docs"
    qa_dir = output_dir / "qa"
    for directory in (data_dir, docs_dir, qa_dir):
        directory.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, object]] = []
    map_rows: list[dict[str, object]] = []
    dictionary_rows: list[dict[str, object]] = []
    privacy_rows: list[dict[str, object]] = []
    suppression_log: list[dict[str, object]] = []

    for spec in TABLES:
        source = source_dir / spec.source_file
        frame = pd.read_csv(source, encoding="utf-8-sig")
        expected_schema = EXPECTED_SCHEMAS[spec.object_name]
        if frame.columns.tolist() != expected_schema:
            raise ValueError(f"{spec.source_file}: schema mismatch")
        if frame.empty:
            raise ValueError(f"{spec.source_file}: empty table")
        suppression_log.extend(apply_complementary_suppression(spec, frame))

        stem = f"Table_{spec.number}_{spec.object_name}"
        csv_path = data_dir / f"{stem}.csv"
        rda_path = data_dir / f"{stem}.rda"
        frame.to_csv(csv_path, index=False, encoding="utf-8-sig", lineterminator="\n")
        roundtrip = write_and_roundtrip_rda(frame, rda_path, spec.object_name)
        findings = privacy_findings(spec, frame)
        status = "PASS" if not findings else "FAIL"
        privacy_rows.append({"table": spec.number, "object_name": spec.object_name, "rows": len(frame), "columns": len(frame.columns), "status": status, "findings": "|".join(findings), "csv_sha256": sha256(csv_path), "rda_sha256": sha256(rda_path)})
        manifest_rows.append({
            "table": spec.number,
            "object_name": spec.object_name,
            "title": spec.title,
            "rows": len(frame),
            "columns": len(frame.columns),
            "source_file": spec.source_file,
            "source_sha256": sha256(source),
            "csv_file": csv_path.relative_to(output_dir).as_posix(),
            "csv_sha256": sha256(csv_path),
            "rda_file": rda_path.relative_to(output_dir).as_posix(),
            "rda_sha256": sha256(rda_path),
            "python_rda_roundtrip": roundtrip,
            "release_scope": "public aggregate only",
        })
        map_rows.append({"table": spec.number, "title": spec.title, "object_name": spec.object_name, "source_file": spec.source_file, "rows": len(frame), "columns": len(frame.columns)})
        for column in frame.columns:
            definition, unit = VARIABLE_INFO.get(column, (column.replace("_", " "), "see method"))
            dictionary_rows.append({
                "table": spec.number,
                "object_name": spec.object_name,
                "unit_of_observation": "aggregate table row",
                "variable": column,
                "pandas_dtype": str(frame[column].dtype),
                "definition": definition,
                "unit_or_levels": unit,
                "missing_rule": "NA denotes unavailable, non-estimable, or disclosure-controlled; see suppressed/status fields",
                "privacy_class": "public aggregate",
                "source_file": spec.source_file,
            })

    manifest = pd.DataFrame(manifest_rows)
    table_map = pd.DataFrame(map_rows)
    dictionary = pd.DataFrame(dictionary_rows)
    privacy = pd.DataFrame(privacy_rows)
    suppressions = pd.DataFrame(suppression_log, columns=["table", "group", "secondary_cell", "reason"])
    manifest.to_csv(output_dir / "release_manifest.csv", index=False, encoding="utf-8-sig", lineterminator="\n")
    table_map.to_csv(docs_dir / "supplementary_table_map.csv", index=False, encoding="utf-8-sig", lineterminator="\n")
    dictionary.to_csv(docs_dir / "data_dictionary.csv", index=False, encoding="utf-8-sig", lineterminator="\n")
    privacy.to_csv(qa_dir / "privacy_scan.csv", index=False, encoding="utf-8-sig", lineterminator="\n")
    suppressions.to_csv(qa_dir / "complementary_suppression_log.csv", index=False, encoding="utf-8-sig", lineterminator="\n")

    release_info = {
        "status": "PASS" if privacy["status"].eq("PASS").all() else "FAIL",
        "scope": "11 privacy-controlled aggregate supplementary tables; no row-level records",
        "analysis_branch": "minimal_corrected",
        "repository_target": "https://github.com/lianjiabian/Respiratory_Pathogen_Dynamics",
        "python": platform.python_version(),
        "pandas": version("pandas"),
        "pyreadr": version("pyreadr"),
        "tables": len(TABLES),
        "all_python_roundtrips_pass": bool(manifest["python_rda_roundtrip"].eq("PASS").all()),
        "secondary_suppression_cells": int(len(suppressions)),
    }
    (qa_dir / "release_summary.json").write_text(json.dumps(release_info, ensure_ascii=False, indent=2), encoding="utf-8")
    if release_info["status"] != "PASS":
        raise ValueError("Public privacy scan failed; see qa/privacy_scan.csv")
    print(json.dumps(release_info, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
