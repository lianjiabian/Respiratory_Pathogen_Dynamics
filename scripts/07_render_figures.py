"""Conservative re-rendering of the submitted Figure 1-4 and Figure S1-S5.

The module has a deliberate two-step safety gate.  Its default action is static
input validation only.  Figure files are created only when BOTH ``--render``
and ``--pdf-marker-confirmed`` are supplied after the parent workflow has run
the required PDF marker.

Inputs are the frozen CSV files produced by ``04_prepare_figure_inputs.py``.
The two accepted analysis branches are:

* ``legacy_exact``: machine key for the internal legacy reference; it is a
  current-environment rerun and is not claimed to be numerically exact.
* ``minimal_corrected``: conservative submission branch.

For every figure and branch, rendering exports editable SVG/EPS/PDF plus a
600-dpi RGB LZW TIFF.  SVG text remains text (``svg.fonttype = none``), while
PDF/EPS use TrueType fonts (font type 42).  No raster layer is placed in a
vector file.
"""

from __future__ import annotations

import argparse
import io
import json
import math
import re
import string
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


BRANCHES = ("legacy_exact", "minimal_corrected")
OUTPUT_SET = {
    "legacy_exact": "legacy_reference_internal",
    "minimal_corrected": "minimal_corrected_submission",
}
FIGURE_ORDER = (
    "Figure_1",
    "Figure_2",
    "Figure_3",
    "Figure_4",
    "Figure_S1",
    "Figure_S2",
    "Figure_S3",
    "Figure_S4",
    "Figure_S5",
)
FIGURE_SIZES_INCH = {
    "Figure_1": (6.8, 5.2),
    "Figure_2": (6.8, 9.0),
    # Author-approved continuous monthly curves with a compact 3-by-3 grid.
    "Figure_3": (7.0, 7.5),
    "Figure_4": (6.65, 8.8),
    "Figure_S1": (6.8, 8.5),
    "Figure_S2": (6.8, 7.4),
    "Figure_S3": (6.8, 7.4),
    "Figure_S4": (6.8, 7.4),
    "Figure_S5": (6.8, 4.0),
}

VIRUSES_PHASE2 = ("FluA", "FluB", "RSV", "ADV", "PIV", "HRV", "MP", "HMPV", "HBOV")
VIRUSES_PHASE1 = ("FluA", "FluB", "RSV", "ADV", "PIV", "PIV1", "PIV2", "PIV3")
S5_VIRUS_ORDER = VIRUSES_PHASE2
AGE_ORDER = ("<1", "1-3", "4-12", "13-17", "18-44", "45-59", "60-74", "75-89", ">=90")
ENVIRONMENT_ORDER = ("windlevel", "temp", "humi", "PM2.5", "PM10", "NO2", "SO2", "O3", "CO", "AQI")
ENVIRONMENT_LABELS = {
    "windlevel": "Wind",
    "temp": "Temp",
    "humi": "Humidity",
    "PM2.5": "PM2.5",
    "PM10": "PM10",
    "NO2": "NO2",
    "SO2": "SO2",
    "O3": "O3",
    "CO": "CO",
    "AQI": "AQI",
}
DISPLAY_VIRUS = {
    "FluA": "FluA",
    "FluB": "FluB",
    "RSV": "RSV",
    "ADV": "ADV",
    "PIV": "HPIV",
    "PIV1": "HPIV-1",
    "PIV2": "HPIV-2",
    "PIV3": "HPIV-3",
    "HRV": "HRV",
    "MP": "MP",
    "HMPV": "HMPV",
    "HBOV": "HBoV",
}

# Exact seaborn ``husl`` colours used by the legacy monthly panels.
VIRUS_COLORS = {
    "FluA": "#f77189",
    "FluB": "#d58c32",
    "RSV": "#a4a031",
    "ADV": "#50b131",
    "PIV": "#34ae91",
    "HRV": "#37abb5",
    "MP": "#3ba3ec",
    "HMPV": "#bb83f4",
    "HBOV": "#f564d4",
    "PIV1": "#3ba3ec",
    "PIV2": "#8d7cf3",
    "PIV3": "#f564d4",
}
AGE_COLORS = {
    age: colour
    for age, colour in zip(
        AGE_ORDER,
        ("#482475", "#414487", "#355f8d", "#2a788e", "#21918c", "#22a884", "#44bf70", "#7ad151", "#bddf26"),
        strict=True,
    )
}


@dataclass(frozen=True)
class InputContract:
    filename: str
    required_columns: tuple[str, ...]
    both_branches: bool = True


INPUT_CONTRACTS = (
    InputContract(
        "Fig1_monthly_phase2.csv",
        ("Year", "Month", "tested", "positive", "rate_display_pct", "virus", "branch", "phase"),
    ),
    InputContract(
        "Fig2_age_rates.csv",
        ("Age_Group", "tested", "positive", "rate_display_pct", "virus", "branch", "phase"),
    ),
    InputContract(
        "Fig2_monthly_age_phase2.csv",
        ("Year", "Month", "Age_Group", "rate_display_pct", "virus", "branch", "phase"),
    ),
    InputContract(
        "Fig3_monthly_age_codetection.csv",
        ("Year", "Month", "Age_Group", "valid_events", "rate_display_pct", "branch"),
    ),
    InputContract(
        "Fig3_age_codetection.csv",
        (
            "Age_Group",
            "codetected",
            "valid_events",
            "codetection_rate_pct",
            "ci95_low_pct",
            "ci95_high_pct",
            "branch",
        ),
    ),
    InputContract(
        "Fig3_codetection_matrix.csv",
        (
            "virus_a",
            "virus_b",
            "joint_tested",
            "both_positive",
            "expected_count",
            "observed_expected_ratio",
            "branch",
        ),
    ),
    InputContract(
        "Fig3_panel_denominators.csv",
        ("n_targets_tested", "events", "codetected", "codetection_rate_pct", "branch"),
    ),
    InputContract(
        "Fig4_predictions.csv",
        (
            "branch",
            "virus",
            "lag",
            "date",
            "actual",
            "predicted",
            "persistence",
            "seasonal_naive",
            "selected_for_plot",
        ),
    ),
    InputContract(
        "Fig4_selected_feature_importance.csv",
        ("branch", "virus", "lag", "feature", "importance", "rank"),
    ),
    InputContract(
        "Fig4_selected_performance.csv",
        (
            "branch",
            "virus",
            "selected_lag",
            "best_baseline",
            "delta_rmse_rf_minus_baseline",
            "rf_rmse_common",
            "rf_r2_common",
            "baseline_rmse_common",
            "baseline_r2_common",
            "ci95_low",
            "ci95_high",
            "rf_win_ci95",
        ),
        both_branches=False,
    ),
    InputContract(
        "FigS1_monthly_rates.csv",
        ("Year", "Month", "rate_display_pct", "virus", "branch", "phase"),
    ),
    InputContract(
        "FigS2_environment_cells.csv",
        (
            "branch",
            "period",
            "age_group",
            "virus",
            "environment",
            "n_months",
            "rho_display",
            "q_period",
            "fdr05_period",
            "status",
        ),
    ),
    InputContract(
        "FigS3_environment_cells.csv",
        (
            "branch",
            "period",
            "age_group",
            "virus",
            "environment",
            "n_months",
            "rho_display",
            "q_period",
            "fdr05_period",
            "status",
        ),
    ),
    InputContract(
        "FigS4_environment_cells.csv",
        (
            "branch",
            "period",
            "age_group",
            "virus",
            "environment",
            "n_months",
            "rho_display",
            "q_period",
            "fdr05_period",
            "status",
        ),
    ),
    InputContract(
        "FigS5_rf_lag_metrics.csv",
        (
            "branch",
            "virus",
            "lag",
            "rf_r2_common",
            "rf_rmse_common",
            "best_baseline",
            "best_baseline_rmse_common",
            "best_baseline_r2_common",
            "delta_r2_rf_minus_best",
            "selected_by_training_cv",
            "legacy_test_selected",
        ),
    ),
)

# The public GitHub bundle deliberately keeps only columns consumed by the
# renderer.  Numerator/denominator and confidence-interval fields that are not
# drawn remain required for the internal two-branch audit, but are not needed
# to reproduce the conservative minimal_corrected figures.  The standalone
# Fig. 3 denominator table is likewise audit-only and is not read by plotting
# or validation code.
PUBLIC_MINIMAL_REQUIRED_COLUMNS: dict[str, tuple[str, ...]] = {
    "Fig1_monthly_phase2.csv":
        ("Year", "Month", "rate_display_pct", "virus", "branch", "phase"),
    "Fig2_age_rates.csv":
        ("Age_Group", "rate_display_pct", "virus", "branch", "phase"),
    "Fig3_monthly_age_codetection.csv":
        ("Year", "Month", "Age_Group", "rate_display_pct", "branch"),
    "Fig3_age_codetection.csv":
        ("Age_Group", "codetection_rate_pct", "branch"),
}
PUBLIC_MINIMAL_OPTIONAL_FILES = frozenset({"Fig3_panel_denominators.csv"})


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    mapping = {"true": True, "false": False, "1": True, "0": False, "yes": True, "no": False}
    normalised = series.astype("string").str.strip().str.lower()
    converted = normalised.map(mapping)
    if converted.isna().any():
        bad = sorted(normalised.loc[converted.isna()].dropna().unique().tolist())
        raise ValueError(f"Unrecognised Boolean values: {bad}")
    return converted.astype(bool)


def _display_age(value: str) -> str:
    return "≥90" if value == ">=90" else value.replace("-", "–")


def _require_set(actual: Iterable[str], expected: Iterable[str], label: str) -> None:
    actual_set = set(actual)
    expected_set = set(expected)
    if actual_set != expected_set:
        raise ValueError(f"{label}: expected {sorted(expected_set)}, got {sorted(actual_set)}")


def _assert_unique(frame: pd.DataFrame, keys: list[str], label: str) -> None:
    duplicate = frame.duplicated(keys, keep=False)
    if duplicate.any():
        example = frame.loc[duplicate, keys].head(5).to_dict("records")
        raise ValueError(f"{label}: duplicate keys {keys}; examples={example}")


def load_inputs(
    input_dir: Path,
    requested_branches: Iterable[str] = BRANCHES,
) -> dict[str, pd.DataFrame]:
    requested = tuple(dict.fromkeys(requested_branches))
    unknown = sorted(set(requested) - set(BRANCHES))
    if not requested or unknown:
        raise ValueError(f"Invalid requested branches: {requested}; unknown={unknown}")
    public_minimal_mode = set(requested) == {"minimal_corrected"}
    tables: dict[str, pd.DataFrame] = {}
    for contract in INPUT_CONTRACTS:
        path = input_dir / contract.filename
        if not path.is_file():
            if public_minimal_mode and contract.filename in PUBLIC_MINIMAL_OPTIONAL_FILES:
                continue
            raise FileNotFoundError(f"Missing input CSV: {path}")
        frame = pd.read_csv(path)
        required_columns = (
            PUBLIC_MINIMAL_REQUIRED_COLUMNS.get(contract.filename, contract.required_columns)
            if public_minimal_mode
            else contract.required_columns
        )
        missing = sorted(set(required_columns) - set(frame.columns))
        if missing:
            raise ValueError(f"{contract.filename}: missing columns {missing}")
        if frame.empty:
            raise ValueError(f"{contract.filename}: empty input")
        branches = set(frame["branch"].dropna().astype(str))
        if contract.both_branches:
            if set(requested) == set(BRANCHES):
                valid = branches == set(BRANCHES)
                expected_description = sorted(BRANCHES)
            else:
                valid = set(requested).issubset(branches) and branches.issubset(set(BRANCHES))
                expected_description = [
                    f"at least {sorted(requested)}",
                    f"and no branches outside {sorted(BRANCHES)}",
                ]
        else:
            valid = branches == {"minimal_corrected"}
            expected_description = ["minimal_corrected"]
        if not valid:
            raise ValueError(
                f"{contract.filename}: expected branches {expected_description}, "
                f"got {sorted(branches)}"
            )
        tables[contract.filename] = frame
    return tables


def validate_inputs(
    tables: dict[str, pd.DataFrame],
    requested_branches: Iterable[str] = BRANCHES,
) -> dict[str, object]:
    """Validate frozen inputs without creating a Matplotlib Figure."""

    active_branches = tuple(dict.fromkeys(requested_branches))
    unknown = sorted(set(active_branches) - set(BRANCHES))
    if not active_branches or unknown:
        raise ValueError(f"Invalid requested branches: {active_branches}; unknown={unknown}")

    fig1 = tables["Fig1_monthly_phase2.csv"]
    _require_set(fig1["virus"].unique(), VIRUSES_PHASE2, "Figure 1 pathogens")
    _assert_unique(fig1, ["branch", "virus", "Year", "Month"], "Figure 1")
    if set(fig1["phase"].astype(str)) != {"phase2"}:
        raise ValueError("Figure 1 must contain phase2 only")

    age_rates = tables["Fig2_age_rates.csv"]
    _assert_unique(age_rates, ["branch", "phase", "virus", "Age_Group"], "Figure 2 age rates")
    for branch in active_branches:
        for phase, pathogens in (("phase1", VIRUSES_PHASE1), ("phase2", VIRUSES_PHASE2)):
            subset = age_rates.loc[(age_rates["branch"] == branch) & (age_rates["phase"] == phase)]
            _require_set(subset["virus"].unique(), pathogens, f"Figure 2 {branch} {phase} pathogens")
            _require_set(subset["Age_Group"].unique(), AGE_ORDER, f"Figure 2 {branch} {phase} ages")

    fig2_monthly = tables["Fig2_monthly_age_phase2.csv"]
    _assert_unique(
        fig2_monthly,
        ["branch", "virus", "Age_Group", "Year", "Month"],
        "Figure 2 monthly age trends",
    )
    _require_set(fig2_monthly["virus"].unique(), VIRUSES_PHASE2, "Figure 2 monthly pathogens")
    _require_set(fig2_monthly["Age_Group"].unique(), AGE_ORDER, "Figure 2 monthly ages")

    fig3_monthly = tables["Fig3_monthly_age_codetection.csv"]
    _assert_unique(fig3_monthly, ["branch", "Age_Group", "Year", "Month"], "Figure 3 monthly co-detection")
    _require_set(fig3_monthly["Age_Group"].unique(), AGE_ORDER, "Figure 3 monthly ages")
    fig3_age = tables["Fig3_age_codetection.csv"]
    _assert_unique(fig3_age, ["branch", "Age_Group"], "Figure 3 age co-detection")
    _require_set(fig3_age["Age_Group"].unique(), AGE_ORDER, "Figure 3 age bars")

    pair_table = tables["Fig3_codetection_matrix.csv"]
    _assert_unique(pair_table, ["branch", "virus_a", "virus_b"], "Figure 3 joint-tested matrix")
    for branch in active_branches:
        subset = pair_table.loc[pair_table["branch"] == branch]
        if len(subset) != len(VIRUSES_PHASE2) * (len(VIRUSES_PHASE2) - 1):
            raise ValueError(f"Figure 3 {branch}: expected 72 directed off-diagonal pathogen pairs")
        if (subset["virus_a"] == subset["virus_b"]).any():
            raise ValueError(f"Figure 3 {branch}: matrix diagonal must be absent")
        reciprocal = subset.merge(
            subset,
            left_on=["virus_a", "virus_b"],
            right_on=["virus_b", "virus_a"],
            suffixes=("", "_reciprocal"),
            validate="one_to_one",
        )
        for column in ("joint_tested", "both_positive", "expected_count", "observed_expected_ratio"):
            a = pd.to_numeric(reciprocal[column], errors="coerce").to_numpy(float)
            b = pd.to_numeric(reciprocal[f"{column}_reciprocal"], errors="coerce").to_numpy(float)
            if not np.allclose(a, b, equal_nan=True, rtol=1e-10, atol=1e-10):
                raise ValueError(f"Figure 3 {branch}: asymmetric {column}")

    predictions = tables["Fig4_predictions.csv"].copy()
    predictions["selected_for_plot"] = _as_bool(predictions["selected_for_plot"])
    predictions["date"] = pd.to_datetime(predictions["date"], errors="raise")
    _assert_unique(predictions, ["branch", "virus", "lag", "date"], "Figure 4 predictions")
    for branch in active_branches:
        selected = predictions.loc[(predictions["branch"] == branch) & predictions["selected_for_plot"]]
        _require_set(selected["virus"].unique(), VIRUSES_PHASE2, f"Figure 4 {branch} selected pathogens")
        lags_per_virus = selected.groupby("virus", observed=True)["lag"].nunique()
        if not lags_per_virus.eq(1).all():
            raise ValueError(f"Figure 4 {branch}: exactly one selected lag is required per pathogen")

    importance = tables["Fig4_selected_feature_importance.csv"]
    _assert_unique(importance, ["branch", "virus", "feature"], "Figure 4 feature importance")
    for branch in active_branches:
        subset = importance.loc[importance["branch"] == branch]
        _require_set(subset["virus"].unique(), VIRUSES_PHASE2, f"Figure 4 {branch} importance pathogens")
        if not subset.groupby("virus", observed=True)["lag"].nunique().eq(1).all():
            raise ValueError(f"Figure 4 {branch}: importance must use one fitted lag per pathogen")
        if (pd.to_numeric(subset["importance"], errors="coerce") < 0).any():
            raise ValueError(f"Figure 4 {branch}: negative feature importance")

    performance = tables["Fig4_selected_performance.csv"]
    _require_set(performance["virus"].unique(), VIRUSES_PHASE2, "Figure 4 corrected performance pathogens")
    _assert_unique(performance, ["branch", "virus"], "Figure 4 corrected performance")
    if "minimal_corrected" in active_branches:
        selected_minimal = (
            predictions.loc[
                (predictions["branch"] == "minimal_corrected")
                & predictions["selected_for_plot"]
            ]
            .groupby("virus", observed=True)["lag"]
            .first()
        )
        perf_lags = performance.set_index("virus")["selected_lag"].astype(int)
        if not selected_minimal.astype(int).sort_index().equals(perf_lags.sort_index()):
            raise ValueError("Figure 4 minimal_corrected selected lags do not match performance table")

    s1 = tables["FigS1_monthly_rates.csv"]
    _assert_unique(s1, ["branch", "phase", "virus", "Year", "Month"], "Figure S1")
    for phase, pathogens in (("phase1", VIRUSES_PHASE1), ("phase2", VIRUSES_PHASE2)):
        _require_set(s1.loc[s1["phase"] == phase, "virus"].unique(), pathogens, f"Figure S1 {phase}")

    environment_periods = {
        "FigS2_environment_cells.csv": "2019",
        "FigS3_environment_cells.csv": "2020-2022",
        "FigS4_environment_cells.csv": "2023-2025",
    }
    environment_summary: dict[str, dict[str, int]] = {}
    for filename, period in environment_periods.items():
        frame = tables[filename].copy()
        if set(frame["period"].astype(str)) != {period}:
            raise ValueError(f"{filename}: expected period {period}")
        _assert_unique(frame, ["branch", "age_group", "virus", "environment"], filename)
        _require_set(frame["age_group"].unique(), AGE_ORDER, f"{filename} ages")
        _require_set(frame["virus"].unique(), VIRUSES_PHASE2, f"{filename} pathogens")
        _require_set(frame["environment"].unique(), ENVIRONMENT_ORDER, f"{filename} exposures")
        rho = pd.to_numeric(frame["rho_display"], errors="coerce")
        if ((rho < -1) | (rho > 1)).fillna(False).any():
            raise ValueError(f"{filename}: Spearman rho outside [-1, 1]")
        q = pd.to_numeric(frame["q_period"], errors="coerce")
        if ((q < 0) | (q > 1)).fillna(False).any():
            raise ValueError(f"{filename}: BH-FDR q outside [0, 1]")
        environment_summary[filename] = {
            "rows": int(len(frame)),
            "blank_cells": int(rho.isna().sum()),
            "fdr05_cells": int(_as_bool(frame["fdr05_period"]).sum()),
        }

    s5 = tables["FigS5_rf_lag_metrics.csv"].copy()
    s5["selected_by_training_cv"] = _as_bool(s5["selected_by_training_cv"])
    s5["legacy_test_selected"] = _as_bool(s5["legacy_test_selected"])
    _assert_unique(s5, ["branch", "virus", "lag"], "Figure S5")
    _require_set(s5["virus"].unique(), VIRUSES_PHASE2, "Figure S5 pathogens")
    if set(pd.to_numeric(s5["lag"], errors="raise").astype(int)) != {0, 1, 2, 3, 4}:
        raise ValueError("Figure S5 must include forecast lags 0-4")
    for branch in active_branches:
        subset = s5.loc[s5["branch"] == branch]
        if not subset.groupby("virus", observed=True).size().eq(5).all():
            raise ValueError(f"Figure S5 {branch}: expected five lag cells per pathogen")
    if "minimal_corrected" in active_branches and not s5.loc[
        s5["branch"] == "minimal_corrected"
    ].groupby("virus")["selected_by_training_cv"].sum().eq(1).all():
        raise ValueError("Figure S5 minimal_corrected: one training-CV-selected lag required per pathogen")
    if "legacy_exact" in active_branches and not s5.loc[
        s5["branch"] == "legacy_exact"
    ].groupby("virus")["legacy_test_selected"].sum().eq(1).all():
        raise ValueError("Figure S5 legacy_exact: one legacy test-selected lag required per pathogen")

    panel_letters = {
        "Figure_1": "A-I",
        "Figure_2": "A-K",
        "Figure_3": "A-K",
        "Figure_4": "A-R",
        "Figure_S1": "A-J",
        "Figure_S2": "A-I",
        "Figure_S3": "A-I",
        "Figure_S4": "A-I",
        "Figure_S5": "single matrix",
    }
    for figure, (width, height) in FIGURE_SIZES_INCH.items():
        if width > 7 or height > 9:
            raise ValueError(f"{figure}: declared size {width}x{height} exceeds 7x9 inches")

    source = Path(__file__).read_text(encoding="utf-8")
    literal_sizes = [
        float(value)
        for value in re.findall(r"(?:fontsize|title_fontsize|labelsize)\s*=\s*([0-9]+(?:\.[0-9]+)?)", source)
    ]
    literal_sizes.extend(
        float(value)
        for value in re.findall(
            r'"(?:font\.size|axes\.titlesize|axes\.labelsize|xtick\.labelsize|ytick\.labelsize|legend\.fontsize)"\s*:\s*([0-9]+(?:\.[0-9]+)?)',
            source,
        )
    )
    source_min_font = min(literal_sizes)
    if source_min_font < 7:
        raise ValueError(f"Visible-font source gate failed: minimum literal size is {source_min_font} pt")

    return {
        "status": "PASS",
        "mode": "static_validation_only_no_figure_created",
        "branches": list(active_branches),
        "input_files": {name: int(len(frame)) for name, frame in tables.items()},
        "environment": environment_summary,
        "panel_letters": panel_letters,
        "figure_sizes_inches": {key: list(value) for key, value in FIGURE_SIZES_INCH.items()},
        "expected_pdf_count": len(active_branches) * len(FIGURE_ORDER),
        "expected_total_files": len(active_branches) * len(FIGURE_ORDER) * 4,
        "vector_text_contract": {"svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42},
        "tiff_contract": {"dpi": 600, "mode": "RGB", "compression": "tiff_lzw"},
        "font_contract": {
            "minimum_visible_font_pt": 7.0,
            "source_literal_minimum_pt": source_min_font,
            "source_static_scan": "PASS",
        },
    }


def _activate_matplotlib():
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors
    from matplotlib import patches

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.0,
            "axes.titlesize": 7.5,
            "axes.labelsize": 7.0,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 7.0,
            "axes.linewidth": 0.6,
            "lines.linewidth": 1.05,
            "grid.linewidth": 0.35,
            "grid.color": "#d9d9d9",
            "grid.linestyle": ":",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "pdf.compression": 9,
        }
    )
    return plt, mdates, mcolors, patches


def _panel_label(ax, label: str) -> None:
    ax.text(-0.10, 1.035, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.2, fontweight="bold")


def _shade_npi(ax, start="2020-01-01", end="2022-12-31") -> None:
    ax.axvspan(pd.Timestamp(start), pd.Timestamp(end), color="#eeeeee", zorder=0)


def _style_line_axis(ax, show_y: bool = True, show_x: bool = True) -> None:
    ax.grid(True, axis="both")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if not show_y:
        ax.set_ylabel("")
    if not show_x:
        ax.tick_params(axis="x", labelbottom=False)


def _dates(frame: pd.DataFrame) -> pd.Series:
    return pd.to_datetime(
        {"year": pd.to_numeric(frame["Year"], errors="raise").astype(int), "month": pd.to_numeric(frame["Month"], errors="raise").astype(int), "day": 1}
    )


def _plot_figure_1(branch: str, tables: dict[str, pd.DataFrame]):
    plt, mdates, _, _ = _activate_matplotlib()
    data = tables["Fig1_monthly_phase2.csv"].loc[lambda x: x["branch"] == branch].copy()
    data["date"] = _dates(data)
    fig, axes = plt.subplots(3, 3, figsize=FIGURE_SIZES_INCH["Figure_1"], sharex=True)
    for index, (ax, virus) in enumerate(zip(axes.flat, VIRUSES_PHASE2, strict=True)):
        subset = data.loc[data["virus"] == virus].sort_values("date")
        _shade_npi(ax)
        ax.plot(subset["date"], pd.to_numeric(subset["rate_display_pct"], errors="coerce"), color=VIRUS_COLORS[virus])
        ax.set_title(DISPLAY_VIRUS[virus], fontweight="bold", pad=2)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("Positive rate (%)" if index % 3 == 0 else "")
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.tick_params(axis="x", rotation=45)
        _style_line_axis(ax, show_y=index % 3 == 0, show_x=index >= 6)
        _panel_label(ax, string.ascii_uppercase[index])
    fig.suptitle("Monthly positive rates, 2019–2025 (PCR)", y=0.992, fontsize=9)
    fig.text(0.5, 0.008, "Grey shading: 2020–2022; months with fewer than 20 tests are not plotted.", ha="center", fontsize=7)
    fig.subplots_adjust(left=0.08, right=0.985, top=0.925, bottom=0.10, hspace=0.45, wspace=0.31)
    return fig


def _plot_age_lines(ax, subset: pd.DataFrame, viruses: Iterable[str], broken_component: str | None = None) -> list:
    handles = []
    for virus in viruses:
        virus_data = subset.loc[subset["virus"] == virus].copy()
        virus_data["Age_Group"] = pd.Categorical(virus_data["Age_Group"], AGE_ORDER, ordered=True)
        virus_data = virus_data.sort_values("Age_Group")
        line = ax.plot(
            np.arange(len(AGE_ORDER)),
            pd.to_numeric(virus_data["rate_display_pct"], errors="coerce"),
            marker="o",
            markersize=2.2,
            linewidth=0.9,
            color=VIRUS_COLORS[virus],
            label=DISPLAY_VIRUS[virus],
        )[0]
        handles.append(line)
    ax.set_xlim(-0.25, len(AGE_ORDER) - 0.75)
    ax.set_xticks(np.arange(len(AGE_ORDER)), [_display_age(value) for value in AGE_ORDER], rotation=35, ha="right")
    ax.grid(True)
    ax.spines["right"].set_visible(False)
    if broken_component != "upper":
        ax.set_xlabel("Age group")
    return handles


def _broken_limits(values: pd.Series) -> tuple[tuple[float, float], tuple[float, float]]:
    finite = np.sort(pd.to_numeric(values, errors="coerce").dropna().unique())
    if len(finite) < 3:
        return (0.0, 1.0), (1.1, 2.0)
    gaps = np.diff(finite)
    index = int(np.argmax(gaps))
    lower_value, upper_value = float(finite[index]), float(finite[index + 1])
    gap = upper_value - lower_value
    lower = (0.0, lower_value + 0.15 * gap)
    upper = (upper_value - 0.15 * gap, float(finite[-1]) * 1.08 + 0.05)
    return lower, upper


def _plot_figure_2(branch: str, tables: dict[str, pd.DataFrame]):
    plt, mdates, _, _ = _activate_matplotlib()
    age_rates = tables["Fig2_age_rates.csv"].loc[lambda x: x["branch"] == branch].copy()
    monthly = tables["Fig2_monthly_age_phase2.csv"].loc[lambda x: x["branch"] == branch].copy()
    monthly["date"] = _dates(monthly)
    fig = plt.figure(figsize=FIGURE_SIZES_INCH["Figure_2"])
    grid = fig.add_gridspec(
        6,
        6,
        height_ratios=[0.44, 1.24, 0.26, 1.0, 1.0, 1.0],
        hspace=0.62,
        wspace=0.55,
    )

    split = grid[0:2, 0:3].subgridspec(2, 1, height_ratios=[0.42, 1.0], hspace=0.045)
    ax_a_upper = fig.add_subplot(split[0, 0])
    ax_a_lower = fig.add_subplot(split[1, 0], sharex=ax_a_upper)
    phase1 = age_rates.loc[age_rates["phase"] == "phase1"]
    handles_a = _plot_age_lines(ax_a_upper, phase1, VIRUSES_PHASE1, broken_component="upper")
    _plot_age_lines(ax_a_lower, phase1, VIRUSES_PHASE1, broken_component="lower")
    lower_ylim, upper_ylim = _broken_limits(phase1["rate_display_pct"])
    ax_a_lower.set_ylim(*lower_ylim)
    ax_a_upper.set_ylim(*upper_ylim)
    ax_a_upper.spines["bottom"].set_visible(False)
    ax_a_lower.spines["top"].set_visible(False)
    ax_a_upper.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    diagonal = 0.012
    kwargs = dict(color="black", clip_on=False, linewidth=0.6)
    ax_a_upper.plot((-diagonal, +diagonal), (-diagonal, +diagonal), transform=ax_a_upper.transAxes, **kwargs)
    ax_a_upper.plot((1 - diagonal, 1 + diagonal), (-diagonal, +diagonal), transform=ax_a_upper.transAxes, **kwargs)
    ax_a_lower.plot((-diagonal, +diagonal), (1 - diagonal, 1 + diagonal), transform=ax_a_lower.transAxes, **kwargs)
    ax_a_lower.plot((1 - diagonal, 1 + diagonal), (1 - diagonal, 1 + diagonal), transform=ax_a_lower.transAxes, **kwargs)
    ax_a_lower.set_ylabel("Positive rate (%)")
    ax_a_upper.set_title("2013–2018 (IgM IFA): positive rate by age group", pad=3)
    _panel_label(ax_a_upper, "A")

    ax_b = fig.add_subplot(grid[0:2, 3:6])
    phase2 = age_rates.loc[age_rates["phase"] == "phase2"]
    handles_b = _plot_age_lines(ax_b, phase2, VIRUSES_PHASE2)
    ax_b.set_ylim(bottom=0)
    ax_b.set_ylabel("Positive rate (%)")
    ax_b.set_title("2019–2025 (PCR): positive rate by age group", pad=3)
    _panel_label(ax_b, "B")

    # Dedicated legend strips prevent the pathogen legends from covering the
    # age plotted age-rate curves in panels A and B.
    ax_legend_a = fig.add_subplot(grid[2, 0:3])
    ax_legend_a.legend(
        handles=handles_a,
        labels=[handle.get_label() for handle in handles_a],
        loc="center",
        frameon=False,
        ncol=4,
        fontsize=7,
        handlelength=1.2,
        columnspacing=0.7,
    )
    ax_legend_a.set_axis_off()
    ax_legend_b = fig.add_subplot(grid[2, 3:6])
    ax_legend_b.legend(
        handles=handles_b,
        labels=[handle.get_label() for handle in handles_b],
        loc="center",
        frameon=False,
        ncol=5,
        fontsize=7,
        handlelength=1.2,
        columnspacing=0.7,
    )
    ax_legend_b.set_axis_off()

    bottom_axes = []
    for index, virus in enumerate(VIRUSES_PHASE2):
        row = 3 + index // 3
        col = (index % 3) * 2
        ax = fig.add_subplot(grid[row, col : col + 2])
        bottom_axes.append(ax)
        _shade_npi(ax)
        virus_data = monthly.loc[monthly["virus"] == virus]
        for age in AGE_ORDER:
            subset = virus_data.loc[virus_data["Age_Group"] == age].sort_values("date")
            ax.plot(subset["date"], pd.to_numeric(subset["rate_display_pct"], errors="coerce"), color=AGE_COLORS[age], linewidth=0.75)
        ax.set_title(DISPLAY_VIRUS[virus], pad=1.5)
        ax.set_ylim(bottom=0)
        ax.set_xlim(pd.Timestamp("2019-01-01"), pd.Timestamp("2026-01-01"))
        ax.set_ylabel("Positive rate (%)" if index % 3 == 0 else "")
        ax.set_xticks(pd.date_range("2020-01-01", "2026-01-01", freq="2YS"))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.tick_params(axis="x", rotation=35)
        _style_line_axis(ax, show_y=index % 3 == 0, show_x=index >= 6)
        _panel_label(ax, string.ascii_uppercase[index + 2])
    age_handles = [bottom_axes[0].plot([], [], color=AGE_COLORS[age], label=_display_age(age))[0] for age in AGE_ORDER]
    fig.legend(age_handles, [_display_age(age) for age in AGE_ORDER], title="Age group", loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.012), fontsize=7, title_fontsize=7)
    fig.suptitle("Age-stratified pathogen positivity", y=0.995, fontsize=9)
    fig.text(0.5, 0.078, "Grey shading: 2020–2022; monthly cells with fewer than 20 tests are omitted.", ha="center", fontsize=7)
    fig.subplots_adjust(left=0.085, right=0.980, top=0.95, bottom=0.145)
    return fig


def _plot_figure_3_composition(branch: str, tables: dict[str, pd.DataFrame]):
    """Rebuild Fig. 3 in the exact visual grammar of the designated prior PDF.

    Only values come from the corrected frozen inputs.  The tall page, purple
    3 x 3 trend grid, equal-width orange lower panels, and count heatmap are
    intentionally retained.  Terminology is updated from co-infection to
    co-detection so that the artwork remains consistent with the revised
    Methods.
    """

    plt, mdates, _, _ = _activate_matplotlib()
    monthly = tables["Fig3_monthly_age_codetection.csv"].loc[lambda x: x["branch"] == branch].copy()
    monthly["date"] = _dates(monthly)
    by_age = tables["Fig3_age_codetection.csv"].loc[lambda x: x["branch"] == branch].copy()
    matrix = tables["Fig3_codetection_matrix.csv"].loc[lambda x: x["branch"] == branch].copy()

    fig = plt.figure(figsize=FIGURE_SIZES_INCH["Figure_3"])
    # The reference PDF contains a deliberate lower blank margin and a
    # 1.65:1 upper-to-lower visual height ratio.  Fixed normalized rectangles
    # reproduce that page architecture without embedding either source raster.
    upper = fig.add_gridspec(
        3,
        3,
        left=0.110,
        right=0.977,
        bottom=0.500,
        top=0.885,
        hspace=0.20,
        wspace=0.15,
    )
    for index, age in enumerate(AGE_ORDER):
        row = index // 3
        col = index % 3
        ax = fig.add_subplot(upper[row, col])
        subset = monthly.loc[monthly["Age_Group"] == age].sort_values("date")
        ax.axvspan(pd.Timestamp("2020-01-01"), pd.Timestamp("2022-12-31"), color="gray", alpha=0.1, zorder=0)
        ax.plot(
            subset["date"],
            pd.to_numeric(subset["rate_display_pct"], errors="coerce"),
            color="#800080",
            linewidth=1.15,
        )
        ax.set_ylim(0, 100)
        ax.set_xlim(pd.Timestamp("2019-01-01"), pd.Timestamp("2026-01-01"))
        ax.set_title(f"Age: {_display_age(age)}", pad=1.5, fontsize=7)
        ax.set_ylabel("Co-detection rate (%)" if col == 0 else "")
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.tick_params(axis="x", rotation=0, labelbottom=row == 2)
        ax.grid(True, linestyle="--", alpha=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if row != 2:
            ax.set_xlabel("")
        else:
            ax.set_xlabel("Date")
        ax.text(
            -0.16,
            1.03,
            string.ascii_uppercase[index],
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=12,
            fontweight="bold",
        )

    lower = fig.add_gridspec(
        1,
        2,
        left=0.06,
        right=0.920,
        bottom=0.18,
        top=0.425,
        wspace=0.15,
        width_ratios=[1.0, 1.0],
    )
    ax_bar = fig.add_subplot(lower[0, 0])
    by_age["Age_Group"] = pd.Categorical(by_age["Age_Group"], AGE_ORDER, ordered=True)
    by_age = by_age.sort_values("Age_Group")
    rate = pd.to_numeric(by_age["codetection_rate_pct"], errors="coerce").to_numpy(float)
    ax_bar.bar(np.arange(len(AGE_ORDER)), rate, color="#df8920", width=0.78)
    ax_bar.set_xticks(
        np.arange(len(AGE_ORDER)),
        [_display_age(value) for value in AGE_ORDER],
        rotation=45,
        ha="right",
    )
    ymax = max(5.0, math.ceil(float(np.nanmax(rate)) * 1.1 / 5.0) * 5.0)
    ax_bar.set_ylim(0, ymax)
    ax_bar.set_ylabel("Rate (% of total events)")
    ax_bar.set_xlabel("Age group")
    ax_bar.set_title("Phase II (2019–2025): co-detection rate by age", pad=2, fontsize=7)
    ax_bar.grid(True, axis="y", linestyle="--", alpha=0.6)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.text(-0.09, 1.04, "J", transform=ax_bar.transAxes, ha="left", va="bottom", fontsize=12, fontweight="bold")

    ax_matrix = fig.add_subplot(lower[0, 1])
    both_positive = matrix.pivot(index="virus_a", columns="virus_b", values="both_positive").reindex(index=VIRUSES_PHASE2, columns=VIRUSES_PHASE2)
    values = both_positive.reindex(index=VIRUSES_PHASE2, columns=VIRUSES_PHASE2).to_numpy(float)
    np.fill_diagonal(values, 0.0)
    mesh = ax_matrix.pcolormesh(values, cmap="Oranges", vmin=0, vmax=float(np.nanmax(values)), edgecolors="none")
    ax_matrix.set_xticks(
        np.arange(len(VIRUSES_PHASE2)) + 0.5,
        [DISPLAY_VIRUS[v] for v in VIRUSES_PHASE2],
        rotation=45,
        ha="right",
    )
    ax_matrix.set_yticks(np.arange(len(VIRUSES_PHASE2)) + 0.5, [DISPLAY_VIRUS[v] for v in VIRUSES_PHASE2])
    for row in range(len(VIRUSES_PHASE2)):
        for col in range(len(VIRUSES_PHASE2)):
            value = values[row, col]
            colour = "white" if value >= 0.60 * float(np.nanmax(values)) else "black"
            ax_matrix.text(col + 0.5, row + 0.5, f"{int(value)}", ha="center", va="center", fontsize=7, color=colour)
    ax_matrix.invert_yaxis()
    ax_matrix.set_title("Phase II: co-detection count matrix", pad=2, fontsize=7)
    colorbar = fig.colorbar(mesh, ax=ax_matrix, fraction=0.045, pad=0.045)
    colorbar.solids.set_rasterized(False)
    colorbar.set_label("Co-detected events, n", fontsize=7)
    colorbar.ax.tick_params(labelsize=7)
    ax_matrix.text(-0.09, 1.04, "K", transform=ax_matrix.transAxes, ha="left", va="bottom", fontsize=12, fontweight="bold")

    fig.text(0.015, 0.963, "Fig 3", ha="left", va="top", fontsize=18)
    fig.text(
        0.52,
        0.915,
        "2019–2025: Monthly Co-detection Rate Trends by Age Group",
        ha="center",
        va="center",
        fontsize=10,
    )
    fig._preserve_fixed_canvas = True
    return fig


def _plot_figure_3(branch: str, tables: dict[str, pd.DataFrame]):
    monthly = tables["Fig3_monthly_age_codetection.csv"].loc[lambda x: x["branch"] == branch].copy()
    fig = _plot_figure_3_composition(branch, tables)
    fig.set_size_inches(7, 7.5)
    # Reference composition: nine purple trend panels above two orange panels.
    for text in list(fig.texts):
        text.remove()
    upper = fig.add_gridspec(3, 3, left=.075, right=.984, bottom=.460, top=.969, hspace=.24, wspace=.14)
    for index, ax in enumerate(fig.axes[:9]):
        ax.set_subplotspec(upper[index // 3, index % 3])
        ax.set_ylim(-5, 105)
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        ax.tick_params(axis="y", labelleft=index % 3 == 0)
        ax.set_ylabel("Co-detection rate (%)" if index % 3 == 0 else "")
        ax.texts[0].set_position((-.15, 1.015))
        ax.texts[0].set_fontsize(11)
        for patch in ax.patches:
            patch.set_facecolor("#f2f2f2")
            patch.set_edgecolor("none")
            patch.set_alpha(1)
        for line in ax.get_xgridlines() + ax.get_ygridlines():
            line.set_color("#ededed")
            line.set_alpha(1)
        ax.set_xlim(pd.Timestamp("2019-01-01"), pd.Timestamp("2026-03-01"))
        expected = monthly.loc[monthly.Age_Group.eq(AGE_ORDER[index])].sort_values(["Year", "Month"])
        np.testing.assert_allclose(ax.lines[0].get_ydata(), expected.rate_display_pct)
        if branch == "minimal_corrected":
            assert np.isfinite(ax.lines[0].get_ydata()).all()
    # Equal-width data panels with dedicated room for the matrix colorbar.
    bar, matrix, colorbar = fig.axes[9:12]
    bar.set_position([.075, .10, .365, .280])
    matrix.set_position([.548, .10, .363, .280])
    colorbar.set_position([.925, .10, .012, .280])
    bar.set_title("Co-detection rate by age", fontsize=8, pad=5)
    matrix.set_title("Co-detection count matrix", fontsize=8, pad=5)
    for label in bar.get_xticklabels():
        label.set_rotation(45)
        label.set_ha("right")
    for label in matrix.get_xticklabels():
        label.set_rotation(0)
        label.set_ha("center")
    for line in bar.get_ygridlines():
        line.set_alpha(1)
        line.set_color("#ededed")
    colorbar.set_ylabel("")
    return fig


def _plot_figure_4(branch: str, tables: dict[str, pd.DataFrame]):
    plt, mdates, _, _ = _activate_matplotlib()
    predictions = tables["Fig4_predictions.csv"].copy()
    predictions["selected_for_plot"] = _as_bool(predictions["selected_for_plot"])
    predictions["date"] = pd.to_datetime(predictions["date"], errors="raise")
    predictions = predictions.loc[(predictions["branch"] == branch) & predictions["selected_for_plot"]]
    importance = tables["Fig4_selected_feature_importance.csv"].loc[lambda x: x["branch"] == branch].copy()
    performance = tables["Fig4_selected_performance.csv"].copy()
    s5 = tables["FigS5_rf_lag_metrics.csv"].copy()
    s5["legacy_test_selected"] = _as_bool(s5["legacy_test_selected"])

    fig, axes = plt.subplots(6, 3, figsize=FIGURE_SIZES_INCH["Figure_4"])
    line_handles = []
    for index, virus in enumerate(VIRUSES_PHASE2):
        ax = axes[index // 3, index % 3]
        subset = predictions.loc[predictions["virus"] == virus].sort_values("date")
        lag = int(pd.to_numeric(subset["lag"], errors="raise").iloc[0])
        actual = ax.plot(subset["date"], pd.to_numeric(subset["actual"], errors="coerce"), color="#4a4a4a", linewidth=1.25, label="Actual")[0]
        rf = ax.plot(subset["date"], pd.to_numeric(subset["predicted"], errors="coerce"), color="#2952cc", linestyle="--", linewidth=1.05, label="RF")[0]
        persistence = ax.plot(subset["date"], pd.to_numeric(subset["persistence"], errors="coerce"), color="#3b8f45", linestyle=":", linewidth=0.95, label="Persistence")[0]
        seasonal = ax.plot(subset["date"], pd.to_numeric(subset["seasonal_naive"], errors="coerce"), color="#e07a18", linestyle="-.", linewidth=0.95, label="Seasonal naive")[0]
        if index == 0:
            line_handles = [actual, rf, persistence, seasonal]
        ax.set_title(f"{DISPLAY_VIRUS[virus]} (lag {lag} week{'s' if lag != 1 else ''})", pad=1.5)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("Positive rate (%)" if index % 3 == 0 else "")
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 3, 5, 7, 9, 11]))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        _style_line_axis(ax, show_y=index % 3 == 0, show_x=True)
        _panel_label(ax, string.ascii_uppercase[index])

    for offset, virus in enumerate(VIRUSES_PHASE2):
        index = offset + 9
        ax = axes[3 + offset // 3, offset % 3]
        subset = importance.loc[importance["virus"] == virus].copy()
        subset["rank"] = pd.to_numeric(subset["rank"], errors="coerce")
        subset["importance"] = pd.to_numeric(subset["importance"], errors="coerce")
        subset = subset.sort_values("rank").head(5).sort_values("importance", ascending=True)
        colours = plt.get_cmap("viridis")(np.linspace(0.72, 0.18, len(subset)))
        ax.barh(np.arange(len(subset)), subset["importance"], color=colours)
        ax.set_yticks(np.arange(len(subset)), [str(value).replace("_", " ") for value in subset["feature"]])
        ax.set_xlabel("Training-period importance")
        ax.set_title(f"{DISPLAY_VIRUS[virus]}: training-only RF importance", pad=1.5, fontsize=7)
        ax.grid(True, axis="x")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if branch == "minimal_corrected":
            row = performance.loc[performance["virus"] == virus].iloc[0]
            baseline_name = str(row["best_baseline"]).replace("_", " ")
            baseline_short = "PERS" if baseline_name == "Persistence" else "SNAIVE"
            note = (
                f"RF: {float(row['rf_r2_common']):+.2f}/{float(row['rf_rmse_common']):.2f}\n"
                f"{baseline_short}: {float(row['baseline_r2_common']):+.2f}/{float(row['baseline_rmse_common']):.2f}\n"
                f"Δ: {float(row['delta_rmse_rf_minus_baseline']):+.2f} "
                f"[{float(row['ci95_low']):+.2f},{float(row['ci95_high']):+.2f}]"
            )
        else:
            row = s5.loc[(s5["branch"] == branch) & (s5["virus"] == virus) & s5["legacy_test_selected"]].iloc[0]
            baseline_name = str(row["best_baseline"]).replace("_", " ")
            baseline_short = "PERS" if baseline_name == "Persistence" else "SNAIVE"
            note = (
                f"RF: {float(row['rf_r2_common']):+.2f}/{float(row['rf_rmse_common']):.2f}\n"
                f"{baseline_short}: {float(row['best_baseline_r2_common']):+.2f}/{float(row['best_baseline_rmse_common']):.2f}\n"
                "test-selected reference"
            )
        ax.set_ylim(-2.8, len(subset) - 0.35)
        ax.text(0.985, 0.018, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="#333333", linespacing=0.92)
        ax.text(
            -0.14,
            1.035,
            string.ascii_uppercase[index],
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=9.2,
            fontweight="bold",
        )

    fig.legend(line_handles, [handle.get_label() for handle in line_handles], loc="upper center", bbox_to_anchor=(0.5, 0.972), ncol=4, frameon=False)
    if branch == "minimal_corrected":
        footer = (
            "Lag fixed by 2023–2024 training CV; panel pairs are R²/RMSE.\n"
            "Δ is RMSE(RF−best naive) [paired 4-week block-bootstrap 95% CI].\n"
            "No selected RF significantly outperformed its best baseline.\n"
            "RF results are exploratory and not biologically interpreted."
        )
    else:
        footer = (
            "Legacy reference (internal; current-environment rerun, not an exact numerical reproduction).\n"
            "test-selected lags are exploratory/descriptive only, are not biologically interpreted, and are not used for revised claims."
        )
    fig.suptitle("Observed and predicted 2025 weekly positivity with RF diagnostics", y=0.997, fontsize=8.8)
    fig.text(0.5, 0.006, footer, ha="center", fontsize=7, multialignment="center")
    fig.subplots_adjust(left=0.115, right=0.945, top=0.925, bottom=0.115, hspace=0.80, wspace=0.48)
    return fig


def _heatmap_monthly(ax, frame: pd.DataFrame, viruses: Iterable[str], vmax: float, title: str, plt):
    frame = frame.copy()
    frame["date"] = _dates(frame)
    pivot = frame.pivot(index="virus", columns="date", values="rate_display_pct").reindex(index=list(viruses))
    values = pivot.to_numpy(float)
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad("#eeeeee")
    mesh = ax.pcolormesh(np.ma.masked_invalid(values), cmap=cmap, vmin=0, vmax=vmax, edgecolors="none")
    positions = np.arange(0, len(pivot.columns), 12)
    ax.set_xticks(positions + 0.5, [pivot.columns[i].strftime("%Y-%m") for i in positions], rotation=40, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)) + 0.5, [DISPLAY_VIRUS[value] for value in pivot.index])
    ax.invert_yaxis()
    ax.set_title(title, pad=2)
    ax.set_xlabel("Month")
    ax.set_ylabel("Pathogen")
    return mesh


def _plot_figure_s1(branch: str, tables: dict[str, pd.DataFrame]):
    plt, mdates, _, _ = _activate_matplotlib()
    data = tables["FigS1_monthly_rates.csv"].loc[lambda x: x["branch"] == branch].copy()
    fig = plt.figure(figsize=FIGURE_SIZES_INCH["Figure_S1"])
    grid = fig.add_gridspec(5, 6, height_ratios=[0.92, 0.92, 1, 1, 1], hspace=0.72, wspace=0.55)

    ax_a = fig.add_subplot(grid[0, :])
    mesh_a = _heatmap_monthly(ax_a, data.loc[data["phase"] == "phase1"], VIRUSES_PHASE1, 40, "Monthly positive rate, 2013–2018 (IgM IFA)", plt)
    cb_a = fig.colorbar(mesh_a, ax=ax_a, fraction=0.018, pad=0.016)
    cb_a.solids.set_rasterized(False)
    cb_a.set_label("Positive rate (%)", fontsize=7)
    cb_a.ax.tick_params(labelsize=7)
    _panel_label(ax_a, "A")

    ax_b = fig.add_subplot(grid[1, :])
    mesh_b = _heatmap_monthly(ax_b, data.loc[data["phase"] == "phase2"], VIRUSES_PHASE2, 60, "Monthly positive rate, 2019–2025 (PCR)", plt)
    cb_b = fig.colorbar(mesh_b, ax=ax_b, fraction=0.018, pad=0.016)
    cb_b.solids.set_rasterized(False)
    cb_b.set_label("Positive rate (%)", fontsize=7)
    cb_b.ax.tick_params(labelsize=7)
    _panel_label(ax_b, "B")

    phase1 = data.loc[data["phase"] == "phase1"].copy()
    phase1["date"] = _dates(phase1)
    for index, virus in enumerate(VIRUSES_PHASE1):
        row = 2 + index // 3
        col = (index % 3) * 2
        ax = fig.add_subplot(grid[row, col : col + 2])
        subset = phase1.loc[phase1["virus"] == virus].sort_values("date")
        ax.plot(subset["date"], pd.to_numeric(subset["rate_display_pct"], errors="coerce"), color=VIRUS_COLORS[virus])
        ax.set_title(DISPLAY_VIRUS[virus], pad=1.5)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("Positive rate (%)" if index % 3 == 0 else "")
        ax.xaxis.set_major_locator(mdates.YearLocator(2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.tick_params(axis="x", rotation=35)
        _style_line_axis(ax, show_y=index % 3 == 0, show_x=index >= 5)
        _panel_label(ax, string.ascii_uppercase[index + 2])
    blank = fig.add_subplot(grid[4, 4:6])
    blank.axis("off")

    fig.suptitle("Monthly pathogen positivity across the two diagnostic phases", y=0.995, fontsize=9)
    fig.text(0.5, 0.008, "Grey heatmap cells and line gaps denote months with fewer than 20 tests.", ha="center", fontsize=7)
    fig.subplots_adjust(left=0.09, right=0.95, top=0.96, bottom=0.07)
    return fig


def _plot_environment(branch: str, filename: str, period: str, tables: dict[str, pd.DataFrame], figure_name: str):
    plt, _, _, _ = _activate_matplotlib()
    data = tables[filename].loc[lambda x: x["branch"] == branch].copy()
    data["fdr05_period"] = _as_bool(data["fdr05_period"])
    fig, axes = plt.subplots(3, 3, figsize=FIGURE_SIZES_INCH[figure_name])
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e2e2e2")
    mesh = None
    for index, (ax, age) in enumerate(zip(axes.flat, AGE_ORDER, strict=True)):
        subset = data.loc[data["age_group"] == age]
        rho = subset.pivot(index="virus", columns="environment", values="rho_display").reindex(index=VIRUSES_PHASE2, columns=ENVIRONMENT_ORDER)
        q = subset.pivot(index="virus", columns="environment", values="q_period").reindex(index=VIRUSES_PHASE2, columns=ENVIRONMENT_ORDER)
        n = subset.pivot(index="virus", columns="environment", values="n_months").reindex(index=VIRUSES_PHASE2, columns=ENVIRONMENT_ORDER)
        values = rho.to_numpy(float)
        mesh = ax.pcolormesh(np.ma.masked_invalid(values), cmap=cmap, vmin=-1, vmax=1, edgecolors="white", linewidth=0.25)
        for row in range(len(VIRUSES_PHASE2)):
            for col in range(len(ENVIRONMENT_ORDER)):
                value = values[row, col]
                if np.isfinite(value) and np.isfinite(q.iloc[row, col]) and float(q.iloc[row, col]) <= 0.05:
                    colour = "white" if abs(value) >= 0.58 else "black"
                    ax.text(col + 0.5, row + 0.5, "*", ha="center", va="center", fontsize=7, fontweight="bold", color=colour)
        ax.set_xticks(np.arange(len(ENVIRONMENT_ORDER)) + 0.5, [ENVIRONMENT_LABELS[value] for value in ENVIRONMENT_ORDER], rotation=60, ha="right")
        ax.set_yticks(np.arange(len(VIRUSES_PHASE2)) + 0.5, [DISPLAY_VIRUS[value] for value in VIRUSES_PHASE2])
        ax.invert_yaxis()
        finite_n = pd.to_numeric(n.where(rho.notna()).stack(future_stack=True), errors="coerce").dropna()
        n_text = "no cells" if finite_n.empty else f"n={int(finite_n.min())}–{int(finite_n.max())}"
        ax.set_title(f"Age {_display_age(age)}; {n_text}", pad=1.5, fontsize=7)
        _panel_label(ax, string.ascii_uppercase[index])
    # Apply final subplot geometry before adding the shared colour bar.  Adding
    # the colour bar first and then resetting the axes made it cover the AQI
    # column in panel F.
    fig.subplots_adjust(left=0.095, right=0.90, top=0.945, bottom=0.13, hspace=0.68, wspace=0.48)
    if mesh is not None:
        colorbar = fig.colorbar(mesh, ax=list(axes.flat), fraction=0.018, pad=0.018)
        colorbar.solids.set_rasterized(False)
        colorbar.set_label("Spearman ρ", fontsize=7)
        colorbar.ax.tick_params(labelsize=7)
    title_suffix = " (exploratory partial year)" if period == "2019" else ""
    fig.suptitle(f"Age-stratified environmental Spearman correlations: {period}{title_suffix}", y=0.995, fontsize=8.5)
    fig.text(
        0.5,
        0.008,
        "Monthly pathogen cells with tested n<20 were excluded; paired months n<6 are blank; * BH-FDR q<0.05 within period.\nExact ρ, n, p and q values are provided in the supplementary tables.",
        ha="center",
        fontsize=7,
    )
    return fig


def _plot_figure_s5(branch: str, tables: dict[str, pd.DataFrame]):
    plt, _, _, patches = _activate_matplotlib()
    data = tables["FigS5_rf_lag_metrics.csv"].loc[lambda x: x["branch"] == branch].copy()
    data["selected_by_training_cv"] = _as_bool(data["selected_by_training_cv"])
    data["legacy_test_selected"] = _as_bool(data["legacy_test_selected"])
    r2 = data.pivot(index="virus", columns="lag", values="rf_r2_common").reindex(index=S5_VIRUS_ORDER, columns=[0, 1, 2, 3, 4])
    delta = data.pivot(index="virus", columns="lag", values="delta_r2_rf_minus_best").reindex(index=S5_VIRUS_ORDER, columns=[0, 1, 2, 3, 4])
    values = r2.to_numpy(float)
    vmin = min(-1.0, float(np.nanmin(values)))
    vmax = max(1.0, float(np.nanmax(values)))

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCH["Figure_S5"])
    mesh = ax.pcolormesh(np.ma.masked_invalid(values), cmap="RdYlGn", vmin=vmin, vmax=vmax, edgecolors="white", linewidth=0.55)
    for row, virus in enumerate(S5_VIRUS_ORDER):
        for col, lag in enumerate([0, 1, 2, 3, 4]):
            r2_value = r2.loc[virus, lag]
            delta_value = delta.loc[virus, lag]
            if np.isfinite(r2_value):
                ax.text(col + 0.5, row + 0.5, f"R² {r2_value:+.2f}\nΔ {delta_value:+.2f}", ha="center", va="center", fontsize=7)
            selected = data.loc[(data["virus"] == virus) & (pd.to_numeric(data["lag"], errors="coerce") == lag)]
            if branch == "minimal_corrected" and _as_bool(selected["selected_by_training_cv"]).iloc[0]:
                ax.add_patch(patches.Rectangle((col + 0.04, row + 0.04), 0.92, 0.92, fill=False, edgecolor="black", linewidth=1.25))
            if branch == "legacy_exact" and _as_bool(selected["legacy_test_selected"]).iloc[0]:
                ax.add_patch(patches.Rectangle((col + 0.04, row + 0.04), 0.92, 0.92, fill=False, edgecolor="black", linewidth=1.05, linestyle="--"))
    ax.set_xticks(np.arange(5) + 0.5, ["0\nnowcast", "1", "2", "3", "4"])
    ax.set_yticks(np.arange(len(S5_VIRUS_ORDER)) + 0.5, [DISPLAY_VIRUS[value] for value in S5_VIRUS_ORDER])
    ax.invert_yaxis()
    ax.set_xlabel("Forecast lag (weeks)")
    ax.set_ylabel("Pathogen")
    ax.set_title("RF test R² (colour) and ΔR² versus best naive baseline (cell text)", pad=4, fontsize=8)
    colorbar = fig.colorbar(mesh, ax=ax, fraction=0.032, pad=0.025)
    colorbar.solids.set_rasterized(False)
    colorbar.set_label("RF test R²", fontsize=7)
    colorbar.ax.tick_params(labelsize=7)
    selection_text = (
        "Solid outline: lag selected using 2023–2024 training CV; ΔR²>0 favours RF."
        if branch == "minimal_corrected"
        else "Dashed outline: legacy 2025-test-selected lag (internal current-environment reference only); ΔR²>0 favours RF."
    )
    fig.text(0.5, 0.012, selection_text, ha="center", fontsize=7)
    fig.subplots_adjust(left=0.13, right=0.94, top=0.91, bottom=0.17)
    return fig


def build_figure(figure_name: str, branch: str, tables: dict[str, pd.DataFrame]):
    if figure_name == "Figure_1":
        return _plot_figure_1(branch, tables)
    if figure_name == "Figure_2":
        return _plot_figure_2(branch, tables)
    if figure_name == "Figure_3":
        return _plot_figure_3(branch, tables)
    if figure_name == "Figure_4":
        return _plot_figure_4(branch, tables)
    if figure_name == "Figure_S1":
        return _plot_figure_s1(branch, tables)
    if figure_name == "Figure_S2":
        return _plot_environment(branch, "FigS2_environment_cells.csv", "2019", tables, figure_name)
    if figure_name == "Figure_S3":
        return _plot_environment(branch, "FigS3_environment_cells.csv", "2020–2022", tables, figure_name)
    if figure_name == "Figure_S4":
        return _plot_environment(branch, "FigS4_environment_cells.csv", "2023–2025", tables, figure_name)
    if figure_name == "Figure_S5":
        return _plot_figure_s5(branch, tables)
    raise KeyError(figure_name)


def _add_branch_page_label(fig, branch: str) -> None:
    """Make every internal legacy page self-identifying if detached from its folder."""

    if branch == "legacy_exact":
        fig.text(
            0.004,
            0.5,
            "Legacy reference (internal; current-environment rerun; not numerically exact; do not submit)",
            ha="left",
            va="center",
            rotation=90,
            fontsize=7,
            color="#8a1c1c",
        )


def _export_figure(fig, stem: Path) -> list[Path]:
    from PIL import Image

    stem.parent.mkdir(parents=True, exist_ok=True)
    outputs = [stem.with_suffix(ext) for ext in (".svg", ".eps", ".pdf", ".tiff")]
    preserve_canvas = bool(getattr(fig, "_preserve_fixed_canvas", False))
    save_kwargs = {"transparent": False}
    if not preserve_canvas:
        save_kwargs.update({"bbox_inches": "tight", "pad_inches": 0.02})
    for path in outputs[:3]:
        fig.savefig(path, format=path.suffix.lstrip("."), **save_kwargs)
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=600, **save_kwargs)
    buffer.seek(0)
    with Image.open(buffer) as raster:
        raster.convert("RGB").save(outputs[3], format="TIFF", compression="tiff_lzw", dpi=(600, 600))
    return outputs


def _pdf_inches(path: Path) -> tuple[float, float] | None:
    # Matplotlib writes an uncompressed MediaBox dictionary near the page object.
    import re

    match = re.search(rb"/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]", path.read_bytes())
    if not match:
        return None
    x0, y0, x1, y1 = (float(value) for value in match.groups())
    return (x1 - x0) / 72.0, (y1 - y0) / 72.0


def verify_outputs(paths: list[Path]) -> dict[str, object]:
    from PIL import Image

    by_suffix = {path.suffix.lower(): path for path in paths}
    if set(by_suffix) != {".svg", ".eps", ".pdf", ".tiff"}:
        raise ValueError(f"Incomplete format set: {sorted(by_suffix)}")
    svg = by_suffix[".svg"].read_text(encoding="utf-8")
    if "<text" not in svg:
        raise ValueError(f"SVG contains no editable text: {by_suffix['.svg']}")
    if "<image" in svg:
        raise ValueError(f"SVG contains an embedded raster image: {by_suffix['.svg']}")
    eps = by_suffix[".eps"].read_bytes().lower()
    if b"colorimage" in eps:
        raise ValueError(f"EPS contains a raster colour image: {by_suffix['.eps']}")
    pdf_bytes = by_suffix[".pdf"].read_bytes()
    if b"/Subtype /Image" in pdf_bytes:
        raise ValueError(f"PDF contains a raster image XObject: {by_suffix['.pdf']}")
    pdf_size = _pdf_inches(by_suffix[".pdf"])
    if pdf_size is not None and (pdf_size[0] > 7.01 or pdf_size[1] > 9.01):
        raise ValueError(f"PDF exceeds 7x9 inches: {pdf_size} at {by_suffix['.pdf']}")
    with Image.open(by_suffix[".tiff"]) as tiff:
        dpi = tiff.info.get("dpi", (0, 0))
        if tiff.mode != "RGB":
            raise ValueError(f"TIFF mode is {tiff.mode}, expected RGB")
        if not all(abs(float(value) - 600) <= 1 for value in dpi):
            raise ValueError(f"TIFF dpi is {dpi}, expected 600")
        width_inches = tiff.width / float(dpi[0])
        height_inches = tiff.height / float(dpi[1])
        if width_inches > 7.01 or height_inches > 9.01:
            raise ValueError(f"TIFF exceeds 7x9 inches: {(width_inches, height_inches)}")
        compression = str(tiff.info.get("compression", "")).lower()
        if "tiff_lzw" not in compression and compression != "5":
            raise ValueError(f"TIFF compression is {compression}, expected LZW")
    return {
        "pdf_inches": pdf_size,
        "tiff_dpi": [float(value) for value in dpi],
        "tiff_mode": "RGB",
        "vector_text": True,
        "vector_raster_layers": False,
    }


def render_all(
    tables: dict[str, pd.DataFrame],
    output_dir: Path,
    branches: Iterable[str],
    figures: Iterable[str] = FIGURE_ORDER,
) -> dict[str, object]:
    plt, _, _, _ = _activate_matplotlib()
    rendered: list[dict[str, object]] = []
    for branch in branches:
        branch_dir = output_dir / OUTPUT_SET[branch]
        for figure_name in figures:
            figure = build_figure(figure_name, branch, tables)
            _add_branch_page_label(figure, branch)
            paths = _export_figure(figure, branch_dir / figure_name)
            plt.close(figure)
            qa = verify_outputs(paths)
            rendered.append(
                {
                    "branch": branch,
                    "output_set": OUTPUT_SET[branch],
                    "figure": figure_name,
                    "files": [path.relative_to(output_dir).as_posix() for path in paths],
                    "qa": qa,
                }
            )
    return {"status": "PASS", "rendered": rendered, "pdf_count": sum(1 for row in rendered for path in row["files"] if path.endswith(".pdf"))}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--branch", choices=("all",) + BRANCHES, default="all")
    parser.add_argument(
        "--figure",
        action="append",
        choices=FIGURE_ORDER,
        default=None,
        help="Optional repeatable subset selector; omit to render all nine figures.",
    )
    parser.add_argument("--render", action="store_true", help="Create the four figure formats after the PDF marker.")
    parser.add_argument(
        "--pdf-marker-confirmed",
        action="store_true",
        help="Required with --render. Supply only after the parent workflow has completed the PDF marker.",
    )
    parser.add_argument("--qa-json", type=Path, default=None, help="Optional QA JSON path; written only in render mode.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    branches = BRANCHES if args.branch == "all" else (args.branch,)
    tables = load_inputs(args.input_dir, branches)
    static_qa = validate_inputs(tables, branches)
    if not args.render:
        if args.pdf_marker_confirmed:
            raise SystemExit("--pdf-marker-confirmed has no effect without --render")
        print(json.dumps(static_qa, indent=2, ensure_ascii=False))
        return 0
    if not args.pdf_marker_confirmed:
        raise SystemExit("Rendering blocked: run the parent PDF marker, then add --pdf-marker-confirmed.")
    figures = FIGURE_ORDER if args.figure is None else tuple(dict.fromkeys(args.figure))
    result = {"static_qa": static_qa, "output_qa": render_all(tables, args.output_dir, branches, figures)}
    if args.qa_json is not None:
        args.qa_json.parent.mkdir(parents=True, exist_ok=True)
        args.qa_json.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
