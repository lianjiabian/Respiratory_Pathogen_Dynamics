# Respiratory Pathogen Dynamics — public aggregate data

Repository target: <https://github.com/lianjiabian/Respiratory_Pathogen_Dynamics>

This directory is the **public-release candidate**. It contains 11 aggregate
supplementary tables (`Table S1`–`Table S11`) as both CSV and RDA files and the
minimum aggregate inputs needed to redraw Figures 1–4 and S1–S5. It does not
contain patient/sample identifiers, stable patient tokens, exact specimen
dates, exact ages, or event-level pathogen records.

## Contents

- `data/`: one CSV and one single-data-frame RDA per logical supplementary table.
- `figure_inputs/`: final `minimal_corrected` monthly/weekly aggregate plotting
  inputs with unused event-count columns removed.
- `docs/data_dictionary.csv`: object- and variable-level definitions.
- `docs/supplementary_table_map.csv`: mapping from Table S1–S11 to source output.
- `qa/privacy_scan.csv`: identifier-pattern, exact-date, small-cell, and
  complementary-disclosure checks.
- `qa/complementary_suppression_log.csv`: four additional cells hidden to prevent
  a single primary-suppressed cell from being recovered by subtraction.
- `release_manifest.csv`: source/output row counts and SHA-256 values.
- `qa/figure_inputs_manifest.csv`: hashes, dimensions, temporal granularity, and
  privacy classification of every plotting input.
- `scripts/reproduce_public.py`: one-command table validation and nine-figure
  reproduction.
- `scripts/07_render_figures.py`: the final conservative renderer.
- `requirements.txt`: pinned full public-reproduction environment.

`NA` means unavailable, non-estimable, or disclosure-controlled according to
the table's `suppressed` or `status` field. It must never be interpreted as a
negative laboratory result.

The fields ending in `_pct` are numeric **percentage points** (for example,
`0.45` means `0.45%`), not proportions. This convention is retained from the
locked reporting outputs.

## Read an RDA file

Python:

```python
import pyreadr

objects = pyreadr.read_r("data/Table_S1_cohort_pathogen_summary.rda")
cohort_pathogen_summary = objects["cohort_pathogen_summary"]
```

R:

```r
e <- new.env(parent = emptyenv())
load("data/Table_S1_cohort_pathogen_summary.rda", envir = e)
e$cohort_pathogen_summary
```

## Reproduce all reported tables and figures

Create an isolated environment, install the pinned dependencies, and run:

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python scripts\reproduce_public.py --root . --output-dir reproduced_outputs
```

This validates the hashes, object names, dimensions, and columns of S1–S11;
creates seven compact table-check summaries; and redraws all nine figures. Each
figure is written as editable SVG/EPS/PDF plus 600-dpi RGB/LZW TIFF. The run
fails unless exactly 36 figure artifacts pass format, vector-text, raster-layer,
physical-size, and TIFF-resolution checks. Add `--rscript PATH_TO_RSCRIPT` for
an independent R `load()` check of all 11 RDA files.

The weekly `date` field in the Figure 4 input is a Sunday week anchor for an
aggregate time series, not an individual specimen date.

## Rebuild the public data artifacts

The builder is Python-only. Supply a directory containing the 11 locked public
aggregate CSV outputs:

```powershell
python scripts/build_public_release.py --source-dir PATH_TO_PUBLIC_AGGREGATES --output-dir .
python scripts/validate_rda_release.py --root . --manifest release_manifest.csv --output qa/rda_compatibility_validation.csv
```

The second command always performs a Python/pyreadr round trip. Add
`--rscript PATH_TO_RSCRIPT` for an independent R `load()` compatibility check.

## Reproducibility boundary

These aggregate files support verification/regeneration of the 11 supplementary
tables and exact redrawing of the nine reported figures without public release
of clinical event records. They do not reproduce raw-record cleaning, age
parsing, target interpretation, duplicate handling, or same-day event
construction. Those upstream correction steps require the controlled source
workbooks, the controlled RDA layer, and institutional data access.

To preserve the author-designated Figure 3 unchanged, its minimized plotting
input retains the exact aggregate pair counts printed in the figure, including
28 directed cells below 20 (14 reciprocal unordered pairs). This is not
event-level data, but it is an explicit exception to the table-level small-cell
suppression policy. Authors and the institutional data custodian must approve
that disclosure before any GitHub upload; otherwise Figure 3 and its input must
be revised together.

Before the first GitHub push, scan the full staged repository and Git history;
never stage the sibling controlled-data directory.

## Figure 3 continuous monthly curves

Figure 3A–I retains all 726 observed age-by-month estimates, including 213
with fewer than 20 eligible events. Lines connect consecutive observed months;
no missing rates are imputed, interpolated, smoothed, or filled with zero.
The monthly figure input exposes only year, month, age group, displayed rate,
and analysis branch, without event identifiers or individual dates. These
descriptive low-denominator estimates are unstable; see the manuscript legend.
The table-specific suppression rules for Tables S1–S11 remain in place.
The displayed monthly rates and the pair counts are explicit Figure 3 exceptions
to table-level suppression and should be considered together in the existing
author/institutional disclosure review before public release.
