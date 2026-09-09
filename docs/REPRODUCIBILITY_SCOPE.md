# Reproducibility scope of the public release

This public release contains the 11 disclosure-controlled aggregate tables used
for Supplementary Tables S1-S11. Each table is supplied as CSV and as a
single-data-frame RDA file. It also contains minimized monthly/weekly aggregate
inputs for exact redrawing of Figures 1-4 and S1-S5. The files support
independent checking of reported aggregate values, regeneration of the Word
supplementary tables, and nine-figure reproduction.

The public release deliberately excludes event-level clinical records,
patient/sample identifiers, pseudonyms, exact specimen dates, exact ages, and
source clinical workbooks. It therefore does **not** reproduce the upstream
record-cleaning and same-day event-construction stage. That stage requires
controlled institutional access to the source workbooks and the controlled
analysis environment.

The source workbook titled `Supplementary Table 1-3.xlsx` is a legacy derived
compilation rather than a raw-data source for the current minimal-corrected
analysis. Its clinical sheets correspond to earlier derived CSV inputs. It
cannot by itself reproduce current corrections to event construction, age
parsing, pooled parainfluenza handling, or the Phase II *M. pneumoniae* rule.

Small counts in the public supplementary tables are primary-suppressed.
Secondary suppression is also applied where a primary-suppressed value could
otherwise be recovered by subtraction. See
`qa/complementary_suppression_log.csv`. Suppressed values must not be
reverse-engineered or replaced in a public derivative.

The author-designated Figure 3 displays exact aggregate pair counts, including
small cells. Its minimized plotting input therefore retains those values. This
exception requires explicit author and institutional data-custodian approval
before public upload; if approval is not available, the figure and its plotting
input must be revised together.

The immutable build evidence is recorded in `release_manifest.csv`,
`qa/rda_compatibility_validation.csv`, and
`qa/artifact_sha256_manifest.csv`.
