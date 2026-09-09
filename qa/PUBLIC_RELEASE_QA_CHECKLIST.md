# Public release QA checklist

- [x] 11 logical supplementary tables mapped continuously as S1-S11.
- [x] One CSV and one single-data-frame RDA created for every table.
- [x] All 11 RDA files round-tripped with Python/pyreadr.
- [x] All 11 RDA files loaded independently with R.
- [x] No direct identifier, stable pseudonym, exact specimen date, exact age,
  or event-level clinical row is present.
- [x] Primary small-cell suppression retained.
- [x] Secondary suppression applied where subtraction could reveal a primary
  suppressed value.
- [x] Data dictionary and table mapping included.
- [x] Figure inputs limited to the final `minimal_corrected` branch and to
  month/week-level aggregate fields consumed by the renderer.
- [x] One-command public reproduction independently rebuilt all 9 figures in
  SVG, EPS, PDF, and 600-dpi RGB/LZW TIFF (36/36 artifacts passed).
- [x] Figure inputs and reproduced artifacts have hash manifests.
- [x] SHA-256 artifact inventory generated.
- [ ] Authors and institutional data custodian approve the exact aggregate
  pair-count cells displayed in unchanged Figure 3 (28 directed cells below
  20, representing 14 reciprocal unordered pairs).
- [ ] Authors confirm the manuscript callouts and final Word table pagination.
- [ ] Authors confirm the institutional controlled-access body and contact.
- [ ] Before upload, scan the complete Git index and history and confirm that
  the controlled sibling directory has never been staged.
