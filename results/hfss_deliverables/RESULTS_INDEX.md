# HFSS results index

Start in `deliverables/`:

- `spectra_tables/`: one-nanometer T/R/A tables for all 12 physical cases.
- `native_exports/`: native HFSS modal data and interpolation provenance.
- `resonance_analysis/`: marked T/R/A figures and `resonance_candidates.csv`.
- `nearfield_maps/`: direct-field XY and XZ maps, normalized as `|E|² / |E₀|²`.

Near-field naming follows `<case>_<angle>deg_<candidate>_<wavelength>_<plane>`.
Every candidate at 0 and 45 degrees has distinct `_xy` and `_xz` PNG/PDF figures and a compressed `.npz` data file. Raw HFSS field grids remain beside the data files for traceability.

`projects/` holds the solved AEDT projects. `diagnostics/` preserves disk-interrupted, insufficient-mode, and superseded-render artifacts. The latter were retained rather than deleted after the XY/XZ filename correction.

## Interpretation boundary

Resonance identifiers are automated absorption-peak candidates. They are not asserted to be mode-tracked physical resonance labels across different incidence angles. The near-fields use direct saved HFSS field solutions at the candidate frequencies.
