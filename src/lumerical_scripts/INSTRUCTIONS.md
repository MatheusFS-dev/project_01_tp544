# Gold nanoring metasurface, Lumerical FDTD

## Purpose

This project implements Tarefa 3 of the supplied assignment, Metassuperfície de nanoanéis circulares de Au. A square-periodic 3D unit cell is illuminated with TM-polarized light at the six assigned angles. Two structures are modeled independently: A, nanorings suspended in air, and B, nanorings resting on a semi-infinite SiO2 approximation.

The simulations are **not run in this package**. You need a licensed Ansys Lumerical FDTD Solutions installation with the required material database. These LSF files have been checked against documentation and structurally inspected, not executed against Lumerical. Confirm property support on your installed version.

## Files

| Script | Role |
|---|---|
| `settings_nanorings.lsf` | Assignment dimensions and adjustable numerical configuration |
| `build_nanorings.lsf` | Rebuild the geometry, FDTD region, source, mesh, and monitors |
| `run_one.lsf` | Generate, save, simulate, and analyze a selected case and angle |
| `run_all.lsf` | Run 2 media × 6 angles, saving each `.fsp` |
| `analyze_nanorings.lsf` | Export spectra, T/R/A plots, and candidate extrema |
| `export_nearfields.lsf` | Load an existing `.fsp` and export electric-field maps at **user-confirmed** resonances |
| `run_reference.lsf` | Reference model without Au for optical-normalization checks |

## Recommended use in FDTD Solutions

1. Extract the entire folder into a local writable directory, open FDTD Solutions **on an empty/disposable model because `build_nanorings` deletes existing model objects**, and switch the current working directory to that folder. LSF scripts invoke each other by filename without the `.lsf` extension.
2. Start with `run_reference.lsf` using `case_index = 1`, then test `case_index = 2`. Check physical plausibility and monitor signs.
3. Open `run_one.lsf`, select `case_index = 1` and `theta_deg = 0`, and execute it. Inspect the generated `.fsp` geometry, mesh, source, and spectra before any large sweep.
4. Perform at least a mesh and PML convergence comparison. Set numerical values in `settings_nanorings.lsf` and rerun the chosen case. Do not overwrite the baseline `.fsp` until its results are backed up.
5. Run `run_all.lsf`. This saves all 12 requested `.fsp` simulations, spectra in `.tsv` and `.mat`, diagnostic plots in `.jpg`, and `_candidates.tsv` lists.
6. Confirm all resonances **manually** by comparing T/R/A spectral features and the associated fields. The candidate-extremum detector is not a mode classifier.
7. For each case at 0 and 45 degrees, edit `export_nearfields.lsf` to choose the appropriate case, angle, and observed resonance wavelengths. Run the script. Repeat for the four combinations. The initially included `[800,1000]` wavelengths are placeholder examples, not measured resonances.
8. Use the generated spectra and field maps in the final report. Keep a record of all numerical configuration settings and convergence outcomes.

For geometry-only inspection, open `run_one.lsf`, execute `settings_nanorings; build_nanorings;` after defining `case_index`, `theta_deg`, and `record_nearfields`, then inspect the CAD layout without launching `run`.

## Exact assignment parameters

| Parameter | Value | Reason |
|---|---:|---|
| `outer_radius` | 90 nm | Specified outer radial extent of the Au ring |
| `inner_radius` | 50 nm | Specified hole radius, making annulus width 40 nm |
| `ring_height` | 30 nm | Specified axial thickness, modeled from z = 0 to +30 nm |
| `period_x`, `period_y` | 320 nm | Specified square-lattice pitch. Exactly one period forms the periodic simulation cell |
| `gold_material` | `Au (Gold) - Palik` | Requested Palik Au dispersive optical constants as used in Lumerical. If the installed material name differs, inspect the database and change it explicitly |
| `substrate_index` | 1.45 | The assignment allows approximate nondispersive SiO2 or an equivalent dispersive model |
| `angles_deg` | 0, 15, 30, 45, 60, 75 degrees | All mandatory angles |
| `min_wavelength`, `max_wavelength` | 450, 1600 nm | Mandatory optical spectral window |
| Source `polarization angle` | 0 degrees | Lumerical p/TM polarization (not the same as the incidence angle) |
| Source `angle phi` | 0 degrees | Places the oblique beam propagation in the xz plane |
| Source `direction` | `Backward` | Illuminates from air above, traveling toward -z |
| Background index | 1 | Air superstrate, and air in the suspended case |
| Substrate upper face | z = 0 | Ensures that the metal ring touches the substrate without a gap |

## Numerical starting values and rationales

These values **are proposals, not imposed by the assignment and not proof of convergence**.

| Parameter | Initial value | Explanation and necessary test |
|---|---:|---|
| `frequency_points` | 401 | Uniform wavelength sampling gives spacing (1600−450)/(401−1) = 2.875 nm. Reasonable first-pass resonance localization, but narrow features need 801 or more points plus long simulation duration |
| `local_mesh_step` | 4 nm | About 10 samples across the ring's 40 nm width and 7.5 across 30 nm thickness. Start here, then compare 3 and 2 nm meshes, particularly the near-field maxima |
| `mesh_accuracy` | 2 | Moderate global automatic mesh for air/substrate far from Au. The explicit 4 nm mesh override controls the sensitive ring region |
| `simulation_duration` | 2500 fs | Upper limit, long enough as an initial optical decay budget for broadband plasmonic FDTD. It must be increased if the solver stops at the maximum time rather than auto shutoff |
| `auto_shutoff` | 1e-6 | Stops after estimated remaining energy is a millionth of the maximum injected level. More stringent than 1e-5 default, and must be tested against 1e-7 where needed |
| `pml_layers` | 16 | Starting absorbing-boundary thickness. For 60 and 75 degrees, increase to 24 or 32 and compare R/T to check grazing-angle PML reflections |
| `pml profile` | 3 | Uses the steep-angle profile in typical Lumerical versions. Confirm this property/option in your version |
| `fdtd_z_min`, `fdtd_z_max` | -900, +850 nm | Creates vertical clearance on both sides of the ring for source, power monitors and PML absorption. Vary z padding to test independence |
| `source_height` | +300 nm | Plane wave launches from air, separated from the ring, with sufficient room below the reflection monitor |
| `reflection_height` | +520 nm | Placed **behind** an incident wave traveling toward -z. Records reflected +z power without needing to subtract incident power |
| `transmission_height` | -500 nm | Records -z propagating power, in the substrate for case B and in air for case A. Lies sufficiently far from the metal |
| `substrate_z_min` | -2000 nm | Extends SiO2 beyond the bottom of the FDTD/PML region, preventing an artificial second interface inside the computational domain |
| `xy_field_height` | +35 nm | 5 nm above the ring top, deliberately in surrounding air for enhanced fields rather than within the dissipative metal |
| `xz_field_z_min`, `xz_field_z_max` | -160, +180 nm | Vertical cut at y = 0 showing the ring, air and substrate for comparison |
| `export_plots` | 1 | Creates JPG figures in the Lumerical GUI. Set to 0 for lower overhead / reduced dependence on graphics |
| Candidate `contrast_threshold` | 0.03 | Flags local A peaks and T dips exceeding a 3-percentage-point contrast over about 29 nm. A plotting convenience, not a scientific threshold |

### Why the solver uses different periodic treatments

At θ = 0°, the source uses `Bloch/periodic` plane wave type and x/y periodic BCs. At θ > 0°, it uses BFAST, whose transverse boundary conditions supersede the conventional periodic settings and hold the physical incidence angle fixed across the bandwidth. The source's default BFAST time-step multiplier is not overridden in this script. Check stability carefully at 60°/75°, where BFAST can be significantly more expensive than normal incidence.

### Why there is no second interface below SiO2

The substrate rectangle reaches below the lower FDTD PML. This models a half-space instead of a thin glass film. Introducing a bottom air boundary in the computational domain would create Fabry–Pérot fringes unrelated to the requested nanoring response.

### Power convention

The source illuminates from +z toward -z. The transmission monitor normal is +z, so its `transmission()` result is negative. We define `T = -transmission("T_monitor")`. The reflection monitor measures upward flux (+z), hence `R = transmission("R_monitor")` is positive. Absorption is `A = 1-T-R`. **Do not take absolute values or clip negative A to zero**: doing so can hide conservation and injection errors.

This power-monitor approach integrates all propagating diffraction orders crossing the monitor, not just the zero order. If interpreting any particular diffracted order, use Lumerical grating projection separately.

### Field normalization

The plane-wave source amplitude is 1 V/m. Under the default CW-normalized frequency-domain monitor convention, `E.E2` corresponds numerically to the dimensionless enhancement `|E|^2/|E0|^2` when `|E0| = 1 V/m`. Do not switch to nonorm before interpreting the field maps. With BFAST and oblique incidence, check the incident-field normalization in the homogeneous reference case before making quantitative enhancement claims.

The export script saves both the requested and actually sampled wavelength. Since the spectral monitor contains discrete samples, a requested resonance wavelength is mapped to the nearest available monitor point. For accurate peak maps, refine sampling or rerun a narrowed monitor/source range as necessary.

### Candidate peaks are not confirmed resonances

The script flags local absorption maxima and transmission dips. It does not guarantee to detect every physical mode or distinguish plasmonic modes from interband absorption, Rayleigh anomalies, Fabry–Pérot artifacts, and numerical fluctuations. Use reflection, a field map, and mesh convergence together to validate a resonance. The assignment calls for every clearly observable resonance and a summary of its movement with angle, which requires manual mode tracking.

## Validation checks before final reporting

1. Remove Au with `run_reference.lsf`. For suspended air at normal incidence, T should approach 1 and R, A approach 0. For bare air/SiO2 at normal incidence, Fresnel theory for n=1/1.45 gives R ≈ ((1.45−1)/(1.45+1))² ≈ 0.03374 and T ≈ 0.96626. Use a nondispersive reference model to test this.
2. Inspect the mesh and material index cross sections: the gold ring must be 40 nm wide with an empty hole, and the substrate must touch its bottom exactly in case B.
3. Compare spectra with 4, 3 and 2 nm metal-region mesh. Record the change in resonance wavelengths, spectral powers and peak field enhancement. Do not assume that satisfactory convergence in T/R implies converged local hotspots.
4. Re-run selected cases with increased vertical PML separation and 24+ PML layers, especially at 75 degrees. Compare against the original run and inspect the simulation status for BFAST divergence.
5. Confirm that the solver reaches auto shutoff rather than only the maximum simulation duration. Increase simulation time for long-lived narrow resonances.
6. Ensure T, R and A are physical for the whole spectrum, and assess tiny negative A in relation to numerical tolerance. Large negatives, >1 values or strong sensitivity to mesh and boundaries invalidate the run.
7. Repeat near-field analysis for each confirmed peak at 0 and 45 degrees in cases A and B. Use consistent map scales between the cases.
8. Check the effect of dispersion choice for Au and SiO2, using the same data set for all comparable simulations.

## Known constraints

- The attached assignment specifies **no numerical convergence tolerance**. Suggested thresholds must be defended empirically, not reported as professor requirements.
- BFAST can require stringent stability settings at very steep angles. The scripts do not automatically tune its time-step multiplier.
- The `.fsp` files at 0/45 degrees can be large, because spatial fields are retained at all 401 wavelengths. If memory is insufficient, set `record_nearfields = 0` during initial spectra runs and perform selected-field wavelength reruns later.
- The `export_nearfields.lsf` script must be fed measured resonance wavelengths from the real simulated spectra, not its numerical placeholders.
- The scripts do **not** automatically assemble the final resonance-angle tracking table or certify that all physical modes are found.

## Technical references

1. [Ansys, Ring geometry and available properties](https://optics.ansys.com/hc/en-us/articles/360034382194-Ring-Simulation-Object)
2. [Ansys, BFAST for broadband oblique incidence](https://optics.ansys.com/hc/en-us/articles/360034902273-Broadband-Fixed-Angle-Source-Technique-BFAST)
3. [Ansys, Reflection monitor positioning](https://optics.ansys.com/hc/en-us/articles/360034915753-Tips-for-accurately-measuring-reflection-in-an-FDTD-simulation)
4. [Ansys, FDTD plasmonics methodology](https://optics.ansys.com/hc/en-us/articles/360041688154-Plasmonics-simulation-methodology)
5. [Ansys, Source wavelength and polarization](https://optics.ansys.com/hc/en-us/articles/360034382854-Plane-wave-and-beam-source-Simulation-object)
6. [Ansys, Monitor and source normalization](https://optics.ansys.com/hc/en-us/articles/360034397034-Units-and-normalization-conventions-in-Lumerical-solvers)
7. [Ansys, Running nested LSF scripts](https://optics.ansys.com/hc/en-us/articles/360034404054-Running-script-files-Script-command)
8. [Ansys, LSF `transmission()` sign conventions](https://optics.ansys.com/hc/en-us/articles/360034405354-transmission-Script-command)
