# HFSS Investigation Plan: Nonphysical Spectral Spikes and Disagreement with Lumerical

Repository: [MatheusFS-dev/project_01_tp544](https://github.com/MatheusFS-dev/project_01_tp544)  
Scope: HFSS Driven Modal implementation, using Lumerical FDTD as an independent comparison.  
Status: Updated investigation plan incorporating the supplied `Treinar.zip` reference model and the exact-match geometry `G0246`. The geometric overlap is confirmed in source code. Its causal effect on the solved fields, and the mechanism of the isolated spectral spikes, still require controlled HFSS reruns.

## 1. Observed issue

The HFSS transmission, reflection, and absorption spectra contain isolated, tooth-like changes that are difficult to reconcile with a continuous physical response. Some features appear at one frequency sample and largely disappear at the adjacent samples. They are present in the HFSS 401-point native-frequency export, before the separate 1 nm report-grid resampling. In addition, HFSS and Lumerical predict substantially different spectra, including for normal incidence (0°). A user-provided HFSS geometry library (`Treinar.zip`) contains a geometry matching the intended annulus and provides an independent structural comparison.

One concrete example is the HFSS `B_SiO2`, 45° spectrum:

| Wavelength (nm) | Transmission (T) | Reflection (R) | Absorption (A) |
| ---: | ---: | ---: | ---: |
| 559.60 | 0.51178 | 0.11299 | 0.37522 |
| 560.86 | 0.33276 | 0.11590 | 0.55133 |
| 562.12 | 0.50561 | 0.11877 | 0.37562 |

Transmission drops by approximately 0.176 at a single sample while reflection remains continuous. Because the analysis defines `A = 1 - T - R`, the transmission deficit is automatically reported as increased absorption. This pattern suggests either incorrect modal power accounting or a numerical sweep artifact, but does not establish the cause.

Relevant files:

- [`src/hfss_scripts/builder.py`](https://github.com/MatheusFS-dev/project_01_tp544/blob/main/src/hfss_scripts/builder.py)
- [`src/hfss_scripts/analysis.py`](https://github.com/MatheusFS-dev/project_01_tp544/blob/main/src/hfss_scripts/analysis.py)
- [`src/hfss_scripts/export_modal_spectra.py`](https://github.com/MatheusFS-dev/project_01_tp544/blob/main/src/hfss_scripts/export_modal_spectra.py)
- [HFSS `B_SiO2`, 45° exported spectrum](https://github.com/MatheusFS-dev/project_01_tp544/blob/main/results/hfss_deliverables/exports/spectrum_B_SiO2_45deg_native.csv)

## 2. Expected corrected behavior

A valid result should satisfy all of the following:

1. **Spectral continuity:** transmission and reflection should vary consistently with wavelength, except where a reproducible physical resonance or diffraction threshold produces a rapid change. Isolated one-point spikes should not survive independent frequency solutions and convergence tests.
2. **Physically credible resonance shape:** the principal nanoring resonance should have a finite, converged linewidth and spatial field enhancement compatible with the geometry and gold losses. It may be relatively broad, as suggested by the Lumerical 0° comparison, but an authentic narrow resonance must not be removed merely because it is narrow.
3. **Correct energy accounting:** propagating reflected and transmitted Floquet orders are included exactly once. Absorbed power independently calculated inside the gold agrees with `1 - R - T` within numerical tolerance. The identity `R + T + A = 1` alone is **not** a validation because `A` is defined from `R` and `T`.
4. **Boundary-condition consistency:** periodic phase shifts, source direction, polarization, and Floquet port scan settings describe the same incident plane wave in both media.
5. **Numerical convergence:** spectra, resonance positions, linewidths, and near-fields remain acceptably stable when the mesh, frequency sampling, number of Floquet modes, and port distance are varied separately.
6. **Agreement between solvers and geometry controls:** after matching geometry, optical constants, illumination, and observable definitions, HFSS and Lumerical should agree in the main spectral structure within documented numerical error. The colleague's `G0246` response offers a geometry-level control but is **not** an exact electromagnetic benchmark, because it uses a different gold model and excitation convention. Neither solver is treated as ground truth until independently validated.
7. **Unambiguous material regions:** the volume occupied by the Au annulus should not also remain assigned as air, unless a deliberately configured and verified HFSS material-overlap policy produces the correct effective domain. Confirm Au bulk losses are enabled and the metal is solved internally.

The physical model is a 320 × 320 nm periodic cell containing an Au nanoring with 50 nm inner radius, 90 nm outer radius, and 30 nm height. Cases: air on both sides (`A_air`) or air/SiO₂ (`B_SiO2`, refractive index 1.45). Incidence is nominally TM in the XZ plane at 0°, 15°, 30°, 45°, 60°, and 75°, over 450–1600 nm.

## 3. Reference-geometry comparison from `Treinar.zip` (new evidence)

### 3.1 Exact geometry match: `G0246`

The supplied archive contains 2,000 spectral CSVs plus `HFSS_Base_Dados.py`, `geometrias_sorteadas.csv`, and `resultados.csv`. In `geometrias_sorteadas.csv`, `G0246` belongs to family `F03` (circular annulus) with `Rout_c = 90 nm` and `w_c = 40 nm`. Thus its inner radius is `90 - 40 = 50 nm`, matching the user's annulus. The colleague's script fixes the thickness at `H_AU = 0.030 µm = 30 nm` and the period at `P = 0.320 µm = 320 nm`, also matching the intended device. It simulates an air-surrounded cell at 0°.

| Quantity | Colleague, `G0246` | User, `A_air`, 0° |
| --- | --- | --- |
| Outer / inner radius | 90 / 50 nm | 90 / 50 nm |
| Au thickness / period | 30 / 320 nm | 30 / 320 nm |
| Surroundings / incidence | Air / 0° | Air / 0° |
| Gold material model | Johnson–Christy analytic fit (`Au_JC`) | Externally imported Palik tabulation |
| Excitation | Average of HFSS incident Floquet modes 1 and 2 | `TopFloquet:2`, nominal TM00 |
| Adaptive frequency | 666.2 THz, approximately 450 nm | 374.7405725 THz, approximately 800 nm |
| Floquet modes | 4 per port in air | 10 per port in air |
| Domain extent | Single air box, approximately 1 µm each side | Two air boxes, 800 nm each side |
| Sweep | 401-point interpolating, tolerance 0.5% | 401-point interpolating, tolerance 0.5% |

The similar geometries make this comparison diagnostically useful, but the differences in optical constants, excitation processing, and solver settings prevent treating the spectra as strictly equivalent.

### 3.2 Confirmed source-code construction difference

In the colleague's `HFSS_Base_Dados.py`, the annulus is created and then the air volume is cut using the preserved Au body:

```python
ouro = pecas[0]
ouro.material_name = "Au_JC"
ouro.solve_inside = True
ar = m.create_box(..., "Air_Domain", "vacuum")
m.subtract(ar, [ouro], keep_originals=True)
```

In the user's [`builder.py`](https://github.com/MatheusFS-dev/project_01_tp544/blob/main/src/hfss_scripts/builder.py), `UpperMedium` spans `0 ≤ z ≤ 800 nm`, and `AuRing` spans `0 ≤ z ≤ 30 nm`. The ring is hollowed out by subtracting its inner cylinder, but **the remaining Au annulus is never subtracted from `UpperMedium`**. Therefore the code constructs overlapping Au and vacuum solids over the annular volume. It also does not explicitly set `AuRing.solve_inside = True`.

**What is established:** geometric overlap is demonstrable from the construction code, and the colleague explicitly avoids it. **What is not established:** whether the generated AEDT project resolves that overlap by an implicit material-priority rule, and whether the overlap is sufficient to explain every anomalous feature. Inspect the effective material assignment in AEDT and make the model unambiguous before drawing a definitive causal conclusion.

A minimal *conceptual* correction, to be inserted after constructing and hollowing the ring, but **before** assigning Floquet ports, lattice boundaries, or mesh, is:

```python
outer.name = "AuRing"
outer.solve_inside = True
if not hfss.modeler.subtract(top, [outer], keep_originals=True):
    raise RuntimeError("could not subtract AuRing from UpperMedium")
```

`top` denotes `UpperMedium` and `outer` denotes the surviving hollow Au ring. Preserve both bodies and revalidate all face IDs after the boolean operation. This snippet illustrates the necessary construction change, not a patch already applied to the repository.

### 3.3 Measured spectral discrepancy near 800 nm

| Solver / source | Actual sample wavelength (nm) | T | R | A |
| --- | ---: | ---: | ---: | ---: |
| Colleague, HFSS, `G0246.csv` | 801.301 | 0.044210 | 0.676938 | 0.278852 |
| User, HFSS, `A_air_0deg_native.csv` | 798.780 | 0.786581 | 0.063367 | 0.150052 |
| User, Lumerical, `A_air_0deg.tsv` | 800.000 | 0.094476 | 0.906298 | -0.000773 |

These are **nearby, not identical**, wavelengths. The colleague's HFSS and the user's Lumerical both show strong reflection and suppressed transmission around 800 nm, whereas the user's HFSS shows the opposite trend. The colleague's `resultados.csv` reports a transmission minimum near 814 nm for `G0246`. The quantitative absorption differs substantially between the colleague's HFSS and the user's Lumerical, reinforcing the material-model and solver-normalization caveats.

**Working diagnosis:** unresolved geometry/material overlap is the leading *structural* suspect for the large 0° discrepancy. The isolated `B_SiO2`, 45° tooth near 560.86 nm may additionally involve Floquet modal bookkeeping or interpolating-sweep behavior. One mechanism must not be assumed to explain both symptoms without controlled evidence.

**Archive traceability:** `Treinar.zip` → `HFSS_Base_Dados.py` (construction around lines 176–200), `geometrias_sorteadas.csv` (`G0246`), `Espectros_TRA\G0246.csv`, and `resultados.csv` (`G0246`). No repository files or AEDT results have been modified or rerun.

## 4. Possible problems, ranked by diagnostic priority

| Priority | Suspected problem | Evidence or risk | Test that distinguishes it |
| --- | --- | --- | --- |
| P0 | **Overlapping Au and vacuum solids / ambiguous material precedence** | The user's upper air box and Au annulus occupy the same volume. The colleague explicitly subtracts Au from air and sets `solve_inside=True`. | Rebuild `A_air`, 0° with the Au void cut from `UpperMedium`. Hold all other settings fixed. Inspect AEDT intersections/materials and compare the 750–850 nm response before and after. |
| P0 | **Incorrect identification or summation of propagating Floquet modes** | `analysis.py` applies `np.isclose(real(Gamma), 0, atol=1e-6)` to classify modes. A false mode transition could make power disappear from T or R. | Export `Gamma` and complex S for every mode around each spike. Compare the per-mode power and propagation mask at neighboring frequencies. |
| P0 | **Interpolating-sweep artifact** | `builder.py` uses `InterpolatingSweep` (`interpolation_tol=0.5`, 401 output frequencies). An output grid is not equivalent to 401 independent FEM solutions. | Run a local **Discrete** sweep using independently solved frequencies and compare against the interpolated trace. |
| P1 | **Incorrect Floquet port setup or mode ordering** | Port creation specifies lattice vectors and number of modes but not explicit scan-angle arguments. Mode 2 is assumed to be the incident TM₀₀ channel. | Inspect the generated AEDT port tables, scan angles, mode labels, power normalization, reference planes, and de-embedding. |
| P1 | **Inconsistent lateral periodic phase / incident wave vector** | Both X and Y lattice pairs use `UseScanAngle` with `phi=0`, `theta=request.angle_deg`. Correctness must be checked against the actual port and source conventions. | Verify the imposed phase across X/Y against `k_parallel · a`, plus the reflected/transmitted wave-vector directions. |
| P1 | **Mesh or adaptive-frequency error** | The adaptive solution is centered at about 800 nm, while some suspicious features are near 500–600 nm. Gold has strong field gradients and a 30 nm thickness. | Independently refine the Au/interface mesh and adapt at a wavelength near the suspect feature. |
| P1 | **Different gold dielectric functions / different reference excitation** | User HFSS imports external Palik, user Lumerical uses built-in Palik, and colleague HFSS uses an analytic Johnson–Christy fit. The colleague also averages two Floquet incident modes while the user excites one. | First match material ε(λ) and polarization definitions. Treat `G0246` as a shape/construction control, not a pointwise quantitative target. |
| P2 | **Unphysical or incomplete absorption measurement** | `A = 1 - R - T` assumes every propagating output is counted correctly. | Integrate material dissipation in Au and compare with S-parameter absorption. |
| P2 | **Too few Floquet modes or port too close to near fields** | The baseline uses 10 modes in air, 14 in SiO₂, and ports 800 nm from the ring plane. | Increase mode count and port distance independently, keeping all other settings fixed. |
| P2 | **Lumerical reference not fully converged or not configured equivalently** | Current Lumerical scripts use `Bloch/periodic`; some repository documentation describes BFAST. Broadband oblique illumination and FDTD shutoff require verification. | Validate source angular behavior over wavelength, actual boundary type, simulation end condition, PMLs, and mesh before comparing curves. |

**Important physical distinction:** a Rayleigh anomaly or a narrow lattice resonance can be real. Verify diffraction cutoffs and field distributions before rejecting a sharp spectral feature. The spike near 560.86 nm for SiO₂ at 45° is a diagnostic target, not proof of a particular failure mechanism.

## 5. Investigation procedure

### Phase 1. Preserve baselines and audit the geometry before new solves

- [ ] Preserve the original AEDT projects, exported S-parameters, T/R/A CSV files, and solver convergence logs. Do not overwrite the original baseline.
- [ ] Record the matched `G0246` design from `Treinar.zip` and reproduce the tabulated 800 nm comparison in §3.3.
- [ ] Open one generated user AEDT project and inspect the 3D solids, model/non-model flags, intersections, Au `Solve Inside`, and actual material assignment in the overlapping annulus.
- [ ] Check the colleague's geometry construction and its single continuous air domain against the user's two-box air domain. Compare the *effective* material topology, not only nominal radii and thicknesses.
- [ ] Select two diagnostic cases: `A_air`, 0° for the missing broad spectral response, and `B_SiO2`, 45° near **560.86 nm** for the isolated tooth.
- [ ] Plot native HFSS, Lumerical, and `G0246` T/R/A without smoothing. Align sample wavelengths explicitly and annotate the different gold models and source conventions.

**Decision:** if the generated AEDT geometry has unresolved overlap or Au is not solved internally, repair that structural defect first. The fact of overlap is confirmed by the scripts, but its effect on the numerical solution remains to be tested.

### Phase 2. Controlled geometry correction, then modal and sweep diagnostics

- [ ] Create a separate copy of the `A_air`, 0° HFSS model. Change **only** the geometry boolean operation to remove `AuRing` from `UpperMedium` and explicitly verify Au `Solve Inside`. Retain original Palik material, source, Floquet modes, mesh, adaptive setup, and sweep settings.
- [ ] Run the original and geometry-corrected models on the same frequency grid. Compare T/R/A around **750–850 nm**, especially the expected reflected-band structure near 800 nm. Record whether the missing response reappears. Do not force a quantitative match to `G0246` while the gold models differ.
- [ ] If geometry alone resolves the large 0° mismatch, apply the same boolean construction correction to a copy of `B_SiO2`, 45°. Measure whether the 560.86 nm tooth changes or remains.
- [ ] Export full complex modal `S(output, TopFloquet:2)` and `Gamma(output)` from both ports around 559.60, 560.86, and 562.12 nm. Record per-mode `|S|²`, propagation masks, and total outgoing power.
- [ ] Determine whether a change in the `real(Gamma)` inclusion mask explains the transmission deficit. Use continuity of mode identification and physical modal power accounting, not a looser arbitrary threshold, to evaluate this.
- [ ] Repeat **555–570 nm** using a Discrete HFSS sweep, initially at **0.25 nm spacing** (convert wavelengths to frequency as required by AEDT).
- [ ] Save an actual field solution at the suspect point, its two neighbors, and an off-resonance reference.
- [ ] Compare discrete-sweep results with the interpolating sweep on the same wavelength axis. Track S-parameter phase continuity, not just magnitude.
- [ ] If the tooth disappears in the discrete sweep, attribute it to the previous frequency-sweep/result-extraction pipeline and confirm which component is responsible.
- [ ] If it persists, perform finer discrete sampling and assess material loss, modal continuity, and local field enhancement before classifying it as physical.

**Decision:** first determine whether the missing normal-incidence resonance responds to the geometric fix. Separately, if the 45° tooth disappears in discrete solves, investigate sweep interpolation. If its missing power matches an excluded Floquet mode, correct the modal accounting. A genuine resonance must reproduce under direct solves and converge with refinement.

### Phase 3. Validate the electromagnetic boundaries and ports

- [ ] Solve **no-Au** controls first: homogeneous air and a bare air/SiO₂ interface at 0° and 45°.
- [ ] Check the air reference against `R ≈ 0`, `T ≈ 1`; check the interface against analytic Fresnel TM reflectance and transmittance, including angular dependence.
- [ ] Verify that the two sets of X/Y periodic pairs match geometrically and implement the correct Bloch phase for propagation in XZ (`phi = 0`).
- [ ] Inspect the actual port definitions in AEDT: lattice origins/vectors, orientation, scan angles, incident-mode label and polarization, and port reference/de-embedding settings.
- [ ] Recalculate the Floquet mode lists using the official **Modes Calculator** at the highest frequency (450 nm) across the relevant scan angles, on **both** ports.
- [ ] Check diffraction-order cutoff wavelengths, especially near the short-wavelength edge in SiO₂ and at oblique incidence. Distinguish a genuine order opening from a discontinuity produced by bookkeeping.
- [ ] Test additional Floquet modes and port offsets independently of mesh changes.

**Decision:** no-Au reference discrepancies indicate a source, port, boundary, normalization, or geometry problem that must be corrected before investigating nanoring resonances.

### Phase 4. Validate geometry, gold dispersion, and mesh

- [ ] Confirm the corrected effective Au annulus volume, material precedence, thickness, radii, and z placement by inspecting the AEDT geometry and cross-sections.
- [ ] Ensure no unintended air, metal, or substrate overlap changes the physical interface. Check the colleague's subtraction operation as a structural precedent, not a validated replacement for auditing your own AEDT project.
- [ ] Export gold optical constants from HFSS and Lumerical and compare ε′ and ε″ throughout 450–1600 nm. Inspect the HFSS conductivity implementation and passivity convention.
- [ ] Run mesh convergence studies in the suspect region, varying only local mesh resolution, with a separate adaptive setup nearer the suspect wavelength.
- [ ] Compare Au volumetric loss to S-parameter-derived absorption for a few representative frequencies.
- [ ] Recheck the Lumerical reference: mesh refinement, PML sensitivity, auto-shutoff, and true broadband oblique incidence settings.

**Decision:** if the isolated spike is fixed but the main 0° resonance still disagrees, prioritize geometry/material/model equivalence rather than further Floquet-mode adjustments.

### Phase 5. Confirm corrections before regenerating results

- [ ] Recompute the two diagnostic cases with the identified change, varying one factor at a time. Keep separate records for geometry-only, port/mode-only, sweep-only, and material-only changes.
- [ ] Confirm no unexplained one-point teeth remain and genuine narrow features have reproducible near-fields and stable linewidths.
- [ ] Confirm the no-Au benchmarks, independent energy balance, and mesh/sweep/mode convergence.
- [ ] Document any remaining HFSS–Lumerical differences and their tested causes. Do not claim agreement based only on visually similar plots.
- [ ] Only then regenerate the full set of 12 HFSS spectra and the resonance/near-field deliverables.

## 6. Evidence and acceptance criteria

For each tested hypothesis, save: **original observation → single change → resulting spectrum → quantitative difference → decision (confirmed, rejected, or unresolved)**.

Minimum deliverables:

1. `geometry_overlap_audit.md`: screenshots or AEDT object/region inspection, overlap detection, Au `Solve Inside`, Material Override status, and verified boolean result.
2. `G0246_reference_comparison.csv`: the colleague's geometry parameters and spectrum compared to the user's 0° HFSS and Lumerical, with explicit wavelength alignment and different-material caveats.
3. `A_air_0deg_geometry_before_after.csv`: controlled geometry-only experiment covering at least 750–850 nm, including T/R/A and resonance metrics.
4. `modal_diagnostics.csv`: wavelength, port, mode, complex S, complex Gamma, inclusion mask, and `|S|²`.
5. `discrete_vs_interpolating.csv` and an unsmoothed comparison plot.
6. `no_Au_Fresnel_validation.csv`: HFSS versus analytic TM Fresnel R/T.
7. `convergence_summary.csv`: mesh, adaptive frequency, Floquet mode count, port distance, principal resonance wavelength, linewidth, and R/T/A differences.
8. `root_cause_summary.md`: confirmed mechanism, changed files/settings, before-and-after evidence, and remaining limitations.

Provisional targets to agree before running: no-Au absolute R/T error ≤ 0.01 and numerical convergence of R/T within 0.01 absolute away from genuine diffraction thresholds. Resonance center and linewidth must stabilize under refinement. For the geometry-only A_air/0° test, evaluate whether a strong resonance appears near the 750–850 nm region, but do **not** require it to reproduce the colleague's `G0246` absorption or center wavelength exactly until gold dispersion and source definitions match. These are **investigation targets**, not evidence that current results pass.

## 7. References

1. [HFSS implementation and datasets, project_01_tp544](https://github.com/MatheusFS-dev/project_01_tp544)
2. [Ansys HFSS, Floquet Ports](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v261/en/Subsystems/HFSS/Content/HFSS/FloquetPorts.htm)
3. [Ansys HFSS, Floquet Port Modes Calculator](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v261/en/Subsystems/HFSS/Content/HFSS/FloquetPortModesCalculator.htm)
4. [Ansys HFSS, Floquet Port Mode Setup](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/HFSS/FloquetPortModeSetup.htm)
5. [Ansys HFSS, Selecting the Sweep Type](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v261/en/Subsystems/HFSS/Content/HFSS/SelectingtheSweepType.htm)
6. [Ansys HFSS, Interpolating Frequency Sweeps](https://ansyshelp.ansys.com/public/views/secured/electronics/v252/en/subsystems/hfss/Content/HFSS/InterpolatingFrequencySweepsinHFSS.htm)
7. [Ansys Optics, Convergence Testing for FDTD](https://optics.ansys.com/hc/en-us/articles/360034915833-Convergence-testing-process-for-FDTD-simulations)
8. Ansys HFSS, [Design Settings / Material Override](https://ansyshelp.ansys.com/public/views/secured/electronics/v252/en/subsystems/hfss/Content/Modeler/DesignSettingsforHFSS.htm).
9. Ansys HFSS, [Solving Inside or on the Surface](https://ansyshelp.ansys.com/public/views/secured/electronics/v261/en/subsystems/hfss/Content/Materials/SolvingInsideorontheSurface.htm).
10. User-supplied reference archive, `Treinar.zip`, specifically `HFSS_Base_Dados.py`, `geometrias_sorteadas.csv`, `Espectros_TRA\G0246.csv`, and `resultados.csv`.
11. [User's normal-incidence HFSS reference spectrum](https://github.com/MatheusFS-dev/project_01_tp544/blob/main/results/hfss_deliverables/exports/spectrum_A_air_0deg_native.csv).
12. [User's normal-incidence Lumerical reference spectrum](https://github.com/MatheusFS-dev/project_01_tp544/blob/main/results/spectra/A_air_0deg.tsv).
