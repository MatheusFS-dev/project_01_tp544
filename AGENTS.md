# Project Instructions

Global task gates and reusable skills are installed at user scope. Codex and Antigravity read this file directly. Claude reads it through `CLAUDE.md`.

- Treat repository files and executable checks as the source of truth when prose context is stale.
- Inspect existing build, test, lint, and type-check configuration before inventing commands.
- Keep changes scoped to the requested behavior and preserve established architecture and naming.
- Do not modify installed skills unless the user explicitly requests a skill change.
- Do not commit datasets, model weights, credentials, environment files, generated logs, or large artifacts unless explicitly requested.

Add only durable repository-specific information below, such as non-obvious architecture constraints, canonical commands, protected paths, and required validation.


# On the problem

The authoritative assignment is `docs/problem.pdf`; `docs/professor_instructions.md` contains supplementary professor guidance. This repository addresses **Tarefa 3: a 2D plasmonic metasurface of circular Au nanorings**.

## Fixed physical configuration

- Unit cell: circular Au ring with outer radius `Rout = 90 nm`, inner radius `Rin = 50 nm`, ring width `w = 40 nm`, and metal thickness `h = 30 nm`.
- Lattice: square periodic array with `Px = Py = 320 nm` and periodicity along `x` and `y`.
- Compare two separate cases: (A) the structure suspended in air and (B) the structure on SiO2.
- Use dispersive Au optical properties from the simulator/database. The professor specifically permits the Palik Au data supplied with Lumerical and wants the same Au dataset used consistently.
- The upper medium is air (`n = 1`). Use `n ≈ 1.45` or an equivalent dispersive model for SiO2.
- In case B, model SiO2 as semi-infinite or sufficiently thick to prevent the lower interface from producing artificial resonances.
- Do not include an adhesion layer.

## Illumination and required calculations

- Simulate wavelengths from `450 nm` to `1600 nm`.
- Use TM (p) polarization in the `xz` incidence plane at `theta = 0, 15, 30, 45, 60, 75 degrees`, where `theta` is measured from the surface normal (`z`). At normal incidence use `E || x`.
- For every angle and both substrate cases, calculate and plot transmission `T(lambda)`, reflection `R(lambda)`, and absorption `A(lambda) = 1 - R(lambda) - T(lambda)`.
- Identify and mark every clearly observable resonance at each angle and report its resonant wavelength `lambda_r`.
- Produce near-field maps for every identified resonance only at `theta = 0 degrees` and `theta = 45 degrees`. Each resonance needs a representative horizontal plane and a vertical `xz` cut through the unit-cell center, showing either normalized field magnitude `|E|/|E0|` or intensity `|E|^2/|E0|^2`.
- Near-field maps are not required at `theta = 15, 30, 60, 75 degrees`.
- Keep plotting scales and conventions consistent between the suspended and SiO2 cases whenever possible.

## Minimum report outputs

- T/R/A spectra at all six angles for both substrate cases.
- Resonant wavelengths at every angle for both cases.
- Required horizontal and `xz` near-field maps at 0 and 45 degrees.
- A table or plot tracking each resonance versus incidence angle, including at minimum the angle, resonance identifier, `lambda_r` without substrate, and `lambda_r` with substrate.
