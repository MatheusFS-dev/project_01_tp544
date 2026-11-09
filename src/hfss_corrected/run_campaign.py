"""Run the assignment with the colleague's HFSS setup, without solver tuning."""

import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import subprocess
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

import reference as ref


# ---------------------------------------------------------------------------- #
#                                    Config                                    #
# ---------------------------------------------------------------------------- #

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results/hfss_corrected"
REFERENCE = Path.home() / "Downloads/Treinar/HFSS_Base_Dados.py"
REFERENCE_CSV = REFERENCE.parent / "Espectros_TRA/G0246.csv"
ANGLES = (0, 15, 30, 45, 60, 75)
CASES = ("A_air", "B_SiO2")

# ---------------------------------------------------------------------------- #

sys.path.insert(0, str(ROOT / "src/hfss_scripts"))
from publication_plot_style import apply_publication_style


def save_json(path, value):
    """Write campaign metadata atomically.

    Args:
        path: Destination path.
        value: JSON-serializable metadata.
    """
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, default=str) + "\n")
    temporary.replace(path)


def plot_spectrum(data, name, invalid_below_nm=None):
    """Plot the native T/R/A samples without smoothing or clipping.

    Args:
        data: Array with wavelength, transmission, reflection and absorption.
        name: Case and angle identifier used for the title and filenames.
        invalid_below_nm: Optional cutoff below which the saved diagnostic
            omits a propagating port mode and cannot be validated.
    """
    apply_publication_style()
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
    for column, label, color in ((1, "T", "#1f4e79"), (2, "R", "#5b6c5d"), (3, "A", "#b54708")):
        ax.plot(data[:, 0], data[:, column], label=label, color=color, linewidth=1.4)
    ax.set(xlabel="Comprimento de onda (nm)", ylabel="Fração de potência", title=name.replace("_", " ") + " · Au JC · TM", xlim=(450, 1600))
    ax.set_ylim(min(-.02, float(data[:, 1:].min()) - .02), max(1.02, float(data[:, 1:].max()) + .02))
    if invalid_below_nm is not None:
        ax.axvspan(450, invalid_below_nm, color="0.4", alpha=.12)
        ax.text(.02, .04, "Diagnóstico: faixa sombreada com modos incompletos", transform=ax.transAxes, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(ncol=3, frameon=False)
    for extension in ("png", "pdf"):
        fig.savefig(OUTPUT / "plots" / f"{name}.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def check_geometry(hfss, case):
    """Verify live AEDT material volumes against the analytic assignment.

    Args:
        hfss: Active design after geometry construction.
        case: Assignment case identifier.

    Returns:
        Dictionary with actual material volumes and analytical errors.

    Raises:
        RuntimeError: If Au is not solved inside or volume checks fail.
    """
    gold = hfss.modeler["Au_anelc"]
    air = hfss.modeler["Air_Domain"]
    expected = math.pi * (.090 ** 2 - .050 ** 2) * ref.H_AU
    air_z = ref.H_AU + ref.H_AR + (ref.H_AR if case == "A_air" else 0)
    expected_air = ref.P ** 2 * air_z - expected
    raw_solve_inside = hfss.modeler.oeditor.GetPropertyValue("Geometry3DAttributeTab", gold.name, "Solve Inside")
    solve_inside = str(raw_solve_inside).lower() == "true"
    result = {"gold_volume_um3": gold.volume, "expected_gold_volume_um3": expected,
              "air_volume_um3": air.volume, "expected_air_volume_um3": expected_air,
              "gold_solve_inside": solve_inside, "gold_material": gold.material_name,
              "pyaedt_cached_solve_inside": gold.solve_inside}
    if not solve_inside or not math.isclose(gold.volume, expected, rel_tol=1e-4):
        raise RuntimeError(f"Au geometry validation failed: {result}")
    if not math.isclose(air.volume, expected_air, rel_tol=1e-4):
        raise RuntimeError(f"Au was not correctly subtracted from air: {result}")
    if case == "B_SiO2":
        result["substrate_volume_um3"] = hfss.modeler["SiO2_Domain"].volume
    return result


def solve_with_progress(hfss, name):
    """Solve the unchanged reference setup while reporting elapsed time.

    Args:
        hfss: Active HFSS design with Setup1.
        name: Case identifier printed in the progress messages.

    Raises:
        RuntimeError: If HFSS reports a failed solve.
    """
    started = time.monotonic()
    # A separate process can report progress while AEDT's RPC holds Python's GIL.
    heartbeat_code = (
        "import time\n"
        f"started = {started!r}\n"
        "while True:\n"
        "    time.sleep(30)\n"
        f"    print('Solving {name}: elapsed %.1f min' % ((time.monotonic() - started) / 60), flush=True)\n"
    )
    heartbeat = subprocess.Popen([sys.executable, "-u", "-c", heartbeat_code])
    try:
        if not hfss.analyze_setup("Setup1"):
            raise RuntimeError(f"HFSS failed to solve {name}")
    finally:
        heartbeat.terminate()
        heartbeat.wait(timeout=5)


def check_mode_coverage(project, case, angle):
    """Reject saved port tables that omit a propagating diffraction channel.

    Args:
        project: Saved AEDT project with its explicit Floquet mode tables.
        case: Assignment case identifier.
        angle: Incident angle in air, in degrees.

    Raises:
        RuntimeError: If a table is missing or omits a propagating TE/TM order.
    """
    source = project.read_text(errors="replace")
    for port, index in (("FloquetCima", 1.0), ("FloquetBaixo", 1.45 if case == "B_SiO2" else 1.0)):
        section = re.search(r"\$begin '" + port + r"'(.*?)\$end '" + port + r"'", source, re.S)
        if not section:
            raise RuntimeError(f"Missing saved mode table: {port}")
        modes = re.findall(r"ModeNumber=(\d+)\s+IndexM=(-?\d+)\s+IndexN=(-?\d+).*?PolarizationState='(TE|TM)'", section.group(1), re.S)
        configured = {(int(m), int(n), polarization) for _, m, n, polarization in modes}
        required = set()
        # All nonzero orders have increasing transverse norm over this band,
        # since lambda/P >= 450/320 > sin(theta). Its short edge is sufficient.
        for m in range(-3, 4):
            for n in range(-3, 4):
                transverse = (math.sin(math.radians(angle)) + m * 450 / 320) ** 2 + (n * 450 / 320) ** 2
                if transverse < index ** 2:
                    required.update((m, n, polarization) for polarization in ("TE", "TM"))
        missing = sorted(required - configured)
        if missing:
            raise RuntimeError(f"{port} omits propagating modes {missing}. Preserving the reference table prevents validation of this case. No table was changed.")


def compare_baselines(data, name):
    """Compare raw results with the old HFSS and matched colleague geometry.

    Args:
        data: Native wavelength and T/R/A array.
        name: Stable case identifier.

    Returns:
        Measured differences, with explicit wavelength interpolation.

    Raises:
        RuntimeError: If the new result duplicates the old HFSS spectrum.
    """
    result = {}
    old_path = ROOT / "results/hfss_deliverables/exports" / f"spectrum_{name}_native.csv"
    old = np.genfromtxt(old_path, delimiter=",", names=True)
    aligned = np.column_stack([np.interp(data[:, 0], old["wavelength_nm"], old[k]) for k in ("T", "R", "A")])
    result["old_hfss_max_absolute_difference"] = float(np.max(np.abs(data[:, 1:] - aligned)))
    result["old_hfss_rms_difference"] = float(np.sqrt(np.mean((data[:, 1:] - aligned) ** 2)))
    if result["old_hfss_max_absolute_difference"] < 1e-5:
        raise RuntimeError("Corrected spectrum duplicates the old HFSS. Investigation required.")
    if name == "A_air_0deg":
        colleague = np.genfromtxt(REFERENCE_CSV, delimiter=",", names=True)
        target = np.column_stack([np.interp(data[:, 0], colleague["lambda_nm"], colleague[k]) for k in ("T", "R", "A")])
        result["G0246_rms_difference"] = float(np.sqrt(np.mean((data[:, 1:] - target) ** 2)))
        result["G0246_note"] = "Reference averages two incident polarizations. Corrected uses TM only."
        np.savetxt(OUTPUT / "exports/G0246_comparison.csv", np.column_stack((data, target, aligned)), delimiter=",", header="lambda_nm,T_new,R_new,A_new,T_G0246,R_G0246,A_G0246,T_old,R_old,A_old", comments="")
    index = int(np.argmin(abs(data[:, 0] - 800)))
    result["sample_near_800nm"] = data[index].tolist()
    return result


def run_one(hfss, case, angle, project, audit_path):
    """Construct, validate, solve and export one assignment case.

    Args:
        hfss: Active design for this case.
        case: Suspended or substrate case name.
        angle: Incidence angle in degrees.
        project: Persistent AEDT project destination.
        audit_path: Persistent audit destination.
    """
    name = f"{case}_{angle}deg"
    ref.CASO = "A" if case == "A_air" else "B"
    audit = json.loads(audit_path.read_text()) if audit_path.exists() else {}
    if not audit.get("built"):
        print(f"Building {name}: reference G0246, Au JC, original solver setup", flush=True)
        ref.materiais(hfss)
        if case == "B_SiO2":
            silica = hfss.materials.add_material("SiO2_reference")
            silica.permittivity = 1.45 ** 2
        air, substrate, ztop, zbottom = ref.construir_caso(hfss, ref.CASO, angle)
        geometry = check_geometry(hfss, case)
        ports = ref.fronteiras_e_portas(hfss, air, substrate, ztop, zbottom)
        if len(ports) != 2 or not all(ports):
            raise RuntimeError("Floquet port creation failed")
        nbottom = ref.MODOS_VIDRO if case == "B_SiO2" else ref.MODOS_AR
        sources = {f"{port.name}:{mode}": ("1W" if index == 0 and mode == 2 else "0W", "0deg")
                   for index, port in enumerate(ports)
                   for mode in range(1, (ref.MODOS_AR if index == 0 else nbottom) + 1)}
        if not hfss.edit_sources(sources):
            raise RuntimeError("TM source selection failed")
        ref.setup_e_varredura(hfss)
        validation, valid = hfss.validate_full_design(output_dir=str(OUTPUT / "audits"), ports=ref.MODOS_AR + nbottom)
        if not valid:
            raise RuntimeError(f"HFSS design validation failed: {validation}")
        hfss.save_project(str(project))
        audit = {"built": True, "solved": False, "case": case, "angle_deg": angle,
                 "reference_source": str(REFERENCE), "reference_sha256": hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
                 "geometry": geometry, "validation": validation,
                 "boundaries": {b.name: dict(b.props) for b in hfss.boundaries},
                 "setup": {"adaptive_frequency_thz": 666.2, "max_passes": 15, "delta_s": .01,
                           "sweep": "Interpolating", "points": 401, "interpolation_tolerance": .5},
                 "material": "Au_JC, unchanged colleague analytical fit", "source": "FloquetCima:2, TM00",
                 "adaptations": ["G0246 circular annulus", "requested angle", "SiO2 half space for B", "TM instead of polarization average"]}
        save_json(audit_path, audit)
    check_mode_coverage(project, case, angle)
    if not audit.get("solved"):
        print(f"Starting solver: {name}", flush=True)
        started = time.monotonic()
        solve_with_progress(hfss, name)
        hfss.save_project(str(project))
        audit["solve_seconds"] = time.monotonic() - started
        audit["solved"] = True
        save_json(audit_path, audit)
    print(f"Exporting native T/R/A: {name}", flush=True)
    nbottom = ref.MODOS_VIDRO if case == "B_SiO2" else ref.MODOS_AR
    wavelength, transmission, reflection = ref.espectro(hfss, "FloquetCima", "FloquetBaixo", 2, ref.MODOS_AR, nbottom)
    data = np.array(sorted(zip(wavelength, transmission, reflection, [1 - r - t for r, t in zip(reflection, transmission)])))
    if data.shape != (401, 4) or not np.isfinite(data).all():
        raise RuntimeError(f"Invalid spectrum size or nonfinite values: {data.shape}")
    audit["comparison"] = compare_baselines(data, name)
    audit["ranges"] = {k: [float(data[:, j].min()), float(data[:, j].max())] for j, k in enumerate(("lambda_nm", "T", "R", "A"))}
    results = project.with_name(project.name + "results")
    records = list(results.rglob("*.EMConvData"))
    audit["adaptive_converged"] = any("IsConverged=true" in p.read_text(errors="replace") for p in records)
    audit["convergence_records"] = [str(p) for p in records]
    hfss.export_convergence("Setup1", output_file=str(OUTPUT / "audits" / f"{name}_convergence.csv"))
    path = OUTPUT / "exports" / f"{name}.csv"
    temporary = path.with_suffix(".tmp")
    np.savetxt(temporary, data, delimiter=",", header="lambda_nm,T,R,A", comments="", fmt="%.12g")
    temporary.replace(path)
    plot_spectrum(data, name)
    audit["exported"] = True
    save_json(audit_path, audit)
    print(f"Finished {name}: convergence={audit['adaptive_converged']}, {audit['comparison']}", flush=True)


def main():
    """Run all twelve spectra, resuming completed cases on the next launch."""
    os.environ["ANSYSEM_ROOT252"] = "/usr/local/AnsysEM25R2/v252/AnsysEM"
    for directory in ("projects", "exports", "plots", "audits"):
        (OUTPUT / directory).mkdir(parents=True, exist_ok=True)
    from ansys.aedt.core import Hfss

    campaign = [(case, angle) for case in CASES for angle in ANGLES]
    for case, angle in tqdm(campaign, desc="HFSS cases", unit="case"):
        name = f"{case}_{angle}deg"
        project = OUTPUT / "projects" / f"{name}.aedt"
        audit_path = OUTPUT / "audits" / f"{name}.json"
        if audit_path.exists():
            previous = json.loads(audit_path.read_text())
            if previous.get("exported") and not previous.get("diagnostic_only"):
                print(f"Already exported: {name}", flush=True)
                continue
        print(f"Opening HFSS for {name}", flush=True)
        hfss = Hfss(project=str(project), design=name, solution_type="Modal", version="2025.2", non_graphical=True, new_desktop=True)
        try:
            run_one(hfss, case, angle, project, audit_path)
        finally:
            hfss.save_project(str(project))
            hfss.release_desktop(close_projects=True, close_desktop=True)
    print(f"All {len(campaign)} requested T/R/A spectra exported to {OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
