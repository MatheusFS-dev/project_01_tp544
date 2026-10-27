"""Solve direct saved-field points for every required resonance map."""

import csv
from pathlib import Path

try:
    from .builder import project_path, project_stem
    from .campaign import campaign_manifest
    from .runner import AEDT_VERSION, configure_aedt_environment
except ImportError:
    from builder import project_path, project_stem
    from campaign import campaign_manifest
    from runner import AEDT_VERSION, configure_aedt_environment


ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = ROOT / "results" / "hfss"
RESONANCES_PATH = RESULTS_ROOT / "resonances" / "resonance_candidates.csv"
FIELD_SWEEP = "FieldSweep"


def field_frequencies_by_case() -> dict[tuple[str, int], list[float]]:
    """Return direct field frequencies for each required case and angle.

    Returns:
        Mapping of physical case and angle to unique resonance frequencies in Hz.
    """
    frequencies = {}
    with RESONANCES_PATH.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            angle = int(row["angle_deg"])
            if angle not in (0, 45):
                continue
            key = (row["case"], angle)
            frequency = 299_792_458 / (float(row["wavelength_nm"]) * 1e-9)
            frequencies.setdefault(key, []).append(frequency)
    return {key: sorted(set(values)) for key, values in frequencies.items()}


def field_sweep_inputs(frequencies_hz: list[float]) -> tuple[list[float], list[bool]]:
    """Create independent mutable inputs required by PyAEDT's sweep builder.

    Args:
        frequencies_hz: Direct resonance frequencies in hertz.

    Returns:
        Fresh frequency and field-save flag lists for one PyAEDT call.
    """
    return list(frequencies_hz), [True] * len(frequencies_hz)


def solve_direct_fields() -> None:
    """Add and solve saved direct field points for physical and reference models."""
    configure_aedt_environment()
    from ansys.aedt.core import Hfss

    manifest = campaign_manifest()
    physical = {(item.case_name, item.angle_deg): item for item in manifest.spectra}
    reference = {(item.case_name, item.angle_deg): item for item in manifest.references}
    for key, frequencies_hz in field_frequencies_by_case().items():
        for request in (physical[key], reference[key]):
            project = project_path(RESULTS_ROOT, request)
            hfss = Hfss(
                project=str(project),
                design=project_stem(request),
                solution_type="Modal",
                version=AEDT_VERSION,
                non_graphical=True,
                new_desktop=True,
                close_on_exit=True,
            )
            try:
                if FIELD_SWEEP not in hfss.get_sweeps("Setup1"):
                    sweep_frequencies, save_single_field = field_sweep_inputs(frequencies_hz)
                    sweep = hfss.create_single_point_sweep(
                        setup="Setup1",
                        unit="Hz",
                        freq=sweep_frequencies,
                        name=FIELD_SWEEP,
                        save_single_field=save_single_field,
                        save_fields=True,
                    )
                    if not sweep:
                        raise RuntimeError(f"could not create {FIELD_SWEEP} for {project}")
                if not hfss.analyze_setup("Setup1", cores=4):
                    raise RuntimeError(f"direct field solve failed for {project}")
                if not hfss.save_project(str(project)):
                    raise RuntimeError(f"could not save direct field result: {project}")
            finally:
                hfss.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    solve_direct_fields()
