"""Export power-complete HFSS Floquet spectra in Lumerical-style tables."""

import csv
import json
from pathlib import Path

import numpy as np

try:
    from .analysis import Traspectrum, propagating_modal_tra_spectrum, resample_spectrum
    from .builder import project_path, project_stem
    from .campaign import campaign_manifest
    from .runner import AEDT_VERSION, configure_aedt_environment
except ImportError:
    from analysis import Traspectrum, propagating_modal_tra_spectrum, resample_spectrum
    from builder import project_path, project_stem
    from campaign import campaign_manifest
    from runner import AEDT_VERSION, configure_aedt_environment


ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = ROOT / "results" / "hfss"
PROJECTS_ROOT = RESULTS_ROOT / "projects"
SPECTRA_ROOT = RESULTS_ROOT / "spectra"
EXPORTS_ROOT = RESULTS_ROOT / "exports"
SETUP_SWEEP = "Setup1 : InterpolatingSweep"
SOURCE_MODE = "TopFloquet:2"


def modal_s_parameter_expressions(port_name: str, floquet_modes: int) -> list[str]:
    """Build S-parameter expressions for every output mode of one port.

    Args:
        port_name: HFSS Floquet-port name.
        floquet_modes: Number of retained output modes.

    Returns:
        Ordered modal S-parameter expressions from the TM incident source.
    """
    return [f"S({port_name}:{mode},{SOURCE_MODE})" for mode in range(1, floquet_modes + 1)]


def modal_gamma_expressions(port_name: str, floquet_modes: int) -> list[str]:
    """Build propagation-constant expressions for every mode of one port.

    Args:
        port_name: HFSS Floquet-port name.
        floquet_modes: Number of retained output modes.

    Returns:
        Ordered modal propagation-constant expressions.
    """
    return [f"Gamma({port_name}:{mode})" for mode in range(1, floquet_modes + 1)]


def read_complex_group(hfss, expressions: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Read aligned complex traces from one HFSS solution report.

    Args:
        hfss: Open PyAEDT HFSS application object.
        expressions: HFSS expressions sharing the same frequency sweep.

    Returns:
        Frequencies in hertz and complex values shaped ``traces, samples``.

    Raises:
        RuntimeError: If HFSS returns missing or misaligned traces.
    """
    data = hfss.post.get_solution_data(expressions=expressions, setup_sweep_name=SETUP_SWEEP)
    if not data:
        raise RuntimeError(f"HFSS returned no solution data for {expressions}")
    reference_frequency = None
    values = []
    for expression in expressions:
        frequency_hz, real = data.get_expression_data(expression, formula="real", convert_to_SI=True)
        _, imaginary = data.get_expression_data(expression, formula="imag", convert_to_SI=True)
        frequency_hz = np.asarray(frequency_hz, dtype=float)
        trace = np.asarray(real, dtype=float) + 1j * np.asarray(imaginary, dtype=float)
        if len(frequency_hz) < 2 or frequency_hz.shape != trace.shape:
            raise RuntimeError(f"incomplete HFSS data for {expression}")
        if reference_frequency is None:
            reference_frequency = frequency_hz
        elif not np.allclose(reference_frequency, frequency_hz, rtol=0, atol=0):
            raise RuntimeError(f"inconsistent HFSS frequency axis for {expression}")
        values.append(trace)
    return reference_frequency, np.asarray(values)


def resample_to_report_grid(spectrum: Traspectrum) -> tuple[Traspectrum, np.ndarray]:
    """Convert native HFSS samples to the assignment's 1 nm wavelength grid.

    Args:
        spectrum: Native-frequency T/R/A result.

    Returns:
        One-nanometer spectrum and flags for samples created by interpolation.
    """
    order = np.argsort(spectrum.wavelength_nm)
    sampled = resample_spectrum(
        spectrum.wavelength_nm[order],
        spectrum.transmission[order],
        spectrum.reflection[order],
        np.arange(450.0, 1601.0, 1.0),
    )
    frequency_hz = 299_792_458 / sampled.wavelength_nm * 1e9
    return (
        Traspectrum(
            frequency_hz=frequency_hz,
            wavelength_nm=sampled.wavelength_nm,
            transmission=sampled.transmission,
            reflection=sampled.reflection,
            absorption=1 - sampled.transmission - sampled.reflection,
        ),
        sampled.interpolated,
    )


def write_lumerical_style_table(path: Path, spectrum: Traspectrum) -> None:
    """Write a tabular 1 nm T/R/A spectrum compatible with the review flow.

    Args:
        path: Destination TSV file.
        spectrum: Wavelength-ordered report spectrum.
    """
    order = np.argsort(spectrum.wavelength_nm)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("wavelength_nm", "T", "R", "A"))
        for index in order:
            writer.writerow(
                (
                    f"{spectrum.wavelength_nm[index]:.12g}",
                    f"{spectrum.transmission[index]:.12g}",
                    f"{spectrum.reflection[index]:.12g}",
                    f"{spectrum.absorption[index]:.12g}",
                )
            )


def write_native_table(path: Path, spectrum: Traspectrum, interpolated: np.ndarray) -> None:
    """Write the original HFSS grid with explicit interpolation provenance.

    Args:
        path: Destination CSV file.
        spectrum: Native or report-grid spectrum.
        interpolated: Boolean flags matching the samples in ``spectrum``.

    Raises:
        ValueError: If the provenance flags have a mismatched shape.
    """
    interpolated = np.asarray(interpolated, dtype=bool)
    if interpolated.shape != spectrum.wavelength_nm.shape:
        raise ValueError("interpolation flags must match the spectrum axis")
    order = np.argsort(spectrum.wavelength_nm)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("wavelength_nm", "frequency_hz", "T", "R", "A", "interpolated"))
        for index in order:
            writer.writerow(
                (
                    f"{spectrum.wavelength_nm[index]:.12g}",
                    f"{spectrum.frequency_hz[index]:.12g}",
                    f"{spectrum.transmission[index]:.12g}",
                    f"{spectrum.reflection[index]:.12g}",
                    f"{spectrum.absorption[index]:.12g}",
                    str(bool(interpolated[index])).lower(),
                )
            )


def export_spectrum(request) -> tuple[Path, Path]:
    """Export one completed physical request with all propagating mode powers.

    Args:
        request: Physical spectrum request from the campaign manifest.

    Returns:
        Report-grid TSV path and native-grid CSV path.

    Raises:
        FileNotFoundError: If the corresponding AEDT project is absent.
        RuntimeError: If HFSS cannot supply its solved modal data.
        ValueError: If the retained Floquet modes are insufficient.
    """
    project = project_path(RESULTS_ROOT, request)
    if not project.is_file():
        raise FileNotFoundError(project)
    configure_aedt_environment()
    from ansys.aedt.core import Hfss

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
        reflected_frequency, reflected_s = read_complex_group(
            hfss, modal_s_parameter_expressions("TopFloquet", request.floquet_modes)
        )
        transmitted_frequency, transmitted_s = read_complex_group(
            hfss, modal_s_parameter_expressions("BottomFloquet", request.floquet_modes)
        )
        reflected_gamma_frequency, reflected_gamma = read_complex_group(
            hfss, modal_gamma_expressions("TopFloquet", request.floquet_modes)
        )
        transmitted_gamma_frequency, transmitted_gamma = read_complex_group(
            hfss, modal_gamma_expressions("BottomFloquet", request.floquet_modes)
        )
    finally:
        hfss.release_desktop(close_projects=True, close_desktop=True)
    frequency_axes = (transmitted_frequency, reflected_gamma_frequency, transmitted_gamma_frequency)
    if any(not np.allclose(reflected_frequency, axis, rtol=0, atol=0) for axis in frequency_axes):
        raise RuntimeError("HFSS modal result axes are inconsistent")
    native = propagating_modal_tra_spectrum(
        reflected_frequency,
        reflected_s,
        transmitted_s,
        reflected_gamma,
        transmitted_gamma,
    )
    report, interpolated = resample_to_report_grid(native)
    stem = project_stem(request)
    report_path = SPECTRA_ROOT / f"{request.case_name}_{request.angle_deg}deg.tsv"
    native_path = EXPORTS_ROOT / f"{stem}_native.csv"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    write_lumerical_style_table(report_path, report)
    report_path.with_suffix(".done").touch()
    write_native_table(native_path, native, np.zeros(len(native.frequency_hz), dtype=bool))
    write_native_table(EXPORTS_ROOT / f"{stem}_450_1600nm.csv", report, interpolated)
    metadata_path = EXPORTS_ROOT / f"{stem}.metadata.json"
    metadata_path.write_text(
        json.dumps(
            {
                "source_mode": SOURCE_MODE,
                "floquet_modes_per_port": request.floquet_modes,
                "power_definition": "sum(|S|^2) over modes with real(Gamma)=0",
                "native_samples": len(native.frequency_hz),
                "report_samples": len(report.frequency_hz),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return report_path, native_path


def main() -> None:
    """Export all physical HFSS spectra after the campaign is complete."""
    for request in campaign_manifest().spectra:
        report, native = export_spectrum(request)
        print(json.dumps({"report": str(report), "native": str(native)}))


if __name__ == "__main__":
    main()
