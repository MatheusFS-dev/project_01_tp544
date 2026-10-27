"""Export the completed air, normal-incidence spectrum from HFSS."""

import csv
import json
from pathlib import Path

import numpy as np

try:
    from .analysis import ResampledSpectrum, Traspectrum, resample_spectrum, specular_tra_spectrum
    from .runner import AEDT_VERSION, configure_aedt_environment
except ImportError:
    from analysis import ResampledSpectrum, Traspectrum, resample_spectrum, specular_tra_spectrum
    from runner import AEDT_VERSION, configure_aedt_environment


ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "results" / "hfss" / "projects" / "spectrum_A_air_0deg.aedt"
OUTPUT_DIR = ROOT / "results" / "hfss" / "exports"
SETUP_SWEEP = "Setup1 : InterpolatingSweep"
INCIDENT_MODE = "TopFloquet:2"
REFLECTION_EXPRESSION = "S(TopFloquet:2,TopFloquet:2)"
TRANSMISSION_EXPRESSION = "S(BottomFloquet:2,TopFloquet:2)"


def complex_s_parameter_data(hfss, expression: str) -> tuple[np.ndarray, np.ndarray]:
    """Read one complex native-frequency S-parameter trace from HFSS.

    Args:
        hfss: Open PyAEDT HFSS Driven Modal application object.
        expression: Modal S-parameter expression recognized by HFSS.

    Returns:
        Strictly increasing frequencies in hertz and the associated complex
        S-parameter values.

    Raises:
        RuntimeError: If HFSS does not return a complete solution trace.
    """
    data = hfss.post.get_solution_data(expressions=expression, setup_sweep_name=SETUP_SWEEP)
    if not data:
        raise RuntimeError(f"HFSS returned no solution data for {expression}")
    frequency_hz, real = data.get_expression_data(expression, formula="real", convert_to_SI=True)
    _, imaginary = data.get_expression_data(expression, formula="imag", convert_to_SI=True)
    frequency_hz = np.asarray(frequency_hz, dtype=float)
    values = np.asarray(real, dtype=float) + 1j * np.asarray(imaginary, dtype=float)
    if len(frequency_hz) < 2 or frequency_hz.shape != values.shape:
        raise RuntimeError(f"incomplete HFSS solution data for {expression}")
    return frequency_hz, values


def write_spectrum_csv(path: Path, spectrum: Traspectrum, interpolated: np.ndarray | None = None) -> None:
    """Write one T/R/A spectrum in ascending wavelength order.

    Args:
        path: CSV output destination.
        spectrum: Native or resampled spectrum to export.
        interpolated: Optional Boolean flags identifying non-native samples.

    Raises:
        ValueError: If provenance flags do not match the spectrum length.
    """
    wavelength_nm = np.asarray(spectrum.wavelength_nm, dtype=float)
    if interpolated is None:
        interpolated = np.zeros(len(wavelength_nm), dtype=bool)
    interpolated = np.asarray(interpolated, dtype=bool)
    if interpolated.shape != wavelength_nm.shape:
        raise ValueError("interpolation flags must match the spectrum axis")
    order = np.argsort(wavelength_nm)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("wavelength_nm", "frequency_hz", "transmission", "reflection", "absorption", "interpolated"),
        )
        writer.writeheader()
        for index in order:
            writer.writerow(
                {
                    "wavelength_nm": f"{wavelength_nm[index]:.12g}",
                    "frequency_hz": f"{spectrum.frequency_hz[index]:.12g}",
                    "transmission": f"{spectrum.transmission[index]:.12g}",
                    "reflection": f"{spectrum.reflection[index]:.12g}",
                    "absorption": f"{spectrum.absorption[index]:.12g}",
                    "interpolated": str(bool(interpolated[index])).lower(),
                }
            )


def resampled_spectrum(spectrum: Traspectrum) -> tuple[Traspectrum, np.ndarray]:
    """Create the required 450 to 1600 nm, one-nanometer report grid.

    Args:
        spectrum: Native HFSS spectrum.

    Returns:
        Resampled T/R/A spectrum and its non-native sample flags.
    """
    order = np.argsort(spectrum.wavelength_nm)
    sampled: ResampledSpectrum = resample_spectrum(
        spectrum.wavelength_nm[order],
        spectrum.transmission[order],
        spectrum.reflection[order],
        np.arange(450.0, 1601.0, 1.0),
    )
    frequency_hz = 299_792_458 / sampled.wavelength_nm * 1e9
    transmission = sampled.transmission
    reflection = sampled.reflection
    return (
        Traspectrum(
            frequency_hz=frequency_hz,
            wavelength_nm=sampled.wavelength_nm,
            transmission=transmission,
            reflection=reflection,
            absorption=1 - reflection - transmission,
        ),
        sampled.interpolated,
    )


def export_air_normal_spectrum(project: Path = PROJECT, output_dir: Path = OUTPUT_DIR) -> tuple[Path, Path]:
    """Export the converged A_air, 0 degree spectrum from its AEDT project.

    The 320 nm period is smaller than every wavelength in the 450 to 1600 nm
    range, so both sides in air have only the specular propagating order. The
    remaining solved Floquet modes are evanescent and are not summed as power.

    Args:
        project: Solved A_air normal-incidence AEDT project.
        output_dir: Directory receiving native and one-nanometer CSV files.

    Returns:
        Paths to the native and resampled CSV outputs.

    Raises:
        FileNotFoundError: If the solved AEDT project is absent.
        RuntimeError: If the two S-parameter traces have inconsistent axes.
    """
    if not project.is_file():
        raise FileNotFoundError(project)
    configure_aedt_environment()
    from ansys.aedt.core import Hfss

    hfss = Hfss(
        project=str(project),
        design="spectrum_A_air_0deg",
        solution_type="Modal",
        version=AEDT_VERSION,
        non_graphical=True,
        new_desktop=True,
        close_on_exit=True,
    )
    try:
        reflection_frequency, reflection_s = complex_s_parameter_data(hfss, REFLECTION_EXPRESSION)
        transmission_frequency, transmission_s = complex_s_parameter_data(hfss, TRANSMISSION_EXPRESSION)
    finally:
        hfss.release_desktop(close_projects=True, close_desktop=True)
    if not np.array_equal(reflection_frequency, transmission_frequency):
        raise RuntimeError("HFSS reflection and transmission axes differ")
    spectrum = specular_tra_spectrum(reflection_frequency, reflection_s, transmission_s)
    native_path = output_dir / "spectrum_A_air_0deg_native.csv"
    report_path = output_dir / "spectrum_A_air_0deg_450_1600nm.csv"
    write_spectrum_csv(native_path, spectrum)
    report_spectrum, interpolated = resampled_spectrum(spectrum)
    write_spectrum_csv(report_path, report_spectrum, interpolated)
    metadata = {
        "case": "A_air",
        "angle_deg": 0,
        "incident_mode": INCIDENT_MODE,
        "reflection_expression": REFLECTION_EXPRESSION,
        "transmission_expression": TRANSMISSION_EXPRESSION,
        "power_channels": "specular modes only, as all higher orders are evanescent over 450-1600 nm in air",
        "native_samples": len(spectrum.frequency_hz),
        "report_samples": len(report_spectrum.frequency_hz),
    }
    (output_dir / "spectrum_A_air_0deg.metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return native_path, report_path


if __name__ == "__main__":
    native, report = export_air_normal_spectrum()
    print(json.dumps({"native": str(native), "report": str(report)}, indent=2))
