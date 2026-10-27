"""Export and normalize HFSS near fields for detected resonances at 0 and 45 degrees."""

import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from .builder import project_path, project_stem
    from .campaign import campaign_manifest
    from .publication_plot_style import create_figure, save_figure, style_axes
    from .runner import AEDT_VERSION, configure_aedt_environment
except ImportError:
    from builder import project_path, project_stem
    from campaign import campaign_manifest
    from publication_plot_style import create_figure, save_figure, style_axes
    from runner import AEDT_VERSION, configure_aedt_environment


ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = ROOT / "results" / "hfss"
RESONANCES_PATH = RESULTS_ROOT / "resonances" / "resonance_candidates.csv"
FIELDS_ROOT = RESULTS_ROOT / "nearfields"
SETUP_SWEEP = "Setup1 : FieldSweep"


def nearest_saved_field_sample(results_path: Path, requested_wavelength_nm: float) -> tuple[float, float]:
    """Select the saved HFSS field basis sample nearest to a candidate.

    Args:
        results_path: HFSS ``*.results`` directory holding field metadata.
        requested_wavelength_nm: Resonance wavelength from the report grid.

    Returns:
        Saved wavelength in nanometers and its HFSS frequency in hertz.

    Raises:
        ValueError: If the result metadata lacks finite field frequencies.
    """
    frequencies = []
    for metadata in results_path.glob("*_SU.txt"):
        match = re.search(r"^\s*Frequency\s+([0-9.eE+-]+)", metadata.read_text(encoding="utf-8", errors="replace"), re.MULTILINE)
        if match:
            frequencies.append(float(match.group(1)))
    frequencies = np.asarray(sorted(set(frequencies)), dtype=float)
    wavelengths = 299_792_458 / frequencies * 1e9
    if not len(wavelengths) or not np.all(np.isfinite(wavelengths)) or not np.all(np.isfinite(frequencies)):
        raise ValueError(f"invalid HFSS field metadata: {results_path}")
    index = int(np.argmin(np.abs(wavelengths - requested_wavelength_nm)))
    return float(wavelengths[index]), float(frequencies[index])


def read_cartesian_plane(path: Path, horizontal_axis: int, vertical_axis: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a scalar Cartesian HFSS field grid and express coordinates in nm.

    Args:
        path: Text field file written by HFSS ``ExportOnGrid``.
        horizontal_axis: Coordinate-column index for the horizontal plot axis.
        vertical_axis: Coordinate-column index for the vertical plot axis.

    Returns:
        Horizontal coordinates, vertical coordinates, and scalar values shaped
        ``horizontal, vertical``.

    Raises:
        ValueError: If the field file has no complete Cartesian scalar grid.
    """
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            values = [float(value) for value in line.split()]
        except ValueError:
            continue
        if len(values) >= 4:
            rows.append(values)
    data = np.asarray(rows, dtype=float)
    if data.ndim != 2 or data.shape[1] < 4 or not np.all(np.isfinite(data[:, :3])):
        raise ValueError(f"invalid HFSS field grid: {path}")
    horizontal = np.unique(data[:, horizontal_axis])
    vertical = np.unique(data[:, vertical_axis])
    output = np.full((len(horizontal), len(vertical)), np.inf)
    horizontal_lookup = {value: index for index, value in enumerate(horizontal)}
    vertical_lookup = {value: index for index, value in enumerate(vertical)}
    for row in data:
        output[horizontal_lookup[row[horizontal_axis]], vertical_lookup[row[vertical_axis]]] = row[-1]
    if np.any(np.isinf(output)):
        raise ValueError(f"incomplete HFSS field grid: {path}")
    return horizontal * 1e9, vertical * 1e9, output


def export_plane(hfss, output_path: Path, start_nm: list[float], stop_nm: list[float], step_nm: list[float], frequency_hz: float) -> None:
    """Export one scalar electric-field magnitude plane from HFSS.

    Args:
        hfss: Open PyAEDT HFSS project.
        output_path: Field-file destination.
        start_nm: Cartesian lower grid point in model nanometers.
        stop_nm: Cartesian upper grid point in model nanometers.
        step_nm: Cartesian grid increment in model nanometers.
        frequency_hz: Field frequency in hertz.

    Raises:
        RuntimeError: If HFSS does not create the requested field file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result = hfss.post.export_field_file_on_grid(
        quantity="E",
        solution=SETUP_SWEEP,
        file_name=str(output_path),
        grid_type="Cartesian",
        grid_start=start_nm,
        grid_stop=stop_nm,
        grid_step=step_nm,
        is_vector=True,
        intrinsics={"Freq": f"{frequency_hz:.15g}Hz", "Phase": "0deg"},
        export_with_sample_points=True,
        export_in_si_system=True,
    )
    if not result or not output_path.is_file():
        raise RuntimeError(f"HFSS could not export field grid: {output_path}")


def open_hfss(request):
    """Open one solved campaign project for read-only field export.

    Args:
        request: Physical or no-Au reference request.

    Returns:
        Open PyAEDT HFSS application object.
    """
    project = project_path(RESULTS_ROOT, request)
    if not project.is_file():
        raise FileNotFoundError(project)
    configure_aedt_environment()
    from ansys.aedt.core import Hfss

    return Hfss(
        project=str(project),
        design=project_stem(request),
        solution_type="Modal",
        version=AEDT_VERSION,
        non_graphical=True,
        new_desktop=True,
        close_on_exit=True,
    )


def _export_request_planes(request, stem: Path, plane_specs: dict, suffix: str, frequency_hz: float) -> None:
    """Export all requested planes from one project before opening another."""
    hfss = open_hfss(request)
    try:
        for plane, (start, stop, step, _, _) in plane_specs.items():
            export_plane(
                hfss,
                stem.with_name(stem.name + f"_{plane}_{suffix}.fld"),
                start,
                stop,
                step,
                frequency_hz,
            )
    finally:
        hfss.release_desktop(close_projects=True, close_desktop=True)


def normalized_planes(physical_request, reference_request, frequency_hz: float, stem: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Export physical and homogeneous-reference planes and normalize by E0.

    Args:
        physical_request: Au-containing spectrum request.
        reference_request: Matching no-Au reference request.
        frequency_hz: Exact saved HFSS field frequency in hertz.
        stem: Common field-output stem.

    Returns:
        XY coordinates and normalized magnitude, followed by XZ coordinates
        and normalized magnitude.

    Raises:
        ValueError: If the homogeneous-reference amplitude is not positive.
    """
    plane_specs = {
        "xy": ([-160, -160, 35], [160, 160, 35], [2, 2, 0], 0, 1),
        "xz": ([-160, 0, -200], [160, 0, 200], [2, 0, 2], 0, 2),
    }
    _export_request_planes(physical_request, stem, plane_specs, "physical", frequency_hz)
    _export_request_planes(reference_request, stem, plane_specs, "reference", frequency_hz)
    output = []
    for plane, (_, _, _, horizontal_axis, vertical_axis) in plane_specs.items():
        x, y, physical_field = read_cartesian_plane(stem.with_name(stem.name + f"_{plane}_physical.fld"), horizontal_axis, vertical_axis)
        reference_x, reference_y, reference_field = read_cartesian_plane(stem.with_name(stem.name + f"_{plane}_reference.fld"), horizontal_axis, vertical_axis)
        if not np.array_equal(x, reference_x) or not np.array_equal(y, reference_y):
            raise ValueError("physical and reference field grids differ")
        incident_field = float(np.median(reference_field))
        if not np.isfinite(incident_field) or incident_field <= 0:
            raise ValueError("homogeneous-reference incident field is invalid")
        output.extend((x, y, physical_field / incident_field))
    return tuple(output)


def normalized_intensity(field: np.ndarray) -> np.ndarray:
    """Convert normalized electric-field magnitude to normalized intensity.

    Args:
        field: Normalized electric-field magnitude samples.

    Returns:
        ``|E|² / |E0|²`` with any metal masks preserved as ``NaN``.
    """
    return np.square(np.asarray(field, dtype=float))


def field_figure_paths(stem: Path) -> tuple[Path, Path]:
    """Create collision-proof PNG and PDF paths for a field-map stem.

    Args:
        stem: Logical map stem, which may contain a decimal wavelength.

    Returns:
        PNG and PDF output paths preserving the complete stem name.
    """
    return stem.parent / f"{stem.name}.png", stem.parent / f"{stem.name}.pdf"


def plot_plane(horizontal: np.ndarray, vertical: np.ndarray, field: np.ndarray, stem: Path, plane: str, title: str, limit: float) -> None:
    """Save one normalized-field plane with shared color scaling.

    Args:
        horizontal: Horizontal coordinates in nm.
        vertical: Vertical coordinates in nm.
        field: Normalized electric-field magnitude array.
        stem: Output figure stem.
        plane: Plane identifier, either ``xy`` or ``xz``.
        title: Scientific title for the plot.
        limit: Common normalized-intensity upper color limit for comparable maps.
    """
    figure, axes = create_figure(width=5.8, height=4.7)
    intensity = normalized_intensity(field)
    image = axes.pcolormesh(horizontal, vertical, intensity.T, shading="auto", cmap="cividis", vmin=0, vmax=limit)
    figure.colorbar(image, ax=axes, label="|E|² / |E₀|²")
    vertical_label = "y" if plane == "xy" else "z"
    style_axes(axes, "x (nm)", f"{vertical_label} (nm)", title=title, grid_axis="both")
    png_path, pdf_path = field_figure_paths(stem)
    save_figure(figure, png_path)
    save_figure(figure, pdf_path)
    plt.close(figure)


def export_all_fields() -> None:
    """Export and plot fields for every detected 0 and 45 degree candidate.

    Raises:
        FileNotFoundError: If resonance candidates have not been generated.
    """
    if not RESONANCES_PATH.is_file():
        raise FileNotFoundError(RESONANCES_PATH)
    manifest = campaign_manifest()
    physical = {(request.case_name, request.angle_deg): request for request in manifest.spectra}
    references = {(request.case_name, request.angle_deg): request for request in manifest.references}
    FIELDS_ROOT.mkdir(parents=True, exist_ok=True)
    records = []
    with RESONANCES_PATH.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if int(row["angle_deg"]) in (0, 45):
                records.append(row)
    maps = []
    for row in records:
        key = (row["case"], int(row["angle_deg"]))
        stem = FIELDS_ROOT / f"{key[0]}_{key[1]}deg_{row['candidate_id']}_{row['wavelength_nm']}nm"
        field_wavelength_nm = float(row["wavelength_nm"])
        field_frequency_hz = 299_792_458 / (field_wavelength_nm * 1e-9)
        xy_x, xy_y, xy_field, xz_x, xz_z, xz_field = normalized_planes(
            physical[key], references[key], field_frequency_hz, stem
        )
        np.savez_compressed(
            stem.with_suffix(".npz"),
            xy_x_nm=xy_x,
            xy_y_nm=xy_y,
            e_xy=xy_field,
            intensity_xy=normalized_intensity(xy_field),
            xz_x_nm=xz_x,
            xz_z_nm=xz_z,
            e_xz=xz_field,
            intensity_xz=normalized_intensity(xz_field),
            candidate_wavelength_nm=float(row["wavelength_nm"]),
            field_wavelength_nm=field_wavelength_nm,
            field_frequency_hz=field_frequency_hz,
        )
        maps.append((row, stem, xy_x, xy_y, xy_field, xz_x, xz_z, xz_field, field_wavelength_nm))
    xy_limit = max(float(np.nanmax(normalized_intensity(item[4]))) for item in maps)
    xz_limit = max(float(np.nanmax(normalized_intensity(item[7]))) for item in maps)
    for row, stem, xy_x, xy_y, xy_field, xz_x, xz_z, xz_field, field_wavelength_nm in maps:
        title = f"{row['case']}, {row['angle_deg']} degrees, {row['candidate_id']}, field at {field_wavelength_nm:.1f} nm"
        plot_plane(xy_x, xy_y, xy_field, stem.with_name(stem.name + "_xy"), "xy", title, xy_limit)
        plot_plane(xz_x, xz_z, xz_field, stem.with_name(stem.name + "_xz"), "xz", title, xz_limit)


if __name__ == "__main__":
    export_all_fields()
