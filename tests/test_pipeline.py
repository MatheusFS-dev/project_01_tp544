"""Synthetic checks for the offline spectrum and near-field pipeline."""

import csv
import importlib.util
from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "src" / "lumerical_scripts"


def load_script(name):
    """Import a numbered Python entry point for isolated testing."""
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def synthetic_spectrum(path):
    """Write one valid synthetic spectrum in reverse wavelength order."""
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("wavelength_nm", "T", "R", "A"))
        for wavelength in range(1600, 449, -1):
            writer.writerow((wavelength, .7, .2, .1))


def test_spectrum_selection_and_outputs(tmp_path):
    """Completed spectra alone generate marked plots and 0/45 requests."""
    review = load_script("03_review_resonances")
    spectrum_dir = tmp_path / "spectra"
    output_dir = tmp_path / "resonances"
    spectrum_dir.mkdir()
    for angle in (0, 15, 45):
        stem = spectrum_dir / f"A_air_{angle}deg"
        synthetic_spectrum(stem.with_suffix(".tsv"))
        if angle != 15:
            stem.with_suffix(".done").write_text("complete")
    spectra = review.completed_spectra(spectrum_dir)
    assert set(spectra) == {("A_air", 0), ("A_air", 45)}
    assert spectra[("A_air", 0)][0, 0] == 450
    manual = output_dir / "resonances.csv"
    assert review.read_resonances(manual, spectra) == []
    original = "case,angle_deg,resonance_id,wavelength_nm\nA_air,0,M1,800\nA_air,45,M1,850\n"
    manual.write_text(original, encoding="utf-8")
    selections = review.read_resonances(manual, spectra)
    review.style()
    review.write_outputs(spectra, selections, output_dir)
    assert manual.read_text(encoding="utf-8") == original
    assert (output_dir / "A_air_0deg_TRA.png").exists()
    assert (output_dir / "A_air_0deg_TRA.pdf").exists()
    assert (output_dir / "resonance_vs_angle.png").exists()
    lines = (output_dir / "nearfield_requests.tsv").read_text().splitlines()
    assert lines == ["request_id\tcase_index\tangle_deg\twavelength_nm", "1\t1\t0\t800", "2\t1\t45\t850"]


@pytest.mark.parametrize("entry", ["A_air,0,M1,800.5", "A_air,0,M1,1700", "A_air,15,M1,800", "other,0,M1,800"])
def test_invalid_resonance_rejected(tmp_path, entry):
    """Manual choices must match the spectral grid and completed cases."""
    review = load_script("03_review_resonances")
    manual = tmp_path / "resonances.csv"
    manual.write_text("case,angle_deg,resonance_id,wavelength_nm\n" + entry + "\n")
    with pytest.raises(ValueError):
        review.read_resonances(manual, {("A_air", 0): np.empty((1151, 4))})


def test_nearfield_plot_and_static_configuration(tmp_path):
    """Field MAT shapes and key LSF settings remain compatible."""
    nearfield = load_script("04_plot_nearfields")
    mat = tmp_path / "field.mat"
    savemat(mat, {"x_nm": [-1, 1], "y_nm": [-1, 1], "z_nm": [-2, 2],
                  "E2_xy": [[1, 2], [3, 4]], "E2_xz": [[2, 3], [4, 5]],
                  "requested_wavelength_nm": 800, "sampled_wavelength_nm": 800})
    data = nearfield.read_fields(mat)
    nearfield.plot_maps(data, tmp_path / "field", "synthetic", {"xy": 4, "xz": 5})
    assert (tmp_path / "field_xy.png").exists()
    assert (tmp_path / "field_xz.pdf").exists()
    settings = (SCRIPTS / "settings_nanorings.lsf").read_text()
    geometry = (SCRIPTS / "build_nanorings.lsf").read_text()
    power = (SCRIPTS / "configure_power_monitors.lsf").read_text()
    spectra = (SCRIPTS / "02_run_spectra.lsf").read_text()
    assert "frequency_points = 1151" in settings
    assert "fdtd_z_min = -2000e-9" in settings and "fdtd_z_max = 2000e-9" in settings
    assert 'set("z min",-200e-9)' in geometry and 'set("z max",200e-9)' in geometry
    assert 'set("dimension",3)' in geometry
    assert 'set("output power",1)' in power and 'set("output Ex",0)' in power
    assert "NF_XY" not in spectra and "NF_XZ" not in spectra
