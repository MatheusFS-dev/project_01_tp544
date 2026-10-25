"""Synthetic checks for the offline spectrum and near-field pipeline."""

import csv
import importlib.util
from pathlib import Path

import numpy as np
import pytest
import h5py
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


def synthetic_fast_spectrum(path):
    """Write one valid 5 nm exploratory spectrum in reverse wavelength order."""
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("wavelength_nm", "T", "R", "A"))
        for wavelength in range(1600, 449, -5):
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


def test_spectrum_plot_marks_resonance_with_exact_nm_in_legend(monkeypatch):
    """A selected resonance is a plotted point with an exact-nm legend label."""
    review = load_script("03_review_resonances")
    captured = {}

    def capture(fig, _stem):
        captured["figure"] = fig

    monkeypatch.setattr(review, "save_figure", capture)
    data = np.array([[799., .2, .3, .5], [800., .1, .8, .1], [801., .2, .4, .4]])
    review.plot_spectrum(data, [("A_air", 0, "M1", 800)], Path("synthetic"), "synthetic")
    axis = captured["figure"].axes[0]
    assert "M1: 800 nm" in [text.get_text() for text in axis.get_legend().get_texts()]
    assert any(np.array_equal(collection.get_offsets(), [[800., .8]]) for collection in axis.collections)


def test_fast_spectrum_selection_and_outputs(tmp_path):
    """Fast 5 nm spectra use isolated review and near-field artifacts."""
    assert (SCRIPTS / "03_review_fast_resonances.py").exists()
    assert (SCRIPTS / "04_plot_fast_nearfields.py").exists()
    review = load_script("03_review_fast_resonances")
    spectrum_dir = tmp_path / "spectra"
    output_dir = tmp_path / "resonances_fast"
    spectrum_dir.mkdir()
    for angle in (0, 15, 45):
        stem = spectrum_dir / f"fast_A_air_{angle}deg"
        synthetic_fast_spectrum(stem.with_suffix(".tsv"))
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
    review.base.style()
    review.base.write_outputs(spectra, selections, output_dir)
    assert (output_dir / "A_air_0deg_TRA.png").exists()
    assert (output_dir / "nearfield_requests.tsv").read_text().splitlines() == [
        "request_id\tcase_index\tangle_deg\twavelength_nm", "1\t1\t0\t800", "2\t1\t45\t850"]
    fast_run = (SCRIPTS / "04_run_fast_nearfields.lsf").read_text()
    assert 'simulation_profile = "fast";' in fast_run
    assert "../../results/resonances_fast/nearfield_requests.tsv" in fast_run
    assert "../../results/nearfields_fast/request_" in fast_run
    assert '../../results/geometry/"+output_tag+case_name' in fast_run


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
    geometry_driver = SCRIPTS / "build_geometry_driver.lsf"
    power = (SCRIPTS / "configure_power_monitors.lsf").read_text()
    excitation = (SCRIPTS / "configure_excitation.lsf").read_text()
    spectra_entry = (SCRIPTS / "02_run_spectra.lsf").read_text()
    spectra_driver = SCRIPTS / "run_spectra_driver.lsf"
    spectra = spectra_driver.read_text()
    reference = (SCRIPTS / "run_reference.lsf").read_text()
    nearfield_run = (SCRIPTS / "04_run_nearfields.lsf").read_text()
    assert "frequency_points = 1151" in settings
    assert 'if (exist("simulation_profile") != 1)' in settings
    assert 'simulation_profile = "final";' in settings
    assert 'if (simulation_profile == "fast")' in settings
    assert "frequency_points = 231" in settings
    assert "local_mesh_step = 6e-9" in settings
    assert "fdtd_z_min = -1200e-9" in settings and "fdtd_z_max = 1200e-9" in settings
    assert "local_mesh_step = 3e-9" in settings
    assert "mesh_x_span = 220e-9" in settings
    assert "mesh_z_min = -30e-9" in settings and "mesh_z_max = 60e-9" in settings
    assert "simulation_duration = 2500e-15" in settings
    assert "simulation_duration = 10000e-15" in settings
    assert "auto_shutoff = 1e-6" in settings
    assert "pml_layers = 16" in settings
    assert "script_directory = filedirectory(currentscriptname);" in settings
    assert "cd(script_directory);" in settings
    assert "fdtd_z_min = -800e-9" in settings and "fdtd_z_max = 800e-9" in settings
    assert 'set("x span",mesh_x_span)' in geometry
    assert 'set("z min",mesh_z_min)' in geometry and 'set("z max",mesh_z_max)' in geometry
    assert 'set("dimension","3D")' in geometry
    assert 'set("index",1)' in geometry
    assert 'set("background index",1)' not in geometry
    fdtd_block = geometry.split("addfdtd;", 1)[1].split("# At normal incidence", 1)[0]
    assert 'set("name","FDTD")' not in fdtd_block
    assert 'set("output power",1)' in power and 'set("output Ex",0)' in power
    assert 'getnamednumber("R_monitor")' in power
    assert 'getnamednumber("T_monitor")' in power
    assert 'getnamednumber("TM_incident")' in excitation
    assert 'set("plane wave type","BFAST")' not in excitation
    assert 'set("plane wave type","Bloch/periodic")' in excitation
    assert 'set("x min bc","Bloch")' in excitation
    assert 'set("set based on source angle",1)' in excitation
    for driver in (spectra, reference, nearfield_run):
        assert driver.index('save(output_base+".fsp");') < driver.rindex("run")
        assert driver.count('save(output_base+".fsp");') == 1
        assert 'run("FDTD","CPU")' not in driver
    assert 'run("FDTD","GPU")' in spectra
    assert 'run("FDTD","GPU")' in reference
    assert 'run("FDTD","GPU")' in nearfield_run
    assert 'if (simulation_profile == "final")' in spectra
    assert "} else {\n                run;" in spectra
    assert geometry_driver.exists()
    assert spectra_driver.exists()
    assert "build_geometry_driver;" in (SCRIPTS / "01_build_geometry.lsf").read_text()
    assert "run_spectra_driver;" in spectra_entry
    assert 'output_tag+"A_air"' in geometry_driver.read_text()
    assert 'output_tag+case_name' in spectra_driver.read_text()
    fast_geometry = (SCRIPTS / "01_build_geometry_fast.lsf").read_text()
    fast_spectra = (SCRIPTS / "02_run_spectra_fast.lsf").read_text()
    assert 'simulation_profile = "fast";' in fast_geometry
    assert 'simulation_profile = "fast";' in fast_spectra
    assert "build_geometry_driver;" in fast_geometry
    assert "run_spectra_driver;" in fast_spectra
    assert "01_build_geometry;" not in fast_geometry
    assert "02_run_spectra;" not in fast_spectra
    assert 'simulation_profile = "final";' in reference
    assert 'simulation_profile = "final";' in nearfield_run
    assert "NF_XY" not in spectra and "NF_XZ" not in spectra


def test_nearfield_plot_reads_lumerical_v73_mat_and_fast_wrapper(tmp_path):
    """Both field plot entry points accept Lumerical's HDF5 MAT exports."""
    normal = load_script("04_plot_nearfields")
    fast = load_script("04_plot_fast_nearfields")
    mat = tmp_path / "field_v73.mat"
    with h5py.File(mat, "w") as stream:
        stream["x_nm"] = [[-1, 1]]
        stream["y_nm"] = [[-1, 1]]
        stream["z_nm"] = [[-2, 0, 2]]
        stream["E2_xy"] = [[1, 3], [2, 4]]
        stream["E2_xz"] = [[10, 20], [30, 40], [50, 60]]
        stream["requested_wavelength_nm"] = [[800]]
        stream["sampled_wavelength_nm"] = [[805]]
    expected_xz = np.array([[10, 30, 50], [20, 40, 60]], dtype=float)
    for reader in (normal.read_fields, fast.read_fields):
        data = reader(mat)
        assert data["E2_xy"].shape == (2, 2)
        assert np.array_equal(data["E2_xz"], expected_xz)
        assert data["sampled_wavelength_nm"] == 805
