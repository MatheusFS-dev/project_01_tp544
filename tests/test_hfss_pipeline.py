"""Regression tests for the HFSS nanoring workflow core."""

import hashlib
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "src" / "hfss_scripts"


def load_module(name):
    """Load an HFSS workflow module without requiring package installation."""
    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_configuration_matches_assignment_and_frequency_limits():
    """The fixed nanoring configuration is encoded once in SI units."""
    config = load_module("config")
    model = config.ModelConfig()

    assert model.outer_radius_m == pytest.approx(90e-9)
    assert model.inner_radius_m == pytest.approx(50e-9)
    assert model.thickness_m == pytest.approx(30e-9)
    assert model.period_m == pytest.approx(320e-9)
    assert model.angles_deg == (0, 15, 30, 45, 60, 75)
    assert model.wavelength_nm == (450, 1600)
    assert model.frequency_thz == pytest.approx((187.37028625, 666.20546222))
    assert model.case_names == ("A_air", "B_SiO2")


def test_palik_table_is_validated_and_hashed(tmp_path):
    """Palik import uses positive-loss samples with coverage over the full band."""
    material = load_module("material")
    table = tmp_path / "Au_Palik.tsv"
    table.write_text(
        "frequency_hz\tepsilon_real\tepsilon_imag\n"
        "1.8e14\t-10\t1.0\n"
        "3.0e14\t-8\t1.5\n"
        "6.8e14\t-2\t2.0\n",
        encoding="utf-8",
    )

    data = material.read_palik_table(table)

    assert data.source_hash == hashlib.sha256(table.read_bytes()).hexdigest()
    assert np.array_equal(data.frequency_hz, [1.8e14, 3e14, 6.8e14])
    assert np.array_equal(data.epsilon_real, [-10, -8, -2])
    assert np.array_equal(data.epsilon_imag, [1, 1.5, 2])
    assert np.all(data.conductivity_s_per_m > 0)


@pytest.mark.parametrize(
    "contents",
    [
        "frequency_hz\tepsilon_real\tepsilon_imag\n1.8e14\t-2\t1\n6.8e14\t-1\t-0.2\n",
        "frequency_hz\tepsilon_real\tepsilon_imag\n3e14\t-2\t1\n1.8e14\t-1\t1\n6.8e14\t-1\t1\n",
        "frequency_hz\tepsilon_real\tepsilon_imag\n3e14\t-2\t1\n6e14\t-1\t1\n",
    ],
)
def test_palik_table_rejects_nonpassive_unsorted_or_incomplete_data(tmp_path, contents):
    """The imported table must be monotonic, passive, and span 450 to 1600 nm."""
    material = load_module("material")
    table = tmp_path / "Au_Palik.tsv"
    table.write_text(contents, encoding="utf-8")

    with pytest.raises(ValueError):
        material.read_palik_table(table)


def test_fresnel_tm_reference_and_resampling_track_interpolation():
    """Reference powers and generated 1 nm samples use explicit interpolation flags."""
    analysis = load_module("analysis")

    air = analysis.fresnel_tm(1.0, 1.0, 45)
    silica = analysis.fresnel_tm(1.0, 1.45, 0)
    grid = analysis.resample_spectrum(
        np.array([450.0, 500.0, 550.0]),
        np.array([0.1, 0.2, 0.3]),
        np.array([0.9, 0.8, 0.7]),
        np.array([450.0, 475.0, 500.0]),
    )

    assert air.reflection == pytest.approx(0.0)
    assert air.transmission == pytest.approx(1.0)
    assert silica.reflection == pytest.approx(((1.45 - 1) / (1.45 + 1)) ** 2)
    assert silica.transmission == pytest.approx(1 - silica.reflection)
    assert np.allclose(grid.transmission, [0.1, 0.15, 0.2])
    assert grid.interpolated.tolist() == [False, True, False]


def test_resampling_accepts_hfss_wavelength_roundoff_at_requested_endpoint():
    """Frequency-to-wavelength roundoff does not reject a physical endpoint."""
    analysis = load_module("analysis")

    grid = analysis.resample_spectrum(
        np.array([450.000000002, 1600.0]),
        np.array([0.1, 0.2]),
        np.array([0.3, 0.4]),
        np.array([450.0, 1600.0]),
    )

    assert grid.interpolated.tolist() == [False, False]


def test_specular_s_parameters_produce_unclipped_tra_spectrum():
    """Power-normalized Floquet S parameters give T, R, and A directly."""
    analysis = load_module("analysis")

    spectrum = analysis.specular_tra_spectrum(
        np.array([2.0e14, 3.0e14]),
        np.array([0.3 + 0.4j, 0.8 + 0j]),
        np.array([0.5 + 0j, 0.7 + 0j]),
    )

    assert spectrum.reflection.tolist() == pytest.approx([0.25, 0.64])
    assert spectrum.transmission.tolist() == pytest.approx([0.25, 0.49])
    assert spectrum.absorption.tolist() == pytest.approx([0.5, -0.13])
    assert spectrum.wavelength_nm.tolist() == pytest.approx([1498.96229, 999.308193])


def test_propagating_modal_sum_excludes_evanescent_floquet_orders():
    """T and R sum only modal powers whose HFSS gamma has no attenuation."""
    analysis = load_module("analysis")
    frequencies = np.array([2.0e14, 3.0e14])
    reflected_s = np.array([[0.2 + 0j, 0.3 + 0j], [9 + 0j, 9 + 0j]])
    transmitted_s = np.array([[0.8 + 0j, 0.7 + 0j], [8 + 0j, 8 + 0j]])
    gamma = np.array([[0 + 4j, 0 + 5j], [3 + 0j, 2 + 0j]])

    spectrum = analysis.propagating_modal_tra_spectrum(
        frequencies,
        reflected_s,
        transmitted_s,
        gamma,
        gamma,
    )

    assert spectrum.reflection.tolist() == pytest.approx([0.04, 0.09])
    assert spectrum.transmission.tolist() == pytest.approx([0.64, 0.49])
    assert spectrum.absorption.tolist() == pytest.approx([0.32, 0.42])


def test_modal_sum_rejects_a_propagating_highest_included_mode():
    """A propagating final mode means the Floquet truncation is inconclusive."""
    analysis = load_module("analysis")
    gamma = np.array([[0 + 4j, 0 + 5j], [0 + 3j, 0 + 2j]])

    with pytest.raises(ValueError, match="highest included"):
        analysis.propagating_modal_tra_spectrum(
            np.array([2.0e14, 3.0e14]),
            np.zeros((2, 2), dtype=complex),
            np.zeros((2, 2), dtype=complex),
            gamma,
            gamma,
        )


def test_air_normal_export_resamples_and_labels_non_native_points(tmp_path):
    """The report grid has 1151 samples and retains interpolation provenance."""
    analysis = load_module("analysis")
    exporter = load_module("export_air_normal")
    spectrum = analysis.Traspectrum(
        frequency_hz=np.array([299_792_458 / 1600e-9, 299_792_458 / 450e-9]),
        wavelength_nm=np.array([1600.0, 450.0]),
        transmission=np.array([0.8, 0.6]),
        reflection=np.array([0.1, 0.2]),
        absorption=np.array([0.1, 0.2]),
    )

    resampled, interpolated = exporter.resampled_spectrum(spectrum)
    output = tmp_path / "spectrum.csv"
    exporter.write_spectrum_csv(output, resampled, interpolated)

    rows = output.read_text(encoding="utf-8").splitlines()
    assert len(resampled.wavelength_nm) == 1151
    assert rows[1].endswith(",false")
    assert rows[2].endswith(",true")
    assert rows[-1].endswith(",false")


def test_modal_export_requests_every_port_mode_from_the_tm_source():
    """Power export includes every retained Floquet channel before filtering."""
    exporter = load_module("export_modal_spectra")

    expressions = exporter.modal_s_parameter_expressions("BottomFloquet", 3)

    assert expressions == [
        "S(BottomFloquet:1,TopFloquet:2)",
        "S(BottomFloquet:2,TopFloquet:2)",
        "S(BottomFloquet:3,TopFloquet:2)",
    ]


def test_modal_export_writes_lumerical_style_tra_table(tmp_path):
    """HFSS output is consumable as a 1 nm T/R/A table with no clipped A."""
    analysis = load_module("analysis")
    exporter = load_module("export_modal_spectra")
    spectrum = analysis.Traspectrum(
        frequency_hz=np.array([299_792_458 / 450e-9, 299_792_458 / 451e-9]),
        wavelength_nm=np.array([450.0, 451.0]),
        transmission=np.array([0.7, 0.8]),
        reflection=np.array([0.1, 0.3]),
        absorption=np.array([0.2, -0.1]),
    )
    destination = tmp_path / "A_air_0deg.tsv"

    exporter.write_lumerical_style_table(destination, spectrum)

    assert destination.read_text(encoding="utf-8").splitlines() == [
        "wavelength_nm\tT\tR\tA",
        "450\t0.7\t0.1\t0.2",
        "451\t0.8\t0.3\t-0.1",
    ]


def test_resonance_detection_labels_prominent_absorption_peaks_by_wavelength():
    """Automated candidates retain their peak position and prominence."""
    review = load_module("review_hfss_results")
    wavelength_nm = np.arange(450.0, 457.0)
    absorption = np.array([0.1, 0.2, 0.6, 0.2, 0.15, 0.5, 0.1])

    candidates = review.absorption_resonances(wavelength_nm, absorption, prominence=0.1)

    assert [(item.identifier, item.wavelength_nm) for item in candidates] == [
        ("A01", 452.0),
        ("A02", 455.0),
    ]
    assert [item.prominence for item in candidates] == pytest.approx([0.5, 0.35])


def test_hfss_field_parser_reconstructs_a_cartesian_plane_from_numeric_rows(tmp_path):
    """HFSS grid exports retain coordinates and scalar field magnitude."""
    fields = load_module("export_hfss_fields")
    path = tmp_path / "plane.fld"
    path.write_text(
        "# HFSS field export\n"
        "0 0 3.5e-8 1.5\n"
        "1e-9 0 3.5e-8 2.5\n"
        "0 1e-9 3.5e-8 3.5\n"
        "1e-9 1e-9 3.5e-8 4.5\n",
        encoding="utf-8",
    )

    x, y, values = fields.read_cartesian_plane(path, horizontal_axis=0, vertical_axis=1)

    assert x.tolist() == pytest.approx([0.0, 1.0])
    assert y.tolist() == pytest.approx([0.0, 1.0])
    assert np.allclose(values, [[1.5, 3.5], [2.5, 4.5]])


def test_hfss_field_parser_preserves_masked_nan_cells_inside_metal(tmp_path):
    """HFSS field grids may mask metal interior samples as NaN."""
    fields = load_module("export_hfss_fields")
    path = tmp_path / "plane.fld"
    path.write_text(
        "0 0 0 1\n1 0 0 Nan\n0 1 0 3\n1 1 0 4\n",
        encoding="utf-8",
    )

    _, _, values = fields.read_cartesian_plane(path, horizontal_axis=0, vertical_axis=1)

    assert np.isnan(values[1, 0])


def test_field_figure_paths_preserve_xy_xz_suffixes_with_decimal_wavelengths(tmp_path):
    """XY and XZ figures do not overwrite each other for decimal stems."""
    fields = load_module("export_hfss_fields")

    xy_png, xy_pdf = fields.field_figure_paths(tmp_path / "B_SiO2_0deg_A01_483.0nm_xy")
    xz_png, xz_pdf = fields.field_figure_paths(tmp_path / "B_SiO2_0deg_A01_483.0nm_xz")

    assert xy_png.name == "B_SiO2_0deg_A01_483.0nm_xy.png"
    assert xy_pdf.name == "B_SiO2_0deg_A01_483.0nm_xy.pdf"
    assert xz_png.name == "B_SiO2_0deg_A01_483.0nm_xz.png"
    assert xz_pdf.name == "B_SiO2_0deg_A01_483.0nm_xz.pdf"


def test_field_intensity_squares_normalized_magnitude_without_filling_nan():
    """HFSS maps use the Lumerical normalized-intensity convention."""
    fields = load_module("export_hfss_fields")

    intensity = fields.normalized_intensity(np.array([[2.0, np.nan]]))

    assert intensity[0, 0] == pytest.approx(4.0)
    assert np.isnan(intensity[0, 1])


def test_field_export_uses_the_documented_vector_electric_field(tmp_path):
    """Grid export derives magnitude from the vector E field calculator."""
    fields = load_module("export_hfss_fields")
    captured = {}
    output = tmp_path / "field.fld"

    class FakePost:
        def export_field_file_on_grid(self, **kwargs):
            captured.update(kwargs)
            Path(kwargs["file_name"]).write_text("0 0 0 1\n", encoding="utf-8")
            return kwargs["file_name"]

    class FakeHfss:
        post = FakePost()

    fields.export_plane(
        FakeHfss(),
        output,
        [-160, -160, 35],
        [160, 160, 35],
        [2, 2, 0],
        374.7405725e12,
    )

    assert captured["quantity"] == "E"
    assert captured["is_vector"] is True
    assert captured["intrinsics"]["Freq"] == "374740572500000Hz"


def test_field_normalization_closes_the_physical_project_before_opening_reference(monkeypatch, tmp_path):
    """Field export avoids two project designs sharing one AEDT desktop."""
    fields = load_module("export_hfss_fields")
    events = []

    class FakeHfss:
        def __init__(self, name):
            self.name = name

        def release_desktop(self, **kwargs):
            events.append(f"close:{self.name}")

    def fake_open(request):
        events.append(f"open:{request.case_name}")
        return FakeHfss(request.case_name)

    def fake_export(hfss, output_path, *args):
        events.append(f"export:{hfss.name}:{output_path.name}")

    def fake_read(path, horizontal_axis, vertical_axis):
        coordinates = np.array([0.0, 1.0])
        values = np.full((2, 2), 2.0 if "physical" in path.name else 1.0)
        return coordinates, coordinates, values

    monkeypatch.setattr(fields, "open_hfss", fake_open)
    monkeypatch.setattr(fields, "export_plane", fake_export)
    monkeypatch.setattr(fields, "read_cartesian_plane", fake_read)
    campaign = load_module("campaign")
    physical = campaign.RunRequest("spectrum", "A_air", 0)
    reference = campaign.RunRequest("reference", "A_air", 0, has_gold=False)

    _, _, xy_field, _, _, xz_field = fields.normalized_planes(physical, reference, 800.0, tmp_path / "map")

    assert events.index("close:A_air") < events.index("open:A_air", 1)
    assert np.allclose(xy_field, 2.0)
    assert np.allclose(xz_field, 2.0)


def test_field_sample_selection_uses_the_nearest_saved_hfss_basis_frequency(tmp_path):
    """Near fields use a solved sweep point instead of an interpolated point."""
    fields = load_module("export_hfss_fields")
    results = tmp_path / "case.results"
    results.mkdir()
    (results / "DV_F1_SU.txt").write_text("Frequency 5.99584916e14\n", encoding="utf-8")
    (results / "DV_F2_SU.txt").write_text("Frequency 5.87828349e14\n", encoding="utf-8")
    (results / "DV_F3_SU.txt").write_text("Frequency 5.76523958e14\n", encoding="utf-8")

    wavelength_nm, frequency_hz = fields.nearest_saved_field_sample(results, 514.0)

    assert wavelength_nm == pytest.approx(510.0)
    assert frequency_hz == pytest.approx(5.87828349e14)


def test_direct_field_sweep_inputs_do_not_reuse_a_mutable_frequency_list():
    """Each physical and reference field sweep receives an independent list."""
    direct_fields = load_module("run_direct_fields")
    frequencies = [1.0, 2.0]

    first, first_flags = direct_fields.field_sweep_inputs(frequencies)
    first.pop()
    second, second_flags = direct_fields.field_sweep_inputs(frequencies)

    assert frequencies == [1.0, 2.0]
    assert second == [1.0, 2.0]
    assert first_flags == [True, True]
    assert second_flags == [True, True]


def test_campaign_manifest_covers_all_required_cases_and_sensitivities():
    """The manifest enumerates every requested spectrum and refinement check."""
    campaign = load_module("campaign")
    manifest = campaign.campaign_manifest()

    assert len(manifest.spectra) == 12
    assert {(item.case_name, item.angle_deg) for item in manifest.spectra} == {
        (case_name, angle)
        for case_name in ("A_air", "B_SiO2")
        for angle in (0, 15, 30, 45, 60, 75)
    }
    assert len(manifest.references) == 12
    assert {(item.case_name, item.angle_deg) for item in manifest.sensitivities} == {
        ("A_air", 0),
        ("A_air", 30),
        ("A_air", 45),
        ("A_air", 75),
        ("B_SiO2", 0),
        ("B_SiO2", 45),
        ("B_SiO2", 75),
    }
    assert all(item.mesh_nm == 2 for item in manifest.sensitivities)
    assert all(item.port_offset_nm == 1200 for item in manifest.sensitivities)


def test_silica_baseline_spectra_retain_extra_floquet_modes():
    """SiO2 spectra retain an evanescent mode beyond the diffraction orders."""
    campaign = load_module("campaign")
    manifest = campaign.campaign_manifest()

    air_modes = {item.floquet_modes for item in manifest.spectra if item.case_name == "A_air"}
    silica_modes = {item.floquet_modes for item in manifest.spectra if item.case_name == "B_SiO2"}

    assert air_modes == {10}
    assert silica_modes == {14}


def test_background_campaign_runs_references_before_physical_and_sensitivity_cases():
    """The unattended CPU campaign preserves the validation-first sequence."""
    campaign_runner = load_module("run_campaign")
    requests = campaign_runner.campaign_requests()

    assert len(requests) == 31
    assert all(request.stage == "reference" for request in requests[:12])
    assert all(request.stage == "spectrum" for request in requests[12:24])
    assert all(request.stage == "sensitivity" for request in requests[24:])


def test_solve_and_save_persists_the_aedt_project_after_analysis(tmp_path):
    """A completed solve is saved before the AEDT session can be closed."""
    runner = load_module("runner")
    calls = []

    class FakeHfss:
        def analyze_setup(self, setup_name, cores):
            calls.append(("analyze", setup_name, cores))
            return True

        def save_project(self, destination):
            calls.append(("save", destination))
            return True

    destination = tmp_path / "case.aedt"
    runner.solve_and_save(FakeHfss(), destination)

    assert calls == [("analyze", "Setup1", 4), ("save", str(destination))]


def test_gold_runs_require_a_valid_palik_table_but_references_do_not(tmp_path):
    """A physical run cannot silently use a non-Palik or missing Au model."""
    campaign = load_module("campaign")
    runner = load_module("runner")
    gold_request = campaign.RunRequest("spectrum", "A_air", 0)
    reference_request = campaign.RunRequest("reference", "A_air", 0, has_gold=False)

    with pytest.raises(ValueError, match="Palik"):
        runner.validate_request(gold_request, None)
    runner.validate_request(reference_request, None)
    with pytest.raises(ValueError, match="case"):
        runner.validate_request(campaign.RunRequest("reference", "other", 0, has_gold=False), None)


def test_hfss_project_names_are_stable_and_do_not_overwrite_lumerical_results(tmp_path):
    """Each HFSS request has an isolated, deterministic AEDT project path."""
    campaign = load_module("campaign")
    builder = load_module("builder")
    request = campaign.RunRequest("reference", "B_SiO2", 45, has_gold=False)

    assert builder.project_stem(request) == "reference_B_SiO2_45deg"
    assert builder.project_path(tmp_path, request) == tmp_path / "projects" / "reference_B_SiO2_45deg.aedt"


def test_hfss_validation_counts_floquet_modes_not_physical_port_surfaces():
    """Two ten-mode Floquet ports expose twenty HFSS excitations."""
    campaign = load_module("campaign")
    runner = load_module("runner")

    assert runner.expected_excitation_count(campaign.RunRequest("reference", "A_air", 0, has_gold=False)) == 20


def test_resume_rejects_a_project_created_for_a_different_request():
    """A saved AEDT project is resumed only when its recorded request matches."""
    campaign = load_module("campaign")
    runner = load_module("runner")
    request = campaign.RunRequest("reference", "A_air", 0, has_gold=False)
    matching_audit = {"request": {"stage": "reference", "case_name": "A_air", "angle_deg": 0,
                                  "mesh_nm": 3, "port_offset_nm": 800, "floquet_modes": 10,
                                  "has_gold": False}}

    assert runner.audit_matches_request(matching_audit, request)
    assert not runner.audit_matches_request(matching_audit, campaign.RunRequest("reference", "A_air", 45, has_gold=False))


def test_converged_hfss_result_is_detected_after_an_interrupted_client(tmp_path):
    """A completed solver result can be recovered without another solve."""
    runner = load_module("runner")
    project = tmp_path / "spectrum_A_air_0deg.aedt"
    result = tmp_path / "spectrum_A_air_0deg.aedtresults" / "run.results" / "run.EMConvData"
    result.parent.mkdir(parents=True)
    result.write_text("$begin 'ConvInfo'\nIsConverged=true\n$end 'ConvInfo'\n", encoding="utf-8")

    assert runner.has_converged_result(project)


def test_unconverged_or_missing_hfss_result_is_not_recovered(tmp_path):
    """Only an explicit HFSS convergence marker authorizes a resume shortcut."""
    runner = load_module("runner")
    project = tmp_path / "spectrum_A_air_0deg.aedt"
    result = tmp_path / "spectrum_A_air_0deg.aedtresults" / "run.results" / "run.EMConvData"
    result.parent.mkdir(parents=True)
    result.write_text("$begin 'ConvInfo'\nIsConverged=false\n$end 'ConvInfo'\n", encoding="utf-8")

    assert not runner.has_converged_result(project)
    assert not runner.has_converged_result(tmp_path / "missing.aedt")


def test_recovery_marks_only_a_matching_converged_audit_as_solved(tmp_path):
    """Recovery persists a provenance flag instead of launching a duplicate solve."""
    runner = load_module("runner")
    project = tmp_path / "spectrum_A_air_0deg.aedt"
    result = tmp_path / "spectrum_A_air_0deg.aedtresults" / "run.EMConvData"
    result.parent.mkdir(parents=True)
    result.write_text("IsConverged=true\n", encoding="utf-8")
    audit = {"solved": False}

    assert runner.recover_converged_audit(audit, project)
    assert audit == {"solved": True, "recovered_from_convergence_record": True}


def test_completed_audit_with_convergence_record_skips_a_duplicate_solve(tmp_path):
    """A saved, converged request is complete even when --solve is requested."""
    runner = load_module("runner")
    project = tmp_path / "spectrum_A_air_0deg.aedt"
    convergence = tmp_path / "spectrum_A_air_0deg.aedtresults" / "run.EMConvData"
    convergence.parent.mkdir(parents=True)
    convergence.write_text("IsConverged=true\n", encoding="utf-8")

    assert runner.is_completed_run({"solved": True}, project)
    assert not runner.is_completed_run({"solved": False}, project)


def test_run_case_status_reads_the_persisted_audit_flag(tmp_path):
    """The command status reports a recovered solve rather than CLI intent."""
    runner = load_module("runner")
    audit_path = tmp_path / "case.build.json"
    audit_path.write_text('{"solved": true}', encoding="utf-8")

    assert runner.persisted_solve_status(audit_path)


def test_lumerical_export_script_preserves_palik_sampled_permittivity():
    """The bridge exports the database samples, not a fitted replacement model."""
    script = (SCRIPTS / "export_palik_from_lumerical.lsf").read_text(encoding="utf-8")

    assert 'material_name = "Au (Gold) - Palik"' in script
    assert 'getmaterial(material_name,"sampled data")' in script
    assert 'real(epsilon)' in script
    assert 'imag(epsilon)' in script
    assert "frequency_hz\\tepsilon_real\\tepsilon_imag" in script


def test_silica_material_uses_the_required_n_approximately_1p45():
    """The substrate model is an explicit lossless n=1.45 bulk material."""
    builder = load_module("builder")

    assert builder.silica_permittivity() == pytest.approx(1.45**2)


def test_external_palik_parser_converts_angstrom_nk_to_passive_permittivity():
    """The approved external Palik table is converted without spectral fitting."""
    external = load_module("external_palik")
    raw = """; Optical constants for polycrystalline Au\n; Lambda (A) n k\n4500 1.5 1.8\n8000 0.2 5.0\n16000 0.5 10.0\n"""

    samples = external.parse_palik_nk(raw)

    assert np.all(np.diff(samples.frequency_hz) > 0)
    assert samples.frequency_hz[0] == pytest.approx(299_792_458 / 1600e-9)
    assert samples.frequency_hz[-1] == pytest.approx(299_792_458 / 450e-9)
    assert np.all(samples.epsilon_imag > 0)
    assert samples.epsilon_real[-1] == pytest.approx(1.5**2 - 1.8**2)


def test_external_palik_parser_retains_first_published_value_at_duplicate_wavelength():
    """Duplicated source wavelengths are resolved deterministically without fitting."""
    external = load_module("external_palik")
    raw = """4500 1.5 1.8\n8000 0.2 5.0\n8000 0.4 6.0\n16000 0.5 10.0\n"""

    samples = external.parse_palik_nk(raw)

    assert len(samples.frequency_hz) == 3
    assert samples.epsilon_imag[1] == pytest.approx(2 * 0.2 * 5.0)
