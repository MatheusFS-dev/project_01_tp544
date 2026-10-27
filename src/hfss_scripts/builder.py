"""Create auditable HFSS Driven Modal unit-cell projects."""

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any


AEDT_VERSION = "2025.2"
PROJECT_SUBDIRECTORY = "projects"
SETUP_NAME = "Setup1"
SWEEP_NAME = "InterpolatingSweep"
SILICA_INDEX = 1.45


def silica_permittivity() -> float:
    """Return the assigned lossless SiO2 relative permittivity.

    Returns:
        Relative permittivity corresponding to refractive index 1.45.
    """
    return SILICA_INDEX**2


def project_stem(request: Any) -> str:
    """Return a stable HFSS filename stem for one campaign request.

    Args:
        request: A :class:`campaign.RunRequest` compatible object.

    Returns:
        Filename-safe stage, case, and angle identifier.
    """
    return f"{request.stage}_{request.case_name}_{request.angle_deg}deg"


def project_path(results_root: Path, request: Any) -> Path:
    """Return the isolated AEDT project location for one request.

    Args:
        results_root: Root directory reserved for HFSS outputs.
        request: A :class:`campaign.RunRequest` compatible object.

    Returns:
        Destination ``.aedt`` path below ``results/hfss/projects``.
    """
    return Path(results_root) / PROJECT_SUBDIRECTORY / f"{project_stem(request)}.aedt"


def _face_at(solid: Any, axis: int, coordinate_nm: float) -> int:
    """Locate the outer solid face closest to one coordinate plane."""
    face = min(solid.faces, key=lambda candidate: abs(candidate.center[axis] - coordinate_nm))
    return face.id


def _assign_lattice_pairs(hfss: Any, solid: Any, z_min_nm: float, z_max_nm: float, request: Any) -> None:
    """Apply matching scan-angle lattice pairs to one homogeneous layer."""
    half_period = 160.0
    theta = f"{request.angle_deg}deg"
    pairs = (
        ("X", _face_at(solid, 0, -half_period), _face_at(solid, 0, half_period)),
        ("Y", _face_at(solid, 1, -half_period), _face_at(solid, 1, half_period)),
    )
    for axis, face_low, face_high in pairs:
        boundary = hfss.assign_lattice_pair(
            [face_low, face_high],
            phase_delay="UseScanAngle",
            phase_delay_param1="0deg",
            phase_delay_param2=theta,
            name=f"Lattice_{axis}_{int(z_min_nm)}_{int(z_max_nm)}",
        )
        if not boundary:
            raise RuntimeError(f"could not assign {axis} lattice pair")


def _add_palik_material(hfss: Any, palik_data: Any) -> str:
    """Add Palik permittivity and positive bulk-loss datasets to HFSS."""
    real_dataset = hfss.create_dataset(
        "AuPalikEpsilonReal",
        x=palik_data.frequency_hz.tolist(),
        y=palik_data.epsilon_real.tolist(),
        x_unit="Hz",
    )
    loss_dataset = hfss.create_dataset(
        "AuPalikEpsilonImag",
        x=palik_data.frequency_hz.tolist(),
        y=palik_data.epsilon_imag.tolist(),
        x_unit="Hz",
    )
    if not real_dataset or not loss_dataset:
        raise RuntimeError("HFSS could not create Palik datasets")
    material_name = "Au_Palik_Imported"
    material = hfss.materials.add_material(material_name)
    if not material:
        raise RuntimeError("HFSS could not create Au material")
    material.permittivity = "pwl($AuPalikEpsilonReal, Freq)"
    material.conductivity = "2*pi*Freq*e0*pwl($AuPalikEpsilonImag, Freq)"
    if not material.update():
        raise RuntimeError("HFSS could not update the imported Au material")
    return material_name


def _add_silica_material(hfss: Any) -> str:
    """Add the assignment's explicit lossless SiO2 half-space material."""
    material_name = "SiO2_n1p45"
    material = hfss.materials.add_material(material_name)
    if not material:
        raise RuntimeError("HFSS could not create the SiO2 material")
    material.permittivity = silica_permittivity()
    material.dielectric_loss_tangent = 0.0
    if not material.update():
        raise RuntimeError("HFSS could not update the SiO2 material")
    return material_name


def build_project(hfss: Any, request: Any, palik_data: Any | None = None) -> dict[str, Any]:
    """Populate one HFSS unit-cell design without solving it.

    Args:
        hfss: Active PyAEDT ``Hfss`` design in Driven Modal mode.
        request: A validated campaign request.
        palik_data: Validated Palik samples for Au-containing models.

    Returns:
        Serializable audit metadata for the exact constructed model.

    Raises:
        RuntimeError: If a required AEDT construction operation fails.
    """
    hfss.modeler.model_units = "nm"
    half_period_nm = 160.0
    port_offset_nm = float(request.port_offset_nm)
    top = hfss.modeler.create_box(
        [-half_period_nm, -half_period_nm, 0],
        [320, 320, port_offset_nm],
        name="UpperMedium",
        material="vacuum",
    )
    lower_material = "vacuum" if request.case_name == "A_air" else _add_silica_material(hfss)
    bottom = hfss.modeler.create_box(
        [-half_period_nm, -half_period_nm, -port_offset_nm],
        [320, 320, port_offset_nm],
        name="LowerMedium",
        material=lower_material,
    )
    if not top or not bottom:
        raise RuntimeError("could not create the upper and lower medium volumes")

    _assign_lattice_pairs(hfss, top, 0, port_offset_nm, request)
    _assign_lattice_pairs(hfss, bottom, -port_offset_nm, 0, request)
    lattice_origin_top = [-half_period_nm, -half_period_nm, port_offset_nm]
    lattice_a_top = [half_period_nm, -half_period_nm, port_offset_nm]
    lattice_b_top = [-half_period_nm, half_period_nm, port_offset_nm]
    lattice_origin_bottom = [-half_period_nm, -half_period_nm, -port_offset_nm]
    lattice_a_bottom = [half_period_nm, -half_period_nm, -port_offset_nm]
    lattice_b_bottom = [-half_period_nm, half_period_nm, -port_offset_nm]
    top_port = hfss.create_floquet_port(
        _face_at(top, 2, port_offset_nm),
        lattice_origin=lattice_origin_top,
        lattice_a_end=lattice_a_top,
        lattice_b_end=lattice_b_top,
        modes=request.floquet_modes,
        name="TopFloquet",
    )
    bottom_port = hfss.create_floquet_port(
        _face_at(bottom, 2, -port_offset_nm),
        lattice_origin=lattice_origin_bottom,
        lattice_a_end=lattice_a_bottom,
        lattice_b_end=lattice_b_bottom,
        modes=request.floquet_modes,
        name="BottomFloquet",
    )
    if not top_port or not bottom_port:
        raise RuntimeError("could not create both Floquet ports")

    source_assignment = {}
    for mode in range(1, request.floquet_modes + 1):
        source_assignment[f"TopFloquet:{mode}"] = ("1W", "0deg") if mode == 2 else ("0W", "0deg")
        source_assignment[f"BottomFloquet:{mode}"] = ("0W", "0deg")
    if not hfss.edit_sources(source_assignment):
        raise RuntimeError("could not select the TM00 top-port excitation")

    material_hash = None
    if request.has_gold:
        if palik_data is None:
            raise RuntimeError("Au geometry requested without validated Palik data")
        gold_material = _add_palik_material(hfss, palik_data)
        outer = hfss.modeler.create_cylinder("Z", [0, 0, 0], 90, 30, name="AuRingOuter", material=gold_material)
        inner = hfss.modeler.create_cylinder("Z", [0, 0, 0], 50, 30, name="AuRingInner", material=gold_material)
        if not outer or not inner or not hfss.modeler.subtract(outer, inner, keep_originals=False):
            raise RuntimeError("could not create the Au nanoring")
        outer.name = "AuRing"
        mesh = hfss.mesh.assign_length_mesh(
            "AuRing",
            inside_selection=True,
            maximum_length=f"{request.mesh_nm}nm",
            maximum_elements=100_000,
            name="AuRingLocalMesh",
        )
        if not mesh:
            raise RuntimeError("could not assign local Au mesh")
        material_hash = palik_data.source_hash

    setup = hfss.create_setup(
        SETUP_NAME,
        setup_type="HFSSDriven",
        Frequency="374.7405725THz",
        MaximumPasses=20,
        MinimumPasses=2,
        MinimumConvergedPasses=2,
        MaxDeltaS=0.005,
    )
    if not setup:
        raise RuntimeError("could not create Driven Modal setup")
    sweep = hfss.create_linear_count_sweep(
        SETUP_NAME,
        "THz",
        187.37028625,
        666.20546222,
        num_of_freq_points=401,
        name=SWEEP_NAME,
        save_fields=True,
        sweep_type="Interpolating",
        interpolation_tol=0.5,
        interpolation_max_solutions=500,
    )
    if not sweep:
        raise RuntimeError("could not create interpolating frequency sweep")
    return {
        "aedt_version": AEDT_VERSION,
        "request": asdict(request),
        "ports": {"top": top_port.name, "bottom": bottom_port.name, "tm_source": "TopFloquet:2"},
        "material_hash": material_hash,
        "setup": SETUP_NAME,
        "sweep": SWEEP_NAME,
    }


def write_audit(path: Path, audit: dict[str, Any]) -> None:
    """Write reproducible build metadata next to an AEDT project.

    Args:
        path: Destination JSON file.
        audit: Serializable project construction information.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
