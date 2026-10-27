"""HFSS launch, model construction, and solver-result orchestration."""

import json
import os
from dataclasses import asdict
from pathlib import Path


AEDT_VERSION = "2025.2"
AEDT_ROOT_ENV = "ANSYSEM_ROOT252"
DEFAULT_AEDT_ROOT = Path("/usr/local/AnsysEM25R2/v252/AnsysEM")
VALID_CASE_NAMES = ("A_air", "B_SiO2")
VALID_ANGLES_DEG = (0, 15, 30, 45, 60, 75)


def validate_request(request, palik_table: Path | None) -> None:
    """Validate a solver request before AEDT is started.

    Args:
        request: A :class:`campaign.RunRequest` compatible object.
        palik_table: Exported Palik table required for a physical Au run.

    Raises:
        ValueError: If the requested case, angle, numerical controls, or
            required Palik data are invalid.
    """
    if request.case_name not in VALID_CASE_NAMES:
        raise ValueError(f"unknown case: {request.case_name}")
    if request.angle_deg not in VALID_ANGLES_DEG:
        raise ValueError(f"unsupported incidence angle: {request.angle_deg}")
    if request.mesh_nm <= 0 or request.port_offset_nm <= 0 or request.floquet_modes < 2:
        raise ValueError("mesh, port offset, and Floquet mode count must be positive")
    if request.has_gold and palik_table is None:
        raise ValueError("a physical Au run requires an exported Palik table")


def expected_excitation_count(request) -> int:
    """Return the number of modal excitations exposed by both Floquet ports.

    Args:
        request: A :class:`campaign.RunRequest` compatible object.

    Returns:
        Total source count used by PyAEDT design validation.
    """
    return 2 * request.floquet_modes


def audit_matches_request(audit: dict, request) -> bool:
    """Check whether persisted metadata describes exactly one requested run.

    Args:
        audit: Parsed build-audit JSON object.
        request: A :class:`campaign.RunRequest` compatible object.

    Returns:
        ``True`` only if all physical and numerical request fields match.
    """
    return audit.get("request") == asdict(request)


def has_converged_result(project_path: Path) -> bool:
    """Return whether HFSS recorded an explicitly converged solution.

    This recovery check is deliberately stricter than testing for a results
    directory. HFSS writes ``EMConvData`` only for the adaptive solution, and
    its ``IsConverged`` marker distinguishes a completed solve from a client
    that stopped while AEDT was still running.

    Args:
        project_path: Path to the AEDT project whose result tree is inspected.

    Returns:
        ``True`` if at least one adaptive convergence record states
        ``IsConverged=true``.
    """
    results_dir = project_path.with_name(f"{project_path.name}results")
    for convergence_file in results_dir.rglob("*.EMConvData"):
        try:
            contents = convergence_file.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if "IsConverged=true" in contents:
            return True
    return False


def is_completed_run(audit: dict, project_path: Path) -> bool:
    """Return whether audit metadata and the HFSS result agree on completion.

    Args:
        audit: Persisted build audit for one request.
        project_path: AEDT project whose adaptive convergence record is read.

    Returns:
        ``True`` only when both the audit and HFSS mark the case complete.
    """
    return bool(audit.get("solved")) and has_converged_result(project_path)


def recover_converged_audit(audit: dict, project_path: Path) -> bool:
    """Mark an audit as solved when HFSS completed after its client stopped.

    Args:
        audit: Mutable build-audit JSON object for an already validated run.
        project_path: AEDT project to associate with the audit.

    Returns:
        ``True`` when the audit was changed based on an HFSS convergence
        record, otherwise ``False``.
    """
    if audit.get("solved") or not has_converged_result(project_path):
        return False
    audit["solved"] = True
    audit["recovered_from_convergence_record"] = True
    return True


def persisted_solve_status(audit_path: Path) -> bool:
    """Read the solved status recorded in a build audit.

    Args:
        audit_path: JSON audit written beside an AEDT project.

    Returns:
        The persisted ``solved`` flag.

    Raises:
        ValueError: If the JSON audit is not an object.
    """
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if not isinstance(audit, dict):
        raise ValueError(f"audit must be a JSON object: {audit_path}")
    return bool(audit.get("solved"))


def solve_and_save(hfss, project_path: Path, setup_name: str = "Setup1", cores: int = 4) -> None:
    """Run a setup and save its result database before closing AEDT.

    Args:
        hfss: Open PyAEDT HFSS application object.
        project_path: AEDT project destination to persist after the solve.
        setup_name: Name of the Driven Modal setup to analyze.
        cores: Local CPU core count requested from HFSS.

    Raises:
        RuntimeError: If HFSS does not report a successful solve or save.
    """
    if not hfss.analyze_setup(setup_name, cores=cores):
        raise RuntimeError("HFSS solver did not report a successful analysis")
    if not hfss.save_project(str(project_path)):
        raise RuntimeError(f"could not save solved project: {project_path}")


def configure_aedt_environment() -> Path:
    """Validate and expose the local AEDT 2025.2 installation to PyAEDT.

    Returns:
        The AEDT installation directory containing the launcher.

    Raises:
        RuntimeError: If the configured installation lacks an executable
            ``ansysedt`` launcher.
    """
    root = Path(os.environ.get(AEDT_ROOT_ENV, DEFAULT_AEDT_ROOT)).expanduser()
    launcher = root / "ansysedt"
    if not launcher.is_file() or not os.access(launcher, os.X_OK):
        raise RuntimeError(f"HFSS launcher not found or executable: {launcher}")
    os.environ[AEDT_ROOT_ENV] = str(root)
    return root
