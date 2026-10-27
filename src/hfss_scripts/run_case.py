"""Command entry point for building or solving one auditable HFSS case."""

import json
import sys
from pathlib import Path

try:
    from .builder import build_project, project_path, project_stem, write_audit
    from .campaign import RunRequest
    from .material import read_palik_table
    from .runner import AEDT_VERSION, audit_matches_request, configure_aedt_environment, expected_excitation_count, is_completed_run, persisted_solve_status, recover_converged_audit, solve_and_save, validate_request
except ImportError:
    from builder import build_project, project_path, project_stem, write_audit
    from campaign import RunRequest
    from material import read_palik_table
    from runner import AEDT_VERSION, audit_matches_request, configure_aedt_environment, expected_excitation_count, is_completed_run, persisted_solve_status, recover_converged_audit, solve_and_save, validate_request


ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = ROOT / "results" / "hfss"


def run_case(request: RunRequest, solve: bool, palik_table: Path | None = None) -> Path:
    """Build and optionally solve one HFSS campaign request.

    Args:
        request: Required physical and numerical case definition.
        solve: Whether to run HFSS with four CPU cores after validation.
        palik_table: Validated exported Lumerical Palik data for Au runs.

    Returns:
        Path of the saved AEDT project.

    Raises:
        FileExistsError: If an existing project could be overwritten.
        RuntimeError: If AEDT construction, validation, saving, or solving fails.
    """
    validate_request(request, palik_table)
    material = read_palik_table(palik_table) if request.has_gold else None
    destination = project_path(RESULTS_ROOT, request)
    audit_path = destination.with_suffix(".build.json")
    existing_audit = None
    if destination.exists():
        if not audit_path.is_file():
            raise FileExistsError(f"existing project lacks audit metadata: {destination}")
        existing_audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if not audit_matches_request(existing_audit, request):
            raise FileExistsError(f"existing project does not match requested configuration: {destination}")
        if not existing_audit.get("validation_passed"):
            raise FileExistsError(f"existing project did not pass validation: {destination}")
        if is_completed_run(existing_audit, destination):
            return destination
        if recover_converged_audit(existing_audit, destination):
            write_audit(audit_path, existing_audit)
            return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    configure_aedt_environment()
    from ansys.aedt.core import Hfss

    hfss = Hfss(
        project=str(destination),
        design=project_stem(request),
        solution_type="Modal",
        version=AEDT_VERSION,
        non_graphical=True,
        new_desktop=True,
        close_on_exit=True,
    )
    try:
        if existing_audit is not None:
            if solve:
                solve_and_save(hfss, destination)
            existing_audit["solved"] = bool(existing_audit.get("solved") or solve)
            write_audit(audit_path, existing_audit)
            return destination
        audit = build_project(hfss, request, material)
        messages, valid = hfss.validate_full_design(
            output_dir=str(destination.parent),
            ports=expected_excitation_count(request),
        )
        audit["validation_messages"] = list(messages)
        audit["validation_passed"] = bool(valid)
        if not valid:
            write_audit(destination.with_suffix(".build.json"), audit)
            raise RuntimeError(f"HFSS design validation failed: {messages}")
        if not hfss.save_project(str(destination)):
            raise RuntimeError(f"could not save project: {destination}")
        if solve:
            solve_and_save(hfss, destination)
        audit["solved"] = solve
        write_audit(destination.with_suffix(".build.json"), audit)
        return destination
    finally:
        hfss.release_desktop(close_projects=True, close_desktop=True)


def main(argv: list[str]) -> int:
    """Parse the intentionally small command interface without argparse.

    Args:
        argv: Arguments after the Python executable.

    Returns:
        Zero on success and two on invalid input.
    """
    usage = "usage: run_case.py <reference|spectrum|sensitivity> <A_air|B_SiO2> <angle> [--solve] [--palik PATH]"
    if len(argv) < 3:
        print(usage, file=sys.stderr)
        return 2
    stage, case_name, angle_text, *options = argv
    try:
        angle = int(angle_text)
    except ValueError:
        print(usage, file=sys.stderr)
        return 2
    solve = "--solve" in options
    palik_table = None
    if "--palik" in options:
        index = options.index("--palik")
        if index + 1 == len(options):
            print(usage, file=sys.stderr)
            return 2
        palik_table = Path(options[index + 1])
    request = RunRequest(stage, case_name, angle, has_gold=stage != "reference")
    try:
        destination = run_case(request, solve=solve, palik_table=palik_table)
    except Exception as error:
        print(f"ERROR: {type(error).__name__}: {error}", file=sys.stderr)
        return 1
    print(json.dumps({"project": str(destination), "solved": persisted_solve_status(destination.with_suffix(".build.json"))}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
