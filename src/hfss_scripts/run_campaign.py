"""Run the complete HFSS campaign sequentially in a durable process."""

import json
from datetime import datetime, timezone
from pathlib import Path

try:
    from .campaign import RunRequest, campaign_manifest
    from .run_case import RESULTS_ROOT, run_case
except ImportError:
    from campaign import RunRequest, campaign_manifest
    from run_case import RESULTS_ROOT, run_case


PALIK_TABLE = RESULTS_ROOT / "Au_Palik_external.tsv"
LOG_PATH = RESULTS_ROOT / "campaign_status.jsonl"


def campaign_requests() -> tuple[RunRequest, ...]:
    """Return the validation-first execution order for the complete campaign.

    Returns:
        No-Au references, then physical spectra, then sensitivity cases.
    """
    manifest = campaign_manifest()
    return manifest.references + manifest.spectra + manifest.sensitivities


def write_status(request: RunRequest, state: str, detail: str = "") -> None:
    """Append one machine-readable campaign status event.

    Args:
        request: Case being started, completed, or failed.
        state: Lifecycle state for the request.
        detail: Optional diagnostic text.
    """
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    event = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "request": request.__dict__,
        "state": state,
        "detail": detail,
    }
    with LOG_PATH.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(event, sort_keys=True) + "\n")


def run_campaign() -> None:
    """Solve every required validation, physical, and sensitivity request.

    Raises:
        RuntimeError: If a case cannot be built, solved, or saved. The status
            log identifies the failed request and no later cases are started.
    """
    for request in campaign_requests():
        write_status(request, "started")
        try:
            palik_table = PALIK_TABLE if request.has_gold else None
            destination = run_case(request, solve=True, palik_table=palik_table)
        except Exception as error:
            write_status(request, "failed", f"{type(error).__name__}: {error}")
            raise
        write_status(request, "completed", str(destination))


if __name__ == "__main__":
    run_campaign()
