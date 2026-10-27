"""Explicit, resumable work manifest for the HFSS validation campaign."""

from dataclasses import dataclass


ANGLES_DEG = (0, 15, 30, 45, 60, 75)
CASE_NAMES = ("A_air", "B_SiO2")


@dataclass(frozen=True)
class RunRequest:
    """One solver request with its numerical configuration.

    Attributes:
        stage: Campaign phase that created the request.
        case_name: ``A_air`` or ``B_SiO2``.
        angle_deg: TM incidence angle from the surface normal.
        mesh_nm: Local maximum element length in the Au ring.
        port_offset_nm: Floquet-port distance from the ring plane.
        floquet_modes: Number of modes retained on each port.
        has_gold: Whether the physical Au ring is included.
    """

    stage: str
    case_name: str
    angle_deg: int
    mesh_nm: int = 3
    port_offset_nm: int = 800
    floquet_modes: int = 10
    has_gold: bool = True


@dataclass(frozen=True)
class CampaignManifest:
    """Grouped complete list of production and verification requests."""

    spectra: tuple[RunRequest, ...]
    references: tuple[RunRequest, ...]
    sensitivities: tuple[RunRequest, ...]


def campaign_manifest() -> CampaignManifest:
    """Build the fixed manifest from the assignment and validation plan.

    Returns:
        Twelve physical spectra, twelve no-Au references, and seven
        deliberately stricter sensitivity runs.
    """
    spectra = tuple(
        RunRequest(
            "spectrum",
            case_name,
            angle,
            floquet_modes=14 if case_name == "B_SiO2" else 10,
        )
        for case_name in CASE_NAMES
        for angle in ANGLES_DEG
    )
    references = tuple(
        RunRequest("reference", case_name, angle, has_gold=False)
        for case_name in CASE_NAMES
        for angle in ANGLES_DEG
    )
    sensitivity_pairs = (
        ("A_air", 0),
        ("A_air", 30),
        ("A_air", 45),
        ("A_air", 75),
        ("B_SiO2", 0),
        ("B_SiO2", 45),
        ("B_SiO2", 75),
    )
    sensitivities = tuple(
        RunRequest(
            "sensitivity",
            case_name,
            angle,
            mesh_nm=2,
            port_offset_nm=1200,
            floquet_modes=14,
        )
        for case_name, angle in sensitivity_pairs
    )
    return CampaignManifest(spectra=spectra, references=references, sensitivities=sensitivities)
