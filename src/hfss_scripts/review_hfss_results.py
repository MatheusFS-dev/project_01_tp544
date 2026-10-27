"""Detect resonance candidates and plot completed HFSS T/R/A spectra."""

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

try:
    from .campaign import campaign_manifest
    from .publication_plot_style import PALETTE, create_figure, save_figure, style_axes, style_legend
except ImportError:
    from campaign import campaign_manifest
    from publication_plot_style import PALETTE, create_figure, save_figure, style_axes, style_legend


ROOT = Path(__file__).resolve().parents[2]
SPECTRA_ROOT = ROOT / "results" / "hfss" / "spectra"
RESONANCES_ROOT = ROOT / "results" / "hfss" / "resonances"
PROMINENCE = 0.02


@dataclass(frozen=True)
class ResonanceCandidate:
    """One automatically detected local maximum in absorption."""

    identifier: str
    wavelength_nm: float
    absorption: float
    prominence: float


def absorption_resonances(
    wavelength_nm: np.ndarray,
    absorption: np.ndarray,
    prominence: float = PROMINENCE,
) -> list[ResonanceCandidate]:
    """Detect prominent absorption maxima on an ascending wavelength grid.

    Candidate labels are rank identifiers within one spectrum. They are not a
    claim that peaks at different angles share a common physical mode.

    Args:
        wavelength_nm: Strictly increasing wavelength samples in nanometers.
        absorption: Absorption samples associated with ``wavelength_nm``.
        prominence: Minimum absolute absorption prominence for a candidate.

    Returns:
        Candidates ordered by wavelength.

    Raises:
        ValueError: If the grid, absorption data, or prominence is invalid.
    """
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    absorption = np.asarray(absorption, dtype=float)
    if wavelength_nm.ndim != 1 or absorption.shape != wavelength_nm.shape or len(wavelength_nm) < 3:
        raise ValueError("wavelength and absorption arrays must share at least three samples")
    if not np.all(np.isfinite(wavelength_nm)) or not np.all(np.isfinite(absorption)):
        raise ValueError("resonance input must be finite")
    if np.any(np.diff(wavelength_nm) <= 0) or prominence <= 0:
        raise ValueError("wavelength must increase and prominence must be positive")
    indices, properties = find_peaks(absorption, prominence=prominence)
    return [
        ResonanceCandidate(
            identifier=f"A{rank:02d}",
            wavelength_nm=float(wavelength_nm[index]),
            absorption=float(absorption[index]),
            prominence=float(properties["prominences"][rank - 1]),
        )
        for rank, index in enumerate(indices, 1)
    ]


def read_spectrum(path: Path) -> np.ndarray:
    """Read a complete one-nanometer HFSS T/R/A table.

    Args:
        path: TSV table produced by :mod:`export_modal_spectra`.

    Returns:
        Array with wavelength, transmission, reflection, and absorption.

    Raises:
        ValueError: If the table is incomplete, malformed, or fails power
            balance within serialization precision.
    """
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["wavelength_nm", "T", "R", "A"]:
            raise ValueError(f"invalid spectrum columns: {path}")
        rows = [[float(row[field]) for field in reader.fieldnames] for row in reader]
    data = np.asarray(rows, dtype=float)
    if data.shape != (1151, 4) or not np.all(np.isfinite(data)):
        raise ValueError(f"incomplete or nonfinite spectrum: {path}")
    data = data[np.argsort(data[:, 0])]
    if not np.allclose(data[:, 0], np.arange(450.0, 1601.0), atol=1e-6, rtol=0):
        raise ValueError(f"invalid wavelength grid: {path}")
    if not np.allclose(data[:, 1:].sum(axis=1), 1, atol=1e-8, rtol=0):
        raise ValueError(f"power balance failure: {path}")
    return data


def plot_spectrum(data: np.ndarray, candidates: list[ResonanceCandidate], stem: Path, title: str) -> None:
    """Save matching PNG and PDF T/R/A plots with candidate markers.

    Args:
        data: Validated wavelength, T, R, A table.
        candidates: Detected absorption candidates to mark.
        stem: Output path without a suffix.
        title: Figure title.
    """
    figure, axes = create_figure(width=7.4, height=4.6)
    for column, label, color in ((1, "T", PALETTE["primary"]), (2, "R", PALETTE["highlight"]), (3, "A", PALETTE["secondary"])):
        axes.plot(data[:, 0], data[:, column], color=color, linewidth=1.35, label=label)
    for candidate in candidates:
        axes.axvline(candidate.wavelength_nm, color=PALETTE["muted"], linewidth=0.7, linestyle="--")
        axes.scatter(candidate.wavelength_nm, candidate.absorption, color=PALETTE["selection"], s=22, zorder=3)
    axes.set_xlim(450, 1600)
    style_axes(axes, "Wavelength (nm)", "Power fraction", title=title, grid_axis="y")
    style_legend(axes, location="best")
    save_figure(figure, stem.with_suffix(".png"))
    save_figure(figure, stem.with_suffix(".pdf"))
    plt.close(figure)


def review_all() -> None:
    """Create plots and a resonance-candidate table for all physical spectra.

    Raises:
        FileNotFoundError: If a required completed spectrum is missing.
    """
    RESONANCES_ROOT.mkdir(parents=True, exist_ok=True)
    rows = []
    for request in campaign_manifest().spectra:
        path = SPECTRA_ROOT / f"{request.case_name}_{request.angle_deg}deg.tsv"
        if not path.is_file() or not path.with_suffix(".done").is_file():
            raise FileNotFoundError(f"completed HFSS spectrum missing: {path}")
        data = read_spectrum(path)
        candidates = absorption_resonances(data[:, 0], data[:, 3])
        plot_spectrum(data, candidates, RESONANCES_ROOT / f"{request.case_name}_{request.angle_deg}deg_TRA", f"{request.case_name}, {request.angle_deg} degrees")
        rows.extend(
            (request.case_name, request.angle_deg, candidate.identifier, candidate.wavelength_nm, candidate.absorption, candidate.prominence)
            for candidate in candidates
        )
    with (RESONANCES_ROOT / "resonance_candidates.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("case", "angle_deg", "candidate_id", "wavelength_nm", "absorption", "prominence"))
        writer.writerows(rows)


if __name__ == "__main__":
    review_all()
