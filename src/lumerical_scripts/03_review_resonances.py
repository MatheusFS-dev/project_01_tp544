"""Review completed spectra and prepare manually selected near-field runs."""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SPECTRA = ROOT / "results" / "spectra"
RESONANCES = ROOT / "results" / "resonances"
CASES = {"A_air": 1, "B_SiO2": 2}
ANGLES = (0, 15, 30, 45, 60, 75)
FIELDS = ("case", "angle_deg", "resonance_id", "wavelength_nm")


def style():
    """Apply a restrained, consistent style to all scientific plots."""
    plt.rcParams.update({"font.family": "DejaVu Serif", "figure.facecolor": "white",
                         "axes.facecolor": "white", "savefig.facecolor": "white",
                         "axes.spines.top": False, "axes.spines.right": False})


def read_spectrum(path):
    """Read and validate a complete 1 nm T/R/A spectral export.

    Args:
        path: Path to a completed spectral TSV.

    Returns:
        Array with columns wavelength, T, R, and A in ascending order.

    Raises:
        ValueError: If the grid, columns, or power balance is invalid.
    """
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["wavelength_nm", "T", "R", "A"]:
            raise ValueError(f"Invalid spectrum columns: {path}")
        data = np.array([[float(row[key]) for key in reader.fieldnames] for row in reader])
    if data.shape != (1151, 4) or not np.isfinite(data).all():
        raise ValueError(f"Incomplete or nonfinite spectrum: {path}")
    data = data[np.argsort(data[:, 0])]
    if not np.allclose(data[:, 0], np.arange(450, 1601), atol=1e-5):
        raise ValueError(f"Invalid wavelength grid: {path}")
    if not np.allclose(data[:, 1:].sum(axis=1), 1, atol=1e-8):
        raise ValueError(f"Power balance failure: {path}")
    return data


def completed_spectra(directory):
    """Return only spectra with a completion marker and valid case and angle."""
    found = {}
    for case in CASES:
        for angle in ANGLES:
            stem = f"{case}_{angle}deg"
            if (directory / f"{stem}.done").exists():
                path = directory / f"{stem}.tsv"
                found[(case, angle)] = read_spectrum(path)
    return found


def read_resonances(path, spectra):
    """Validate manual resonance selections against completed spectra."""
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x", newline="", encoding="utf-8") as stream:
            csv.writer(stream).writerow(FIELDS)
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != list(FIELDS):
            raise ValueError("resonances.csv has invalid columns")
        rows = []
        keys = set()
        positions = set()
        for row in reader:
            case = row["case"].strip()
            angle = int(row["angle_deg"])
            identifier = row["resonance_id"].strip()
            wavelength = float(row["wavelength_nm"])
            if case not in CASES or angle not in ANGLES or (case, angle) not in spectra:
                raise ValueError(f"Unknown or unfinished spectrum: {row}")
            if not identifier or not np.isfinite(wavelength) or not 450 <= wavelength <= 1600:
                raise ValueError(f"Invalid resonance: {row}")
            if abs(wavelength - round(wavelength)) > 1e-6:
                raise ValueError(f"Resonance is off the 1 nm grid: {row}")
            key = (case, angle, identifier)
            if key in keys:
                raise ValueError(f"Duplicate resonance: {key}")
            position = (case, angle, int(round(wavelength)))
            if position in positions:
                raise ValueError(f"Duplicate resonance wavelength: {position}")
            keys.add(key)
            positions.add(position)
            rows.append((case, angle, identifier, int(round(wavelength))))
    return sorted(rows, key=lambda row: (row[2], row[1], row[0]))


def save_figure(fig, stem):
    """Save matching PNG and PDF figures."""
    for suffix in ("png", "pdf"):
        fig.savefig(stem.with_suffix(f".{suffix}"), dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_spectrum(data, selections, stem, title):
    """Plot T, R, and A with all selected resonances marked."""
    fig, ax = plt.subplots(figsize=(7, 4))
    for index, (label, color) in enumerate((("T", "#405e78"), ("R", "#a16c54"), ("A", "#667c61")), 1):
        ax.plot(data[:, 0], data[:, index], label=label, color=color, lw=1.25)
    for _, _, identifier, wavelength in selections:
        point_index = int(np.argmin(np.abs(data[:, 0] - wavelength)))
        reflection = data[point_index, 2]
        ax.axvline(wavelength, color="#777777", lw=.7, ls="--")
        ax.scatter(wavelength, reflection, s=34, marker="o", color="#222222",
                   edgecolors="white", linewidths=.6, zorder=4,
                   label=f"{identifier}: {wavelength} nm")
    ax.set(xlabel="Comprimento de onda (nm)", ylabel="Fração de potência", title=title, xlim=(450, 1600))
    ax.legend(frameon=False, ncol=3)
    save_figure(fig, stem)


def write_outputs(spectra, selections, directory):
    """Write tracking table, angle plot, and numeric near-field requests."""
    directory.mkdir(parents=True, exist_ok=True)
    for (case, angle), data in spectra.items():
        chosen = [row for row in selections if row[:2] == (case, angle)]
        plot_spectrum(data, chosen, directory / f"{case}_{angle}deg_TRA", f"{case}, {angle}°")
    grouped = sorted({(angle, identifier) for _, angle, identifier, _ in selections})
    lookup = {(case, angle, identifier): wavelength for case, angle, identifier, wavelength in selections}
    with (directory / "resonance_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("angle_deg", "resonance_id", "lambda_air_nm", "lambda_SiO2_nm"))
        for angle, identifier in grouped:
            writer.writerow((angle, identifier, lookup.get(("A_air", angle, identifier), ""),
                             lookup.get(("B_SiO2", angle, identifier), "")))
    requests = [row for row in selections if row[1] in (0, 45)]
    with (directory / "nearfield_requests.tsv").open("w", newline="", encoding="ascii") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("request_id", "case_index", "angle_deg", "wavelength_nm"))
        for request_id, (case, angle, _, wavelength) in enumerate(requests, 1):
            writer.writerow((request_id, CASES[case], angle, wavelength))
    with (directory / "nearfield_manifest.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("request_id",) + FIELDS)
        for request_id, row in enumerate(requests, 1):
            writer.writerow((request_id,) + row)
    fig, ax = plt.subplots(figsize=(7, 4))
    for identifier in sorted({row[2] for row in selections}):
        for case, color, marker in (("A_air", "#405e78", "o"), ("B_SiO2", "#a16c54", "s")):
            points = sorted((angle, wavelength) for c, angle, i, wavelength in selections if c == case and i == identifier)
            if points:
                ax.plot(*zip(*points), marker=marker, color=color, lw=1, label=f"{identifier}, {case}")
    ax.set(xlabel="Ângulo de incidência (°)", ylabel="Comprimento de onda ressonante (nm)", xlim=(-2, 77))
    if selections:
        ax.legend(frameon=False, fontsize=8)
    save_figure(fig, directory / "resonance_vs_angle")


def main():
    """Process all completed results without modifying manual selections."""
    style()
    spectra = completed_spectra(SPECTRA)
    selections = read_resonances(RESONANCES / "resonances.csv", spectra)
    write_outputs(spectra, selections, RESONANCES)


if __name__ == "__main__":
    main()
