"""Review completed 5 nm fast spectra and prepare fast near-field requests."""

import csv
import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SPECTRA = ROOT / "results" / "spectra"
RESONANCES = ROOT / "results" / "resonances_fast"
CASES = {"A_air": 1, "B_SiO2": 2}
ANGLES = (0, 15, 30, 45, 60, 75)
FIELDS = ("case", "angle_deg", "resonance_id", "wavelength_nm")
FAST_WAVELENGTHS_NM = np.arange(450, 1601, 5)

base_spec = importlib.util.spec_from_file_location(
    "review_resonances", Path(__file__).with_name("03_review_resonances.py"))
base = importlib.util.module_from_spec(base_spec)
base_spec.loader.exec_module(base)


def read_spectrum(path):
    """Read one complete exploratory 5 nm T/R/A spectrum."""
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["wavelength_nm", "T", "R", "A"]:
            raise ValueError(f"Invalid spectrum columns: {path}")
        data = np.array([[float(row[key]) for key in reader.fieldnames] for row in reader])
    if data.shape != (FAST_WAVELENGTHS_NM.size, 4) or not np.isfinite(data).all():
        raise ValueError(f"Incomplete or nonfinite fast spectrum: {path}")
    data = data[np.argsort(data[:, 0])]
    if not np.allclose(data[:, 0], FAST_WAVELENGTHS_NM, atol=1e-5):
        raise ValueError(f"Invalid 5 nm wavelength grid: {path}")
    if not np.allclose(data[:, 1:].sum(axis=1), 1, atol=1e-8):
        raise ValueError(f"Power balance failure: {path}")
    return data


def completed_spectra(directory):
    """Return only complete fast spectra, identified by their fast_ prefix."""
    found = {}
    for case in CASES:
        for angle in ANGLES:
            stem = f"fast_{case}_{angle}deg"
            if (directory / f"{stem}.done").exists():
                found[(case, angle)] = read_spectrum(directory / f"{stem}.tsv")
    return found


def read_resonances(path, spectra):
    """Validate manually selected resonances on the 5 nm fast grid."""
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
            wavelength_int = int(round(wavelength))
            if case not in CASES or angle not in ANGLES or (case, angle) not in spectra:
                raise ValueError(f"Unknown or unfinished spectrum: {row}")
            if (not identifier or not np.isfinite(wavelength)
                    or wavelength_int not in FAST_WAVELENGTHS_NM
                    or abs(wavelength - wavelength_int) > 1e-6):
                raise ValueError(f"Resonance is off the 5 nm fast grid: {row}")
            key = (case, angle, identifier)
            position = (case, angle, wavelength_int)
            if key in keys or position in positions:
                raise ValueError(f"Duplicate resonance: {row}")
            keys.add(key)
            positions.add(position)
            rows.append((case, angle, identifier, wavelength_int))
    return sorted(rows, key=lambda row: (row[2], row[1], row[0]))


def main():
    """Render completed fast spectra and write isolated fast requests."""
    base.style()
    spectra = completed_spectra(SPECTRA)
    selections = read_resonances(RESONANCES / "resonances.csv", spectra)
    base.write_outputs(spectra, selections, RESONANCES)


if __name__ == "__main__":
    main()
