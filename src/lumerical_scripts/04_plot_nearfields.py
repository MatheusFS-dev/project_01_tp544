"""Plot selected electric intensity maps from completed near-field runs."""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

import importlib.util

review_spec = importlib.util.spec_from_file_location("review_resonances", Path(__file__).with_name("03_review_resonances.py"))
review = importlib.util.module_from_spec(review_spec)
review_spec.loader.exec_module(review)
ROOT = Path(__file__).resolve().parents[2]


def _load_mat_arrays(path):
    """Read legacy MAT files and MATLAB v7.3/HDF5 files from Lumerical."""
    try:
        return loadmat(path)
    except (NotImplementedError, ValueError):
        import h5py

        if not h5py.is_hdf5(path):
            raise
        with h5py.File(path, "r") as stream:
            # MATLAB v7.3 stores arrays in column-major order. h5py exposes
            # their dimensions reversed, so transpose them back to MATLAB order.
            return {key: np.asarray(stream[key]).T for key in stream
                    if isinstance(stream[key], h5py.Dataset)}


def read_fields(path):
    """Load and validate one Lumerical MAT field export."""
    data = _load_mat_arrays(path)
    arrays = {}
    for key in ("x_nm", "y_nm", "z_nm", "E2_xy", "E2_xz"):
        arrays[key] = np.asarray(data[key], dtype=float).squeeze()
    x, y, z = (arrays[key] for key in ("x_nm", "y_nm", "z_nm"))
    if arrays["E2_xy"].shape != (x.size, y.size) or arrays["E2_xz"].shape != (x.size, z.size):
        raise ValueError(f"Field dimensions do not match coordinates: {path}")
    if not all(np.isfinite(value).all() for value in arrays.values()):
        raise ValueError(f"Nonfinite field data: {path}")
    arrays["requested_wavelength_nm"] = float(np.asarray(data["requested_wavelength_nm"]).squeeze())
    arrays["sampled_wavelength_nm"] = float(np.asarray(data["sampled_wavelength_nm"]).squeeze())
    return arrays


def plot_maps(data, stem, title, limits):
    """Save XY and XZ intensity maps using supplied shared color limits."""
    for plane, horizontal, vertical in (("xy", "x_nm", "y_nm"), ("xz", "x_nm", "z_nm")):
        fig, ax = plt.subplots(figsize=(5.3, 4.2))
        image = ax.pcolormesh(data[horizontal], data[vertical], data[f"E2_{plane}"].T,
                              cmap="cividis", shading="auto", vmin=0, vmax=limits[plane])
        ax.set(xlabel="x (nm)", ylabel=f"{plane[1]} (nm)",
               title=f"{title}, {plane.upper()}, {data['sampled_wavelength_nm']:.1f} nm")
        fig.colorbar(image, ax=ax, label="|E|²/|E₀|²")
        review.save_figure(fig, stem.with_name(stem.name + f"_{plane}"))


def main():
    """Render complete requests, sharing scale between matched modes."""
    review.style()
    resonances = ROOT / "results" / "resonances"
    nearfields = ROOT / "results" / "nearfields"
    with (resonances / "nearfield_manifest.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    loaded = {}
    for row in rows:
        request_id = int(row["request_id"])
        stem = nearfields / f"request_{request_id}"
        if stem.with_suffix(".done").exists():
            data = read_fields(stem.with_suffix(".mat"))
            if abs(data["requested_wavelength_nm"] - float(row["wavelength_nm"])) > 1e-6:
                raise ValueError(f"Request wavelength mismatch: {stem}")
            loaded[request_id] = data
    for row in rows:
        request_id = int(row["request_id"])
        if request_id not in loaded:
            continue
        peers = [loaded[int(peer["request_id"])] for peer in rows
                 if int(peer["request_id"]) in loaded and peer["angle_deg"] == row["angle_deg"]
                 and peer["resonance_id"] == row["resonance_id"]]
        limits = {plane: max(float(np.max(data[f"E2_{plane}"])) for data in peers) for plane in ("xy", "xz")}
        title = f"{row['case']}, {row['angle_deg']}°, {row['resonance_id']}"
        plot_maps(loaded[request_id], nearfields / f"request_{request_id}", title, limits)


if __name__ == "__main__":
    main()
