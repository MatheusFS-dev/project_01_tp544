"""Palik material-table validation and HFSS-ready loss conversion."""

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np


EPSILON_0_F_PER_M = 8.854_187_8128e-12
LOWER_FREQUENCY_HZ = 299_792_458.0 / (1600e-9)
UPPER_FREQUENCY_HZ = 299_792_458.0 / (450e-9)
EXPECTED_HEADER = ("frequency_hz", "epsilon_real", "epsilon_imag")


@dataclass(frozen=True)
class PalikData:
    """Validated complex relative-permittivity samples from Lumerical.

    Attributes:
        path: Exported Lumerical sampled-data table.
        source_hash: SHA-256 digest of the exact table used in HFSS.
        frequency_hz: Strictly ascending sample frequencies.
        epsilon_real: Real relative permittivity.
        epsilon_imag: Positive loss convention from the Lumerical table.
    """

    path: Path
    source_hash: str
    frequency_hz: np.ndarray
    epsilon_real: np.ndarray
    epsilon_imag: np.ndarray

    @property
    def conductivity_s_per_m(self) -> np.ndarray:
        """Convert positive imaginary permittivity to positive HFSS conductivity.

        Returns:
            Frequency-dependent bulk conductivity in S/m. HFSS represents
            material loss through positive conductivity under its phasor
            convention.
        """
        return 2 * np.pi * self.frequency_hz * EPSILON_0_F_PER_M * self.epsilon_imag


def read_palik_table(path: Path) -> PalikData:
    """Read and validate an exported ``Au (Gold) - Palik`` sampled-data table.

    Args:
        path: Tab-separated file with frequency and complex-permittivity parts.

    Returns:
        Validated material data and a digest of the source bytes.

    Raises:
        ValueError: If the header, values, passivity, ordering, or spectral
            coverage does not satisfy the HFSS campaign requirements.
    """
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"Palik table not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    if not rows or tuple(rows[0]) != EXPECTED_HEADER:
        raise ValueError(f"expected TSV header {EXPECTED_HEADER!r}")
    if len(rows) < 2:
        raise ValueError("Palik table contains no samples")

    try:
        values = np.asarray([[float(value) for value in row] for row in rows[1:]], dtype=float)
    except ValueError as error:
        raise ValueError("Palik table contains a nonnumeric value") from error
    if values.ndim != 2 or values.shape[1] != 3 or not np.isfinite(values).all():
        raise ValueError("Palik table must contain three finite numeric columns")

    frequency_hz = values[:, 0]
    epsilon_real = values[:, 1]
    epsilon_imag = values[:, 2]
    if np.any(frequency_hz <= 0) or np.any(np.diff(frequency_hz) <= 0):
        raise ValueError("Palik frequencies must be positive and strictly ascending")
    if np.any(epsilon_imag <= 0):
        raise ValueError("Palik samples must use strictly positive loss")
    if frequency_hz[0] > LOWER_FREQUENCY_HZ or frequency_hz[-1] < UPPER_FREQUENCY_HZ:
        raise ValueError("Palik table does not cover 450 to 1600 nm")

    return PalikData(
        path=path,
        source_hash=hashlib.sha256(path.read_bytes()).hexdigest(),
        frequency_hz=frequency_hz,
        epsilon_real=epsilon_real,
        epsilon_imag=epsilon_imag,
    )
