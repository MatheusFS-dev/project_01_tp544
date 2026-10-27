"""Retrieve and convert the user-authorized external Palik Au table."""

import csv
import hashlib
import json
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from urllib.request import urlopen

import numpy as np


PALIK_SOURCE_URL = "https://raw.githubusercontent.com/ulfgri/nk/master/windt/Au_palik.nk"
SPEED_OF_LIGHT_M_PER_S = 299_792_458.0


@dataclass(frozen=True)
class ExternalPalikSamples:
    """Complex Au samples converted directly from a Palik n-k table."""

    frequency_hz: np.ndarray
    epsilon_real: np.ndarray
    epsilon_imag: np.ndarray


def parse_palik_nk(raw_text: str) -> ExternalPalikSamples:
    """Convert Palik ``lambda(A), n, k`` rows into increasing-frequency εr.

    Args:
        raw_text: Original text from the approved external table.

    Returns:
        Passivity-preserving relative-permittivity samples sorted by frequency.

    Raises:
        ValueError: If no valid, positive n-k samples are found.
    """
    rows = []
    for line in raw_text.splitlines():
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        fields = line.split()
        if len(fields) != 3:
            continue
        try:
            wavelength_angstrom, refractive_index, extinction = map(float, fields)
        except ValueError:
            continue
        if wavelength_angstrom > 0 and refractive_index > 0 and extinction > 0:
            rows.append((wavelength_angstrom, refractive_index, extinction))
    if len(rows) < 2:
        raise ValueError("external Palik source contains fewer than two valid n-k samples")
    samples = np.asarray(rows, dtype=float)
    frequency_hz = SPEED_OF_LIGHT_M_PER_S / (samples[:, 0] * 1e-10)
    order = np.argsort(frequency_hz, kind="stable")
    frequency_hz = frequency_hz[order]
    refractive_index = samples[:, 1][order]
    extinction = samples[:, 2][order]
    keep = np.concatenate(([True], np.diff(frequency_hz) > 0))
    frequency_hz = frequency_hz[keep]
    refractive_index = refractive_index[keep]
    extinction = extinction[keep]
    return ExternalPalikSamples(
        frequency_hz=frequency_hz,
        epsilon_real=refractive_index**2 - extinction**2,
        epsilon_imag=2 * refractive_index * extinction,
    )


def download_external_palik(url: str = PALIK_SOURCE_URL) -> tuple[str, str]:
    """Download the authorized external source and calculate its SHA-256 hash.

    Args:
        url: Immutable or raw URL of the Palik n-k table.

    Returns:
        Decoded source text and SHA-256 digest of its UTF-8 bytes.
    """
    with urlopen(url, timeout=30) as response:
        raw_bytes = response.read()
    raw_text = raw_bytes.decode("utf-8")
    return raw_text, hashlib.sha256(raw_bytes).hexdigest()


def write_hfss_table(output_path: Path, source_url: str = PALIK_SOURCE_URL) -> Path:
    """Download, convert, and record external Palik data for the HFSS importer.

    Args:
        output_path: TSV target read by :func:`material.read_palik_table`.
        source_url: URL used to retrieve the Palik n-k data.

    Returns:
        Path to a sidecar JSON provenance record.
    """
    raw_text, source_hash = download_external_palik(source_url)
    samples = parse_palik_nk(raw_text)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("frequency_hz", "epsilon_real", "epsilon_imag"))
        writer.writerows(zip(samples.frequency_hz, samples.epsilon_real, samples.epsilon_imag))
    metadata_path = output_path.with_suffix(".provenance.json")
    metadata_path.write_text(
        json.dumps(
            {
                "source_url": source_url,
                "source_sha256": source_hash,
                "retrieved_utc": datetime.now(UTC).isoformat(),
                "conversion": "epsilon_r=(n+i*k)^2 from lambda_angstrom,n,k rows",
                "samples": int(len(samples.frequency_hz)),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return metadata_path
