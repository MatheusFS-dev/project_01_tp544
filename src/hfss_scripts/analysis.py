"""Independent post-processing helpers for HFSS spectra."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FresnelPower:
    """Power coefficients of a lossless TM Fresnel interface."""

    reflection: float
    transmission: float


@dataclass(frozen=True)
class ResampledSpectrum:
    """Spectrum on a requested wavelength grid with provenance flags."""

    wavelength_nm: np.ndarray
    transmission: np.ndarray
    reflection: np.ndarray
    interpolated: np.ndarray


@dataclass(frozen=True)
class Traspectrum:
    """Native-frequency transmission, reflection, and absorption results."""

    frequency_hz: np.ndarray
    wavelength_nm: np.ndarray
    transmission: np.ndarray
    reflection: np.ndarray
    absorption: np.ndarray


def fresnel_tm(index_incident: float, index_transmitted: float, angle_deg: float) -> FresnelPower:
    """Calculate lossless Fresnel TM power coefficients.

    Args:
        index_incident: Refractive index in the incident half space.
        index_transmitted: Refractive index in the transmitting half space.
        angle_deg: Incidence angle from the interface normal in degrees.

    Returns:
        TM reflected and transmitted powers.

    Raises:
        ValueError: If an index is nonpositive or total internal reflection is
            requested for this lossless reference.
    """
    if index_incident <= 0 or index_transmitted <= 0:
        raise ValueError("refractive indices must be positive")
    theta_incident = np.deg2rad(angle_deg)
    sine_transmitted = index_incident * np.sin(theta_incident) / index_transmitted
    if abs(sine_transmitted) > 1:
        raise ValueError("total internal reflection is outside this reference model")
    cosine_incident = np.cos(theta_incident)
    cosine_transmitted = np.sqrt(1 - sine_transmitted**2)
    denominator = index_transmitted * cosine_incident + index_incident * cosine_transmitted
    reflection_amplitude = (
        index_transmitted * cosine_incident - index_incident * cosine_transmitted
    ) / denominator
    transmission_amplitude = 2 * index_incident * cosine_incident / denominator
    reflection = float(abs(reflection_amplitude) ** 2)
    transmission = float(
        index_transmitted * cosine_transmitted / (index_incident * cosine_incident)
        * abs(transmission_amplitude) ** 2
    )
    return FresnelPower(reflection=reflection, transmission=transmission)


def specular_tra_spectrum(
    frequency_hz: np.ndarray,
    reflection_s: np.ndarray,
    transmission_s: np.ndarray,
) -> Traspectrum:
    """Calculate T, R, and A from power-normalized specular S parameters.

    This helper applies only where the specular modes are the sole propagating
    Floquet channels. It intentionally retains negative absorption so that a
    nonphysical energy balance remains visible to validation.

    Args:
        frequency_hz: Strictly increasing native HFSS frequencies in hertz.
        reflection_s: Complex S parameters from the incident mode to the
            specular mode of the incident-side port.
        transmission_s: Complex S parameters from the incident mode to the
            specular mode of the opposite-side port.

    Returns:
        Native T/R/A spectrum, without clipping absorption.

    Raises:
        ValueError: If arrays are incompatible or frequencies are invalid.
    """
    frequency_hz = np.asarray(frequency_hz, dtype=float)
    reflection_s = np.asarray(reflection_s, dtype=complex)
    transmission_s = np.asarray(transmission_s, dtype=complex)
    if frequency_hz.ndim != 1 or len(frequency_hz) < 2:
        raise ValueError("frequency axis must contain at least two samples")
    if np.any(~np.isfinite(frequency_hz)) or np.any(frequency_hz <= 0) or np.any(np.diff(frequency_hz) <= 0):
        raise ValueError("frequencies must be finite, positive, and strictly increasing")
    if reflection_s.shape != frequency_hz.shape or transmission_s.shape != frequency_hz.shape:
        raise ValueError("S-parameter arrays must match the frequency axis")
    reflection = np.abs(reflection_s) ** 2
    transmission = np.abs(transmission_s) ** 2
    return Traspectrum(
        frequency_hz=frequency_hz,
        wavelength_nm=299_792_458 / frequency_hz * 1e9,
        transmission=transmission,
        reflection=reflection,
        absorption=1 - reflection - transmission,
    )


def propagating_modal_tra_spectrum(
    frequency_hz: np.ndarray,
    reflected_s: np.ndarray,
    transmitted_s: np.ndarray,
    reflected_gamma: np.ndarray,
    transmitted_gamma: np.ndarray,
) -> Traspectrum:
    """Sum power in every propagating Floquet mode reported by HFSS.

    HFSS reports the port propagation constant as ``Gamma``. In lossless port
    media a mode is propagating when its attenuation, ``real(Gamma)``, is zero.
    Evanescent modes remain in the S-matrix for port accuracy but carry no
    power to the semi-infinite exterior and are excluded from T/R.

    Args:
        frequency_hz: Strictly increasing native HFSS frequencies in hertz.
        reflected_s: Complex S matrix slice with dimensions ``modes, samples``
            for the incident-side Floquet port.
        transmitted_s: Complex S matrix slice with dimensions ``modes, samples``
            for the opposite-side Floquet port.
        reflected_gamma: Complex modal propagation constants for reflected_s.
        transmitted_gamma: Complex modal propagation constants for
            transmitted_s.

    Returns:
        Native total T/R/A spectrum without clipping absorption.

    Raises:
        ValueError: If modal axes do not agree, a port lacks a propagating
            channel, or the highest included mode propagates. The latter means
            the Floquet modal truncation must be increased before reporting.
    """
    frequency_hz = np.asarray(frequency_hz, dtype=float)
    reflected_s = np.asarray(reflected_s, dtype=complex)
    transmitted_s = np.asarray(transmitted_s, dtype=complex)
    reflected_gamma = np.asarray(reflected_gamma, dtype=complex)
    transmitted_gamma = np.asarray(transmitted_gamma, dtype=complex)
    expected_shape = (None, len(frequency_hz))
    if frequency_hz.ndim != 1 or len(frequency_hz) < 2:
        raise ValueError("frequency axis must contain at least two samples")
    if np.any(~np.isfinite(frequency_hz)) or np.any(frequency_hz <= 0) or np.any(np.diff(frequency_hz) <= 0):
        raise ValueError("frequencies must be finite, positive, and strictly increasing")
    modal_arrays = (reflected_s, transmitted_s, reflected_gamma, transmitted_gamma)
    if any(values.ndim != 2 or values.shape[1] != expected_shape[1] for values in modal_arrays):
        raise ValueError("modal data must have dimensions modes by frequency samples")
    if reflected_s.shape != reflected_gamma.shape or transmitted_s.shape != transmitted_gamma.shape:
        raise ValueError("each S matrix must match its Gamma matrix")
    reflected_propagating = np.isclose(reflected_gamma.real, 0.0, rtol=0, atol=1e-6)
    transmitted_propagating = np.isclose(transmitted_gamma.real, 0.0, rtol=0, atol=1e-6)
    if not np.all(np.any(reflected_propagating, axis=0)) or not np.all(np.any(transmitted_propagating, axis=0)):
        raise ValueError("each port must retain at least one propagating Floquet mode")
    if np.any(reflected_propagating[-1]) or np.any(transmitted_propagating[-1]):
        raise ValueError("highest included Floquet mode propagates, so modal convergence is inconclusive")
    reflection = np.sum(np.abs(reflected_s) ** 2 * reflected_propagating, axis=0)
    transmission = np.sum(np.abs(transmitted_s) ** 2 * transmitted_propagating, axis=0)
    return Traspectrum(
        frequency_hz=frequency_hz,
        wavelength_nm=299_792_458 / frequency_hz * 1e9,
        transmission=transmission,
        reflection=reflection,
        absorption=1 - reflection - transmission,
    )


def resample_spectrum(
    wavelength_nm: np.ndarray,
    transmission: np.ndarray,
    reflection: np.ndarray,
    requested_nm: np.ndarray,
) -> ResampledSpectrum:
    """Linearly resample native HFSS data and label non-native output points.

    Args:
        wavelength_nm: Strictly increasing native wavelength samples.
        transmission: Native power transmission samples.
        reflection: Native power reflection samples.
        requested_nm: Wavelength grid required by the report.

    Returns:
        Resampled powers plus a Boolean interpolation provenance mask.

    Raises:
        ValueError: If dimensions, ordering, or requested coverage are invalid.
    """
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    transmission = np.asarray(transmission, dtype=float)
    reflection = np.asarray(reflection, dtype=float)
    requested_nm = np.asarray(requested_nm, dtype=float)
    if wavelength_nm.ndim != 1 or requested_nm.ndim != 1:
        raise ValueError("wavelength axes must be one-dimensional")
    if len(wavelength_nm) < 2 or len(transmission) != len(wavelength_nm) or len(reflection) != len(wavelength_nm):
        raise ValueError("native spectrum arrays must share at least two samples")
    if np.any(np.diff(wavelength_nm) <= 0):
        raise ValueError("native wavelengths must be strictly increasing")
    endpoint_tolerance_nm = 1e-6
    if (
        requested_nm.min() < wavelength_nm[0] - endpoint_tolerance_nm
        or requested_nm.max() > wavelength_nm[-1] + endpoint_tolerance_nm
    ):
        raise ValueError("requested wavelengths exceed the native HFSS range")
    native_matches = np.isclose(
        requested_nm[:, np.newaxis],
        wavelength_nm[np.newaxis, :],
        rtol=0,
        atol=endpoint_tolerance_nm,
    )
    interpolated = ~np.any(native_matches, axis=1)
    return ResampledSpectrum(
        wavelength_nm=requested_nm,
        transmission=np.interp(requested_nm, wavelength_nm, transmission),
        reflection=np.interp(requested_nm, wavelength_nm, reflection),
        interpolated=interpolated,
    )
