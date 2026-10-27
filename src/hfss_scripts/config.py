"""Immutable physical and numerical configuration for the HFSS campaign."""

from dataclasses import dataclass


SPEED_OF_LIGHT_M_PER_S = 299_792_458.0


@dataclass(frozen=True)
class ModelConfig:
    """Store the fixed assignment configuration in SI units.

    Attributes:
        outer_radius_m: Outer Au ring radius.
        inner_radius_m: Inner Au ring radius.
        thickness_m: Au ring thickness.
        period_m: Square lattice period.
        port_offset_m: Distance from the ring plane to each Floquet port.
        silica_index: Non-dispersive SiO2 refractive index used in case B.
        angles_deg: Incidence angles from the positive z normal.
        wavelength_nm: Inclusive wavelength range for report products.
        case_names: Stable names for the suspended and substrate cases.
    """

    outer_radius_m: float = 90e-9
    inner_radius_m: float = 50e-9
    thickness_m: float = 30e-9
    period_m: float = 320e-9
    port_offset_m: float = 800e-9
    silica_index: float = 1.45
    angles_deg: tuple[int, ...] = (0, 15, 30, 45, 60, 75)
    wavelength_nm: tuple[int, int] = (450, 1600)
    case_names: tuple[str, str] = ("A_air", "B_SiO2")

    @property
    def frequency_thz(self) -> tuple[float, float]:
        """Return the ascending solution-frequency limits in THz.

        Returns:
            Lower and upper frequencies corresponding to 1600 and 450 nm.
        """
        lower_wavelength_m = self.wavelength_nm[1] * 1e-9
        upper_wavelength_m = self.wavelength_nm[0] * 1e-9
        return (
            SPEED_OF_LIGHT_M_PER_S / lower_wavelength_m / 1e12,
            SPEED_OF_LIGHT_M_PER_S / upper_wavelength_m / 1e12,
        )

    @property
    def wavelength_grid_nm(self) -> tuple[int, ...]:
        """Return the required inclusive 1 nm reporting grid.

        Returns:
            Wavelengths from 450 through 1600 nm.
        """
        return tuple(range(self.wavelength_nm[0], self.wavelength_nm[1] + 1))
