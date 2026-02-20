"""
Unit tests for Compton scattering physics functions.
"""
import sys
import os
import numpy as np
import pytest

# Allow importing from parent directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from compton import (
    LAMBDA_C, EV, H_PLANCK, C,
    delta_lambda, scattered_lambda, photon_energy_eV, recoil_angle,
)


# ── delta_lambda ─────────────────────────────────────────

class TestDeltaLambda:
    def test_zero_angle(self):
        """No scattering → zero shift."""
        assert delta_lambda(0) == pytest.approx(0.0, abs=1e-20)

    def test_90_degrees(self):
        """θ = 90° → Δλ = λ_C ≈ 2.426 pm."""
        result = delta_lambda(90)
        assert result == pytest.approx(LAMBDA_C, rel=1e-6)

    def test_180_degrees(self):
        """θ = 180° → Δλ = 2λ_C ≈ 4.852 pm (maximum)."""
        result = delta_lambda(180)
        assert result == pytest.approx(2 * LAMBDA_C, rel=1e-6)

    def test_known_value_90deg_pm(self):
        """90° shift in picometres matches textbook value."""
        dlam_pm = delta_lambda(90) * 1e12
        assert dlam_pm == pytest.approx(2.426, abs=0.001)

    def test_known_value_180deg_pm(self):
        """180° shift in picometres ≈ 4.852 pm."""
        dlam_pm = delta_lambda(180) * 1e12
        assert dlam_pm == pytest.approx(4.852, abs=0.002)

    def test_monotonically_increasing(self):
        """Δλ is non-decreasing from 0° to 180°."""
        angles = np.linspace(0, 180, 181)
        values = np.array([delta_lambda(a) for a in angles])
        diffs = np.diff(values)
        assert np.all(diffs >= -1e-25), "delta_lambda must be non-decreasing"

    def test_symmetry_around_90(self):
        """Δλ(180°−θ) != Δλ(θ) in general — just verify 45° < 90°."""
        assert delta_lambda(45) < delta_lambda(90)
        assert delta_lambda(90) < delta_lambda(135)

    @pytest.mark.parametrize("theta", [0, 45, 90, 135, 180])
    def test_formula_manual(self, theta):
        """Cross-check against formula directly."""
        expected = LAMBDA_C * (1 - np.cos(np.deg2rad(theta)))
        assert delta_lambda(theta) == pytest.approx(expected, rel=1e-9)


# ── scattered_lambda ──────────────────────────────────────

class TestScatteredLambda:
    @pytest.mark.parametrize("theta", np.linspace(0, 180, 19))
    def test_always_geq_lambda0(self, theta):
        """Scattered wavelength is always ≥ incident wavelength."""
        lam0 = 70.8e-12   # Mo Kα
        assert scattered_lambda(lam0, theta) >= lam0 - 1e-25

    def test_zero_angle_returns_lambda0(self):
        lam0 = 50e-12
        assert scattered_lambda(lam0, 0) == pytest.approx(lam0, rel=1e-9)

    def test_90_degrees_shift(self):
        lam0 = 70.8e-12
        lamp = scattered_lambda(lam0, 90)
        assert lamp == pytest.approx(lam0 + LAMBDA_C, rel=1e-6)

    def test_180_degrees_shift(self):
        lam0 = 100e-12
        lamp = scattered_lambda(lam0, 180)
        assert lamp == pytest.approx(lam0 + 2 * LAMBDA_C, rel=1e-6)


# ── photon_energy_eV ──────────────────────────────────────

class TestPhotonEnergyEV:
    def test_visible_red(self):
        """700 nm red light → ~1.77 eV."""
        E = photon_energy_eV(700e-9)
        assert E == pytest.approx(1.771, abs=0.002)

    def test_visible_violet(self):
        """400 nm violet → ~3.10 eV."""
        E = photon_energy_eV(400e-9)
        assert E == pytest.approx(3.101, abs=0.002)

    def test_xray_mo_kalpha(self):
        """Mo Kα at 70.8 pm → ~17.5 keV."""
        E = photon_energy_eV(70.8e-12)
        assert E == pytest.approx(17530, rel=0.005)

    def test_energy_inversely_proportional_to_wavelength(self):
        """Doubling wavelength halves energy."""
        lam1 = 50e-12
        lam2 = 100e-12
        E1 = photon_energy_eV(lam1)
        E2 = photon_energy_eV(lam2)
        assert E1 == pytest.approx(2 * E2, rel=1e-6)

    def test_formula_consistency(self):
        """E = hc/λ in eV matches manual calculation."""
        lam = 0.5e-9   # 0.5 nm
        expected = H_PLANCK * C / (lam * EV)
        assert photon_energy_eV(lam) == pytest.approx(expected, rel=1e-9)


# ── recoil_angle ─────────────────────────────────────────

class TestRecoilAngle:
    def test_zero_scattering_gives_zero_recoil(self):
        assert recoil_angle(70.8e-12, 0) == pytest.approx(0.0, abs=1e-10)

    def test_180_scattering_gives_zero_recoil(self):
        """Backscatter → electron goes forward (φ → 0)."""
        assert recoil_angle(70.8e-12, 180) == pytest.approx(0.0, abs=1e-10)

    def test_90_degree_scattering(self):
        """At 90°, recoil angle is less than 45°."""
        phi = recoil_angle(70.8e-12, 90)
        assert 0 < phi < 45

    def test_recoil_always_forward_hemisphere(self):
        """Electron recoil angle φ is always in [0°, 90°)."""
        lam0 = 70.8e-12
        for theta in range(1, 180):
            phi = recoil_angle(lam0, theta)
            assert 0 <= phi < 90, f"φ={phi:.2f}° out of range at θ={theta}°"
