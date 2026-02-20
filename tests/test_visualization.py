"""
Smoke tests for Compton visualization: figure construction and sliders.

We use a non-interactive Matplotlib backend so the tests work headlessly.
"""
import sys
import os
import pytest

# Force non-interactive backend BEFORE importing matplotlib.pyplot
import matplotlib
matplotlib.use('Agg')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


@pytest.fixture(scope='module')
def fig_and_widgets():
    """Build the figure once for all tests in this module."""
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider
    import numpy as np
    from compton import (
        LAMBDA_C, delta_lambda, scattered_lambda,
        photon_energy_eV, recoil_angle,
    )

    # Reproduce the figure construction from main() without plt.show()
    import importlib, types

    # We run the module-level main() but monkey-patch plt.show to be a no-op
    import compton as comp_module

    original_show = plt.show
    plt.show = lambda: None
    try:
        comp_module.main()
    finally:
        plt.show = original_show

    fig = plt.gcf()
    yield fig
    plt.close('all')


class TestFigureCreation:
    def test_figure_exists(self, fig_and_widgets):
        """main() produces a figure without raising."""
        fig = fig_and_widgets
        assert fig is not None

    def test_figure_size(self, fig_and_widgets):
        """Figure should be 14×9 inches (allow small floating-point slack)."""
        fig = fig_and_widgets
        w, h = fig.get_size_inches()
        assert abs(w - 14) < 0.5
        assert abs(h - 9) < 0.5

    def test_three_axes(self, fig_and_widgets):
        """Figure must have at least three data axes (+ slider axes)."""
        fig = fig_and_widgets
        # Collect axes that are NOT slider axes (slider axes have no title
        # and are very short).  We count axes whose height > 5% of figure.
        data_axes = [ax for ax in fig.axes
                     if ax.get_position().height > 0.05]
        assert len(data_axes) >= 3, (
            f"Expected ≥3 data axes, found {len(data_axes)}"
        )


class TestSliders:
    def _get_sliders(self, fig):
        import compton
        return list(compton._sliders)

    def test_two_sliders_exist(self, fig_and_widgets):
        """Exactly two Slider widgets should be present."""
        sliders = self._get_sliders(fig_and_widgets)
        assert len(sliders) == 2, (
            f"Expected 2 sliders, found {len(sliders)}"
        )

    def test_theta_slider_range(self, fig_and_widgets):
        """θ slider must span 0° – 180°."""
        sliders = self._get_sliders(fig_and_widgets)
        # The θ slider has a smaller valmax (180) vs λ slider (200)
        theta_sl = min(sliders, key=lambda s: s.valmax)
        assert theta_sl.valmin == pytest.approx(0, abs=1e-3)
        assert theta_sl.valmax == pytest.approx(180, abs=1e-3)

    def test_lambda_slider_range(self, fig_and_widgets):
        """λ₀ slider must span 1 – 200 pm."""
        sliders = self._get_sliders(fig_and_widgets)
        lam_sl = max(sliders, key=lambda s: s.valmax)
        assert lam_sl.valmin == pytest.approx(1, abs=1e-3)
        assert lam_sl.valmax == pytest.approx(200, abs=1e-3)

    def test_theta_default(self, fig_and_widgets):
        """θ slider default is 90°."""
        sliders = self._get_sliders(fig_and_widgets)
        theta_sl = min(sliders, key=lambda s: s.valmax)
        assert theta_sl.val == pytest.approx(90, abs=1e-3)

    def test_lambda_default(self, fig_and_widgets):
        """λ₀ slider default is 70.8 pm (Mo Kα)."""
        sliders = self._get_sliders(fig_and_widgets)
        lam_sl = max(sliders, key=lambda s: s.valmax)
        assert lam_sl.val == pytest.approx(70.8, abs=0.3)


class TestAxesTitles:
    def test_panel_titles_present(self, fig_and_widgets):
        """Each data panel should have a non-empty title."""
        fig = fig_and_widgets
        data_axes = [ax for ax in fig.axes
                     if ax.get_position().height > 0.05]
        titles = [ax.get_title() for ax in data_axes]
        non_empty = [t for t in titles if t.strip()]
        assert len(non_empty) >= 2, (
            f"Expected ≥2 titled panels, got: {titles}"
        )
