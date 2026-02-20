# Plan: Compton Effect Visualization

## Context

Build an interactive visualization of Compton scattering for sophomore modern physics (PHY 213). A photon scatters from a free electron at rest; the scattered photon emerges at angle θ with a longer wavelength. The key relationship is:

```
Δλ = λ' − λ₀ = (h / mₑc)(1 − cos θ) = λ_C (1 − cos θ)
```

where λ_C ≈ 2.426 pm is the Compton wavelength. The visualization should be self-contained, slider-driven, and match the dark-themed matplotlib style used in the sibling `_photoelectric` project.

---

## Technology Stack

- **Python 3.10+**, **NumPy**, **Matplotlib ≥ 3.5** — identical to the photoelectric project
- `matplotlib.widgets.Slider` for interactive controls
- No additional dependencies (no scipy, no bokeh, no panel)

---

## File Structure

```
_compton/
├── compton.py          # Main script: physics + layout + interaction
├── PLAN.md             # This file
├── README.md
├── .gitignore
└── tests/
    ├── __init__.py
    ├── test_physics.py        # Unit tests for Compton physics functions
    └── test_visualization.py  # Smoke tests for figure/slider construction
```

---

## Physics Module (top of `compton.py`)

Constants:
```python
H_PLANCK  = 6.626e-34    # J·s
M_E       = 9.109e-31    # kg
C         = 2.998e8      # m/s
LAMBDA_C  = H_PLANCK / (M_E * C)  # ≈ 2.426e-12 m (Compton wavelength)
EV        = 1.602e-19    # J per eV
```

Key functions:
- `delta_lambda(theta_deg)` → Δλ in metres = λ_C * (1 − cos θ)
- `scattered_lambda(lambda0, theta_deg)` → λ' = λ₀ + Δλ
- `photon_energy_eV(wavelength)` → E = hc/λ in eV
- `recoil_angle(lambda0, theta_deg)` → φ (electron recoil angle)

---

## Layout: Three-Panel Figure

Canvas: **14 × 9 inches**, dark background `#0d0d1a`

```
┌───────────────────────────┬──────────────────┐
│                           │                  │
│   Panel A: Scattering     │  Panel C:        │
│   Geometry Diagram        │  Numeric Readout │
│   (static arrows/labels,  │  (8 quantities)  │
│    updated on slider)     │                  │
│                           │                  │
├───────────────────────────┴──────────────────┤
│   Panel B: Δλ vs θ Theory Curve              │
│   (0–180°, live marker at current θ)         │
└──────────────────────────────────────────────┘
        [Angle slider]  [Incident λ slider]
```

GridSpec: 2 rows × 2 cols, height_ratios=[1.2, 1], width_ratios=[1.4, 1]

### Panel A — Scattering Diagram (top-left)

Static reference frame; elements updated on each slider change:

- **Incident photon**: horizontal arrow from left → electron, colored by λ₀ (violet–red X-ray scale)
- **Electron**: filled circle at origin, labeled "eˉ at rest"
- **Scattered photon**: arrow at angle θ above horizontal, colored by λ' (slightly redder)
- **Recoil electron**: arrow at angle φ below horizontal (gray, secondary)
- Arc drawn at origin showing angle θ with label
- Text overlays: "λ₀", "λ'", "θ", "φ"
- Numeric values update when sliders change

### Panel B — Δλ vs θ Curve (bottom, full width)

- Solid cyan curve of Δλ(θ) for θ ∈ [0°, 180°]
- Horizontal dashed line at Δλ = 0 (reference)
- Vertical dashed lines at θ = 90° and θ = 180°
- **Live yellow dot** at (θ_current, Δλ_current) — updates on slider
- Y-axis in picometres (× 10¹²)
- Annotation: "Δλ_max = 2λ_C ≈ 4.85 pm" at θ = 180°

### Panel C — Numeric Readout (top-right)

Monospace text block, 8 quantities:

```
θ (scattering angle)   = XX.X °
φ (recoil angle)       = XX.X °
λ₀ (incident)         = X.XXX pm
λ' (scattered)        = X.XXX pm
Δλ (shift)            = X.XXX pm
E₀ (incident)         = XX.XX keV
E' (scattered)        = XX.XX keV
ΔE (energy to e⁻)     = X.XXX keV
```

---

## Interactive Controls

Two `matplotlib.widgets.Slider` bars below the figure:

| Slider | Color | Range | Default |
|--------|-------|-------|---------|
| Scattering angle θ | `#5577ff` (blue) | 0° – 180° | 90° |
| Incident wavelength λ₀ | `#cc4444` (red) | 1 pm – 200 pm | 70.8 pm (Mo Kα) |

Callback `on_slider_change`:
1. Read θ and λ₀ from sliders
2. Recompute all physics quantities
3. Update Panel A arrow directions and text labels
4. Move the live dot in Panel B
5. Refresh Panel C text block
6. Call `fig.canvas.draw_idle()`

---

## Tests

### `tests/test_physics.py`
- `delta_lambda(0)` == 0
- `delta_lambda(90)` ≈ λ_C ≈ 2.426 pm
- `delta_lambda(180)` ≈ 2λ_C ≈ 4.852 pm
- `scattered_lambda(lambda0, theta)` ≥ `lambda0` for all θ
- `photon_energy_eV` correct for known wavelengths
- Compton formula matches textbook values at θ = 0°, 90°, 180°

### `tests/test_visualization.py`
- Figure creation does not raise
- Two sliders present with correct ranges
- Three axes exist

---

## Verification

Run the visualization:
```bash
cd _compton && python compton.py
```

Manual checks:
- θ = 0°: Δλ = 0, λ' = λ₀, no deflection in diagram
- θ = 90°: Δλ ≈ 2.43 pm, dot at midpoint of curve
- θ = 180°: Δλ ≈ 4.85 pm, scattered photon points back
- Changing λ₀: curve shape unchanged (Δλ independent of λ₀), keV values change

Run tests:
```bash
pytest tests/ -v
```
