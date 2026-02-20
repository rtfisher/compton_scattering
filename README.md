# Compton Scattering — Interactive Visualization

An interactive matplotlib tool for PHY 213 (Modern Physics) that lets students
explore the Compton effect in real time.  Drag the sliders or click a preset
X-ray source button to instantly update the scattering geometry, the Δλ–θ
curve, and the numeric readout.

---

## Physics background

When an X-ray photon scatters off a free electron at rest, its wavelength
increases by an amount that depends only on the scattering angle θ:

$$\Delta\lambda = \lambda' - \lambda_0 = \frac{h}{m_e c}(1 - \cos\theta) = \lambda_C(1 - \cos\theta)$$

where **λ_C ≈ 2.426 pm** is the Compton wavelength of the electron.  Key
features:

- The shift is **independent of λ₀** — it is a purely quantum-mechanical result.
- Maximum shift (θ = 180°, backscatter): **Δλ_max = 2λ_C ≈ 4.85 pm**.
- The recoil electron carries the energy difference ΔE = E₀ − E′.

---

## Interface

```
┌────────────────────────────────────────────────────────┐
│  Panel A: Scattering Geometry   │  Panel C: Readout    │
│  (animated arrows + arcs)       │  (8 quantities)      │
├─────────────────────────────────┴──────────────────────┤
│  Panel B: Δλ vs θ  (plasma-gradient curve + live dot)  │
├────────────────────────────────────────────────────────┤
│  Slider: Angle θ (°)    │  Slider: Incident λ₀ (pm)   │
├────────────────────────────────────────────────────────┤
│  Presets: [Mo Kα] [Cu Kα] [Ag Kα] [Fe Kα] [Co Kα]   │
└────────────────────────────────────────────────────────┘
```

### Panel A — Scattering geometry
Coloured arrows (hue tracks approximate photon energy) show the incident
photon, scattered photon, and recoil electron.  Arc labels update live with
the current θ and φ values.

### Panel B — Δλ vs θ curve
A **plasma-colormap gradient** runs from dark purple (θ = 0°, no shift)
through orange (θ = 90°, Δλ = λ_C) to bright yellow (θ = 180°, maximum
shift).  The gradient colours are fixed — a reminder that the Compton curve is
universal (independent of λ₀).  A yellow dot marks the current (θ, Δλ) point.

### Panel C — Numeric readout
| Symbol | Quantity |
|--------|----------|
| θ | Scattering angle |
| φ | Electron recoil angle |
| λ₀ | Incident wavelength |
| λ′ | Scattered wavelength |
| Δλ | Wavelength shift |
| E₀ | Incident photon energy (keV) |
| E′ | Scattered photon energy (keV) |
| ΔE | Energy transferred to electron (keV) |

---

## Controls

| Control | Range | Default |
|---------|-------|---------|
| Scattering angle θ | 0° – 180° (1° steps) | 90° |
| Incident wavelength λ₀ | 1 – 200 pm (0.2 pm steps) | 70.8 pm (Mo Kα) |

### Preset X-ray source buttons

| Button | Source | λ₀ | E₀ |
|--------|--------|-----|-----|
| Mo Kα | Molybdenum K-alpha | 70.8 pm | 17.5 keV |
| Cu Kα | Copper K-alpha | 154.2 pm | 8.0 keV |
| Ag Kα | Silver K-alpha | 56.1 pm | 22.1 keV |
| Fe Kα | Iron K-alpha | 193.7 pm | 6.4 keV |
| Co Kα | Cobalt K-alpha | 178.9 pm | 6.9 keV |

Clicking a preset snaps the λ₀ slider to that value and highlights the active
button (blue-tinted background).

---

## Requirements

| Package | Version tested |
|---------|---------------|
| Python | 3.11 |
| NumPy | 2.3 |
| Matplotlib | 3.10 |

A `conda` environment called **`npscipy`** is assumed.  To create it from
scratch:

```bash
conda create -n npscipy python=3.11 numpy matplotlib
conda activate npscipy
```

---

## Running

```bash
conda run -n npscipy python compton.py
```

Or activate the environment first:

```bash
conda activate npscipy
python compton.py
```

---

## Tests

52 unit tests cover all physics functions and the visualization layer.

```bash
conda run -n npscipy python -m pytest tests/ -v
```

Expected output: **52 passed**.

### Test coverage summary

| Module | Tests |
|--------|-------|
| `delta_lambda` | Zero angle, 90°, 180°, known pm values, monotonicity, parametrised formula check |
| `scattered_lambda` | λ′ ≥ λ₀ for all θ (parametrised), special cases, formula consistency |
| `photon_energy_eV` | Visible red/violet, Mo Kα, inverse-proportionality, hc/λ formula |
| `recoil_angle` | Zero / 180° edge cases, 90° bounds, forward-hemisphere invariant |
| Visualization | Figure exists, size, three axes, two sliders, ranges, defaults, panel titles |

---

## Project structure

```
_compton/
├── compton.py          # Single-file application (physics + GUI)
├── tests/
│   ├── test_physics.py         # Physics function unit tests
│   └── test_visualization.py   # GUI / slider tests
└── README.md
```

All physics helpers (`delta_lambda`, `scattered_lambda`, `photon_energy_eV`,
`recoil_angle`) live at module level and are imported directly by the test
suite.  The GUI is isolated inside `main()`, which is only called when the
script is run directly.

---

## Related

This tool is a companion to the **Photoelectric Effect** interactive
visualization used in the same course (`_photoelectric`), and follows the same
dark-theme design language (background `#0d0d1a`, panel `#060610`).

---

*PHY 213 — Modern Physics · Spring 2026*
