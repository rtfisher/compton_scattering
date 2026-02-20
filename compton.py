#!/usr/bin/env python3
"""
Compton Effect — Interactive Visualization
==========================================
Panel A:  Scattering geometry diagram (updated by sliders).
Panel B:  Δλ vs θ theory curve (0–180°) with live marker.
Panel C:  Numeric readout of 8 quantities.
Sliders:  Scattering angle θ · Incident wavelength λ₀.

Physics:
  Δλ = λ' − λ₀ = (h / mₑc)(1 − cos θ) = λ_C (1 − cos θ)
  λ_C ≈ 2.426 pm  (Compton wavelength)
"""

import numpy as np

# ═══════════════════════════════════════════════════════════
#  CONSTANTS
# ═══════════════════════════════════════════════════════════
H_PLANCK = 6.626e-34    # J·s
M_E      = 9.109e-31    # kg
C        = 2.998e8      # m/s
LAMBDA_C = H_PLANCK / (M_E * C)   # ≈ 2.426e-12 m  (Compton wavelength)
EV       = 1.602e-19    # J per eV


# ═══════════════════════════════════════════════════════════
#  PHYSICS HELPERS
# ═══════════════════════════════════════════════════════════
def delta_lambda(theta_deg):
    """Compton wavelength shift Δλ [m] for scattering angle θ [degrees]."""
    theta = np.deg2rad(theta_deg)
    return LAMBDA_C * (1.0 - np.cos(theta))


def scattered_lambda(lambda0, theta_deg):
    """Scattered photon wavelength λ' [m]."""
    return lambda0 + delta_lambda(theta_deg)


def photon_energy_eV(wavelength):
    """Photon energy [eV] for given wavelength [m]."""
    return H_PLANCK * C / (wavelength * EV)


def recoil_angle(lambda0, theta_deg):
    """Electron recoil angle φ [degrees] (from momentum conservation).

    cot φ = (1 + λ_C/λ₀) tan(θ/2)
    Returns 0 for θ = 0 (no scattering).
    """
    if theta_deg <= 0.0:
        return 0.0
    if theta_deg >= 180.0:
        return 0.0
    theta = np.deg2rad(theta_deg)
    cot_phi = (1.0 + LAMBDA_C / lambda0) * np.tan(theta / 2.0)
    phi = np.arctan(1.0 / cot_phi)
    return float(np.rad2deg(phi))


# ═══════════════════════════════════════════════════════════
#  COLOUR HELPERS
# ═══════════════════════════════════════════════════════════
def _wavelength_to_rgb(lam_m):
    """Approximate RGB colour for a wavelength in metres.

    Works in the hard X-ray / UV / visible range.
    For wavelengths < 10 nm (hard X-ray) we return violet.
    """
    lam_nm = lam_m * 1e9
    if lam_nm < 10:       # hard X-ray → deep violet
        return (0.55, 0.0, 0.80)
    if lam_nm < 380:      # soft X-ray / UV → bright violet
        return (0.60, 0.0, 0.85)
    if lam_nm < 440:
        t = (lam_nm - 380) / 60
        return (0.5*(1-t), 0.0, 1.0)
    if lam_nm < 490:
        t = (lam_nm - 440) / 50
        return (0.0, t, 1.0)
    if lam_nm < 510:
        t = (lam_nm - 490) / 20
        return (0.0, 1.0, 1.0 - t)
    if lam_nm < 580:
        t = (lam_nm - 510) / 70
        return (t, 1.0, 0.0)
    if lam_nm < 645:
        t = (lam_nm - 580) / 65
        return (1.0, 1.0 - t, 0.0)
    return (1.0, 0.0, 0.0)  # red / IR


# ═══════════════════════════════════════════════════════════
#  GUI
# ═══════════════════════════════════════════════════════════
# Module-level list populated by main() so tests can inspect the sliders.
_sliders = []


def main():
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.gridspec import GridSpec
    from matplotlib.widgets import Slider, Button
    from matplotlib.collections import LineCollection
    from matplotlib import cm
    import matplotlib.patheffects as pe

    # ── Figure & GridSpec ────────────────────────────────
    fig = plt.figure(figsize=(14, 9))
    fig.patch.set_facecolor('#0d0d1a')
    fig.suptitle('Compton Scattering — Interactive Visualization',
                 color='white', fontsize=16, fontweight='bold', y=0.98)

    gs = GridSpec(
        2, 2,
        figure=fig,
        height_ratios=[1.2, 1],
        width_ratios=[1.4, 1],
        left=0.06, right=0.97,
        top=0.92, bottom=0.22,
        hspace=0.35, wspace=0.30,
    )

    ax_geom  = fig.add_subplot(gs[0, 0])   # Panel A — geometry
    ax_curve = fig.add_subplot(gs[1, :])   # Panel B — Δλ vs θ  (full width)
    ax_info  = fig.add_subplot(gs[0, 1])   # Panel C — numeric readout

    for ax in (ax_geom, ax_curve, ax_info):
        ax.set_facecolor('#060610')
        for sp in ax.spines.values():
            sp.set_color('#333355')

    # ── Panel B: static theory curve ─────────────────────
    theta_arr = np.linspace(0, 180, 360)
    dlam_arr  = delta_lambda(theta_arr) * 1e12   # convert to pm

    # proxy line for legend only (zero width)
    ax_curve.plot([], [], color='#e040fb', lw=2.5,
                  label=r'$\Delta\lambda = \lambda_C\,(1 - \cos\theta)$')
    # gradient LineCollection: plasma colormap, fixed colours (universal curve)
    pts      = np.column_stack([theta_arr, dlam_arr]).reshape(-1, 1, 2)
    segs     = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc_curve = LineCollection(segs,
                              colors=cm.plasma(dlam_arr[:-1] / (2 * LAMBDA_C * 1e12)),
                              lw=2.5, zorder=2)
    ax_curve.add_collection(lc_curve)
    ax_curve.axhline(0,  color='#555577', ls='-',  lw=0.8, alpha=0.5)
    ax_curve.axvline(90, color='#ffffff', ls='--', lw=1.0, alpha=0.35,
                     label=r'$\theta = 90°$')
    ax_curve.axvline(180, color='#ffffff', ls=':', lw=1.0, alpha=0.25)
    ax_curve.annotate(
        r'$\Delta\lambda_{\max} = 2\lambda_C \approx 4.85\,\mathrm{pm}$',
        xy=(180, 2 * LAMBDA_C * 1e12),
        xytext=(130, 2 * LAMBDA_C * 1e12 + 0.35),
        color='#aaaacc', fontsize=9,
        arrowprops=dict(arrowstyle='->', color='#aaaacc', lw=1.0),
    )

    ax_curve.set_xlim(0, 180)
    ax_curve.set_ylim(-0.3, 5.5)
    ax_curve.set_xlabel('Scattering Angle θ (degrees)',
                         color='white', fontsize=10)
    ax_curve.set_ylabel('Δλ (pm)', color='white', fontsize=10)
    ax_curve.tick_params(colors='white', labelsize=9)
    ax_curve.grid(True, color='#1a1a33', lw=0.5)
    ax_curve.legend(fontsize=9, loc='upper left',
                    facecolor='#0d0d1a', edgecolor='#555', labelcolor='white')
    ax_curve.set_title('Wavelength Shift vs Scattering Angle',
                        color='#aaaadd', fontsize=10, pad=4)

    # live dot on Panel B
    live_dot, = ax_curve.plot([], [], 'o', ms=12, color='#ffee00',
                              zorder=6, markeredgecolor='white',
                              markeredgewidth=1.2)

    # ── Panel C: numeric readout ──────────────────────────
    ax_info.set_xlim(0, 1)
    ax_info.set_ylim(0, 1)
    ax_info.axis('off')
    ax_info.set_title('Physics Readout', color='#aaaadd', fontsize=10, pad=4)

    _LABELS = [
        ('θ  (scattering angle)',  0.87),
        ('φ  (recoil angle)',      0.77),
        ('λ₀ (incident)',          0.67),
        ("λ' (scattered)",         0.57),
        ('Δλ (shift)',             0.47),
        ('E₀ (incident)',          0.36),
        ("E' (scattered)",         0.26),
        ('ΔE (energy to e⁻)',      0.15),
    ]
    for lab, yy in _LABELS:
        ax_info.text(0.04, yy, lab, color='#7799bb', fontsize=8.5,
                     family='monospace', va='center')

    # value text objects — updated in callback
    _val_keys = ['theta', 'phi', 'lam0', 'lamp', 'dlam', 'E0', 'Ep', 'dE']
    _val_ypos = [y for _, y in _LABELS]
    vtxt = {}
    for k, yy in zip(_val_keys, _val_ypos):
        vtxt[k] = ax_info.text(0.97, yy, '', color='white', fontsize=9,
                               ha='right', family='monospace',
                               fontweight='bold', va='center')
    ax_info.plot([0.03, 0.97], [0.72, 0.72], color='#222244', lw=0.8)
    ax_info.plot([0.03, 0.97], [0.41, 0.41], color='#222244', lw=0.8)

    # ── Panel A: geometry diagram ────────────────────────
    ax_geom.set_xlim(-1.3, 1.3)
    ax_geom.set_ylim(-1.0, 1.15)
    ax_geom.set_aspect('equal')
    ax_geom.axis('off')
    ax_geom.set_title('Scattering Geometry', color='#aaaadd',
                       fontsize=10, pad=4)

    # electron circle (static)
    e_circle = plt.Circle((0, 0), 0.07, color='#33bbff',
                           zorder=5, lw=1.5, ec='white')
    ax_geom.add_patch(e_circle)
    ax_geom.text(0, -0.16, 'eˉ at rest', color='#33bbff',
                 ha='center', va='top', fontsize=8)

    # horizontal reference axis (dashed, static)
    ax_geom.plot([-1.25, 1.25], [0, 0], color='#333355',
                 ls='--', lw=0.8, zorder=1)

    # arrows — FancyArrow objects (we recreate them on each update)
    _arrow_kw = dict(width=0.022, head_width=0.065, head_length=0.05,
                     length_includes_head=True, zorder=4)

    # store mutable containers so callback can reach them
    state = dict(
        incident_arrow=None,
        scattered_arrow=None,
        recoil_arrow=None,
        arc_patch=None,
        arc_patch2=None,
        lam0_txt=None,
        lamp_txt=None,
        theta_txt=None,
        phi_txt=None,
    )

    def _remove(key):
        obj = state.get(key)
        if obj is not None:
            try:
                obj.remove()
            except Exception:
                pass
            state[key] = None

    def _update_geometry(theta_deg, lambda0):
        """Redraw Panel A for given θ and λ₀."""
        phi_deg = recoil_angle(lambda0, theta_deg)
        theta_rad = np.deg2rad(theta_deg)
        phi_rad   = np.deg2rad(phi_deg)

        col0 = _wavelength_to_rgb(lambda0)
        lamp = scattered_lambda(lambda0, theta_deg)
        colp = _wavelength_to_rgb(lamp)

        L = 1.05   # arrow length

        # ── incident photon (left → origin, horizontal) ──
        _remove('incident_arrow')
        arr = ax_geom.annotate(
            '', xy=(0, 0), xytext=(-L, 0),
            arrowprops=dict(arrowstyle='->', color=col0,
                            lw=2.5, mutation_scale=18),
            zorder=4,
        )
        state['incident_arrow'] = arr

        # ── scattered photon (origin → upper right at θ) ──
        _remove('scattered_arrow')
        sx = L * np.cos(theta_rad)
        sy = L * np.sin(theta_rad)
        arr2 = ax_geom.annotate(
            '', xy=(sx, sy), xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', color=colp,
                            lw=2.5, mutation_scale=18),
            zorder=4,
        )
        state['scattered_arrow'] = arr2

        # ── recoil electron (origin → lower right at −φ) ──
        _remove('recoil_arrow')
        rx = L * 0.80 * np.cos(-phi_rad)
        ry = L * 0.80 * np.sin(-phi_rad)
        arr3 = ax_geom.annotate(
            '', xy=(rx, ry), xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', color='#888899',
                            lw=1.8, mutation_scale=14),
            zorder=3,
        )
        state['recoil_arrow'] = arr3

        # ── arc for θ ──
        _remove('arc_patch')
        arc_r = 0.28
        theta1_arc = 0
        theta2_arc = theta_deg
        arc = mpatches.Arc((0, 0), 2*arc_r, 2*arc_r,
                            angle=0,
                            theta1=theta1_arc, theta2=theta2_arc,
                            color='#ffee00', lw=1.4, zorder=5)
        ax_geom.add_patch(arc)
        state['arc_patch'] = arc

        # ── arc for φ ──
        _remove('arc_patch2')
        arc2_r = 0.22
        arc2 = mpatches.Arc((0, 0), 2*arc2_r, 2*arc2_r,
                             angle=0,
                             theta1=-phi_deg, theta2=0,
                             color='#aaaaaa', lw=1.2, zorder=5)
        ax_geom.add_patch(arc2)
        state['arc_patch2'] = arc2

        # ── text labels ──
        for k in ('lam0_txt', 'lamp_txt', 'theta_txt', 'phi_txt'):
            _remove(k)

        # λ₀ label above incident arrow
        state['lam0_txt'] = ax_geom.text(
            -0.55, 0.10,
            f'λ₀ = {lambda0*1e12:.2f} pm',
            color=col0, fontsize=8.5, ha='center', va='bottom',
            fontweight='bold',
        )

        # λ' label along scattered arrow
        mid_sx = sx * 0.55
        mid_sy = sy * 0.55
        perp_x = -np.sin(theta_rad) * 0.13
        perp_y =  np.cos(theta_rad) * 0.13
        state['lamp_txt'] = ax_geom.text(
            mid_sx + perp_x, mid_sy + perp_y,
            f"λ' = {lamp*1e12:.2f} pm",
            color=colp, fontsize=8.5, ha='center', va='bottom',
            fontweight='bold',
        )

        # θ label inside arc
        theta_mid = np.deg2rad(theta_deg / 2)
        state['theta_txt'] = ax_geom.text(
            (arc_r + 0.10) * np.cos(theta_mid),
            (arc_r + 0.10) * np.sin(theta_mid),
            f'θ={theta_deg:.0f}°',
            color='#ffee00', fontsize=8, ha='center', va='center',
        )

        # φ label
        phi_mid = np.deg2rad(-phi_deg / 2)
        if phi_deg > 1.0:
            state['phi_txt'] = ax_geom.text(
                (arc2_r + 0.10) * np.cos(phi_mid),
                (arc2_r + 0.10) * np.sin(phi_mid),
                f'φ={phi_deg:.1f}°',
                color='#aaaaaa', fontsize=7.5, ha='center', va='center',
            )

    # ── Sliders ──────────────────────────────────────────
    skw = dict(facecolor='#0d0d24')
    ax_sl_theta = fig.add_axes([0.12, 0.14, 0.34, 0.025], **skw)
    ax_sl_lam   = fig.add_axes([0.57, 0.14, 0.34, 0.025], **skw)

    sl_theta = Slider(ax_sl_theta, 'Angle θ (°)', 0, 180,
                      valinit=90, valstep=1, color='#5577ff')
    sl_lam   = Slider(ax_sl_lam,   'Incident λ₀ (pm)', 1, 200,
                      valinit=70.8, valstep=0.2, color='#cc4444')

    for sl in (sl_theta, sl_lam):
        sl.label.set_color('white')
        sl.valtext.set_color('white')

    # Expose for testing
    _sliders.clear()
    _sliders.extend([sl_theta, sl_lam])

    fig.text(0.12, 0.17, 'Scattering angle', color='#7799bb',
             fontsize=8, ha='left')
    fig.text(0.57, 0.17, 'Incident wavelength',
             color='#7799bb', fontsize=8, ha='left')

    # ── Preset source buttons ─────────────────────────────
    PRESETS = [
        ('Mo Kα\n70.8 pm',  70.8),
        ('Cu Kα\n154.2 pm', 154.2),
        ('Ag Kα\n56.1 pm',  56.1),
        ('Fe Kα\n193.7 pm', 193.7),
        ('Co Kα\n178.9 pm', 178.9),
    ]

    fig.text(0.10, 0.100, 'Presets:', color='#7799bb', fontsize=8,
             ha='left', va='center')

    preset_axes_list = []
    preset_btn_list  = []
    selected_preset  = [0]

    def _update_preset_highlight(idx):
        for i, pax in enumerate(preset_axes_list):
            pax.set_facecolor('#2a3a6a' if i == idx else '#1a1a3a')

    btn_x_starts = [0.10, 0.26, 0.42, 0.58, 0.74]
    for i, ((label, lam_pm), bx) in enumerate(zip(PRESETS, btn_x_starts)):
        pax = fig.add_axes([bx, 0.055, 0.15, 0.040])
        pax.set_facecolor('#1a1a3a')
        btn = Button(pax, label, color='#1a1a3a', hovercolor='#2a3a6a')
        btn.label.set_color('white')
        btn.label.set_fontfamily('monospace')
        btn.label.set_fontsize(8)
        preset_axes_list.append(pax)
        preset_btn_list.append(btn)

    def make_preset_cb(idx, lam_pm):
        def cb(event):
            selected_preset[0] = idx
            _update_preset_highlight(idx)
            sl_lam.set_val(lam_pm)   # triggers on_slider_change automatically
        return cb

    for i, (_, lam_pm) in enumerate(PRESETS):
        preset_btn_list[i].on_clicked(make_preset_cb(i, lam_pm))

    # ── Callback ─────────────────────────────────────────
    def on_slider_change(_):
        theta_deg = float(sl_theta.val)
        lam0_pm   = float(sl_lam.val)
        lam0      = lam0_pm * 1e-12   # m

        dlam   = delta_lambda(theta_deg)
        lamp   = scattered_lambda(lam0, theta_deg)
        E0_eV  = photon_energy_eV(lam0)
        Ep_eV  = photon_energy_eV(lamp)
        dE_eV  = E0_eV - Ep_eV
        phi_deg = recoil_angle(lam0, theta_deg)

        # Panel B — live dot
        live_dot.set_data([theta_deg], [dlam * 1e12])

        # Panel C — numeric readout
        vtxt['theta'].set_text(f'{theta_deg:.1f} °')
        vtxt['phi'].set_text(f'{phi_deg:.1f} °')
        vtxt['lam0'].set_text(f'{lam0_pm:.3f} pm')
        vtxt['dlam'].set_text(f'{dlam*1e12:.3f} pm')
        vtxt['lamp'].set_text(f'{lamp*1e12:.3f} pm')
        vtxt['E0'].set_text(f'{E0_eV/1e3:.3f} keV')
        vtxt['Ep'].set_text(f'{Ep_eV/1e3:.3f} keV')
        vtxt['dE'].set_text(f'{dE_eV/1e3:.3f} keV')

        # Panel A — geometry
        _update_geometry(theta_deg, lam0)

        fig.canvas.draw_idle()

    sl_theta.on_changed(on_slider_change)
    sl_lam.on_changed(on_slider_change)

    # ── Initial draw ──────────────────────────────────────
    on_slider_change(None)
    _update_preset_highlight(0)   # highlight Mo Kα at startup

    plt.show()


if __name__ == '__main__':
    main()
