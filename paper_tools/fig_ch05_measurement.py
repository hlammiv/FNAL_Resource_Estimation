#!/usr/bin/env python3
"""Ch. 5: the 2033 two-time measurement design, from the Ch. 5 model (markup board Q5-6, 2026-11-09).

Normalized C_{k_min}(t) for the planning decay band tau = 1/k (fast end, T = 300 MeV) to 1/(2 pi T) (slow end,
T = 225 MeV), the two fit points at Gamma t = 0.5 and 1.8 on each edge (rounded up onto the Delta t = 0.01 fm/c grid,
with the first result's relative errors), the dropped early transient Gamma t < 0.5, the hydrodynamic decay time the
chapter says has no basis at this k, and the free-streaming correlator 3 (sin kt - kt cos kt) / (kt)^3, which reaches
1/e only at kt = 2.9. Every number is read from the model; the asserts pin the values the chapter prints.
"""
import numpy as np
import matplotlib.pyplot as plt

import estimates
import fnalfig as ff


def free_streaming(x):
    """Transverse momentum-density correlator of free-streaming massless quanta, x = k t (model, T1 verifier)."""
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    nz = x > 1e-9
    out[nz] = 3 * (np.sin(x[nz]) - x[nz] * np.cos(x[nz])) / x[nz] ** 3
    return out


FS = 7.0     # smallest label; the canvas is the printed size (full text width), so this prints at 7pt


def main():
    ff.setup(truetype=True, tight=False)
    m = estimates.load(5)
    a = m.Assumptions()
    i = m.model(a, "2033").intermediates

    tau_fast, tau_slow = i["tau_band_fm"]                     # 0.0955 (1/k), 0.1396 fm/c (1/(2 pi T) at 225 MeV)
    k_inv_fm = 1.0 / i["tau_free_streaming_fm"]               # k_min in fm^-1
    dt = a.a_t_fm.lo * a.dt_over_at_2033.lo                   # 0.01 fm/c grid
    steps_fast, steps_slow = i["fit_steps_grid"]              # (5, 18), (7, 26)
    x_fast, x_slow = i["fit_x_realized"]                      # Gamma t at the grid times
    err_fast, err_slow = i["corr_rel_err_first_at_fit_times"] # first result: 0.186, 0.364 at the two times
    x_cut = a.fit_x.lo                                        # the fit drops Gamma t < 0.5
    frac_early = a.fit_frac_early.lo                          # 20% of the shots at the early time
    fs_kt, fs_fm = i["free_streaming_1e_kt"], i["free_streaming_1e_fm"]
    hyd_lo, hyd_hi = sorted(i["tau_decay_at_eta_s_assumed_fm"])

    # the values the chapter prints
    # staggered ruling 2026-11-13 (2^3, k_min a_s = pi): tau_fast 0.064 (was 0.0955 on 3^3), grid (4, 12) (was (5, 18)),
    # late-time error 0.33 at the fast end (was 0.36), free-streaming 1/e at 0.19 fm/c (was 0.28), hydro record 0.03-0.04
    assert round(tau_fast, 2) == 0.06 and round(tau_slow, 2) == 0.14
    assert (a.fit_x.lo, a.fit_x.hi) == (0.5, 1.8) and frac_early == 0.2
    assert steps_fast == (4, 12) and steps_slow == (7, 26) and abs(dt - 0.01) < 1e-12
    assert all(abs(x - xp) < 0.15 for xs in (x_fast, x_slow) for x, xp in zip(xs, (0.5, 1.8)))   # 2^3: fast grid lands at 0.63, 1.88
    assert round(err_fast[0], 2) == 0.19 and round(err_fast[1], 2) == 0.33
    assert round(err_slow[0], 2) == 0.19 and round(err_slow[1], 2) == 0.36
    assert round(fs_kt, 1) == 2.9 and round(fs_fm, 2) == 0.19
    assert round(hyd_lo, 2) == 0.03 and round(hyd_hi, 2) == 0.04
    assert round(i["design_A_over_B"], 2) == 0.37 and round(i["design_C_over_B"], 1) == 2.7

    t = np.linspace(0, 0.32, 400)
    c_fast, c_slow = np.exp(-t / tau_fast), np.exp(-t / tau_slow)
    e_inv = np.exp(-1.0)

    fig, ax = plt.subplots(figsize=(ff.TEXTWIDTH_IN, 2.47))

    # dropped early transient: Gamma t < 0.5 is t < 0.5 tau, 0.048 fm/c at the fast end to 0.070 at the slow end
    ax.axvspan(0, x_cut * tau_fast, color=ff.PALEGREY, alpha=0.55, lw=0)
    ax.axvspan(x_cut * tau_fast, x_cut * tau_slow, color=ff.PALEGREY, alpha=0.25, lw=0)
    ax.text(0.037, 0.012, "dropped from the fit\n($\\Gamma t<0.5$)", ha="center", va="bottom",
            fontsize=FS, color=ff.DARKGREY, linespacing=1.1)

    # planning decay band and its two edges
    ax.fill_between(t, c_fast, c_slow, color=ff.BLUE, alpha=0.15, lw=0)
    ax.plot(t, c_fast, color=ff.BLUE, lw=0.9)
    ax.plot(t, c_slow, color=ff.BLUE, lw=0.9, ls=(0, (4, 2)))
    ax.text(0.232, 0.17, "$\\tau=1/(2\\pi T)$, $T=225$ MeV", color=ff.BLUE, fontsize=FS, ha="left", va="bottom",
            rotation=-17)
    ax.text(0.085, 0.03, "$\\tau=1/k$, $T=300$ MeV", color=ff.BLUE, fontsize=FS, ha="left", va="bottom")
    ax.text(0.135, 0.60, "planning band $e^{-t/\\tau}$", color=ff.BLUE, fontsize=FS, ha="center", va="center")

    # free streaming, 1/e only at kt = 2.9
    c_fs = free_streaming(k_inv_fm * t)
    ax.plot(t, c_fs, color=ff.ORANGE, lw=0.9)
    ax.text(0.168, 0.86, "free streaming\n$3(\\sin kt-kt\\cos kt)/(kt)^3$", color=ff.ORANGE, fontsize=FS,
            ha="left", va="bottom", linespacing=1.1)
    ax.plot([fs_fm], [e_inv], "o", color=ff.ORANGE, ms=3.5)
    ax.annotate(f"$kt={fs_kt:.1f}$\n(${fs_fm:.2f}$ fm/$c$)", xy=(fs_fm, e_inv), xytext=(4, 3),
                textcoords="offset points", fontsize=FS, color=ff.ORANGE, ha="left", va="bottom", linespacing=1.1)

    # the 1/e line, with the hydrodynamic tau that has no basis at this k
    ax.axhline(e_inv, color=ff.GREY, lw=0.6, ls=":")
    ax.text(0.318, e_inv + 0.012, "$1/e$", color=ff.GREY, fontsize=FS, ha="right", va="bottom")
    ax.plot([hyd_lo, hyd_hi], [e_inv, e_inv], color=ff.DARKGREY, lw=2.2, solid_capstyle="butt", zorder=4)
    ax.text(0.002, 0.335,
            f"$\\tau_{{\\rm hydro}}=T/((\\eta/s)k^2)$\n$\\eta/s={a.eta_s_assumed.lo}$: "
            f"${hyd_lo:.2f}$--${hyd_hi:.2f}$ fm/$c$\n(no basis at this $k$)",
            fontsize=FS, color=ff.DARKGREY, ha="left", va="top", linespacing=1.0)

    # the two fit points on each edge, rounded up onto the 0.01 fm/c grid, with the first result's errors
    for steps, xs, errs, tau, filled in ((steps_fast, x_fast, err_fast, tau_fast, True),
                                         (steps_slow, x_slow, err_slow, tau_slow, False)):
        tt = np.array(steps) * dt
        cc = np.exp(-np.array(xs))
        assert np.allclose(cc, np.exp(-tt / tau))
        ax.errorbar(tt, cc, yerr=np.array(errs) * cc, fmt="o", color=ff.BLUE, ms=3.8, capsize=2, lw=0.8,
                    mfc=ff.BLUE if filled else "white", mec=ff.BLUE, zorder=5)
    ax.annotate(f"$\\Gamma t=0.5$\n({int(round(100 * frac_early))}% of shots)",
                xy=(steps_slow[0] * dt, np.exp(-x_slow[0])), xytext=(0.084, 0.72), fontsize=FS, color=ff.INK,
                ha="left", va="center", linespacing=1.1,
                arrowprops=dict(arrowstyle="-", color=ff.GREY, lw=0.5, shrinkA=0, shrinkB=3))
    # one line, inside the band left of the free-streaming curve, with a leader as for Gamma t = 0.5
    ax.annotate(f"$\\Gamma t=1.8$ ({int(round(100 * (1 - frac_early)))}%)",
                xy=(steps_slow[1] * dt, np.exp(-x_slow[1])), xytext=(0.226, 0.085),
                fontsize=FS, color=ff.INK, ha="right", va="center",
                arrowprops=dict(arrowstyle="-", color=ff.GREY, lw=0.5, shrinkA=0, shrinkB=3))

    ax.set_xlim(0, 0.32)
    ax.set_ylim(0, 1.04)
    ax.set_xlabel("$t$ (fm/$c$)")
    ax.set_ylabel("$C_{k_{\\min}}(t)\\,/\\,C_{k_{\\min}}(0)$")
    ff.despine(ax)
    ax.tick_params(top=False, right=False)
    fig.tight_layout(pad=0.2)
    ff.save(fig, "ch05_measurement")


if __name__ == "__main__":
    main()
