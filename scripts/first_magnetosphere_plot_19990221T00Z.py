#!/usr/bin/env python3
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

import spacepy.time as spt
import spacepy.omni as om
import spacepy.toolbox as tb

from geopack import geopack as gp


def unix_seconds(t_utc: dt.datetime) -> float:
    """Seconds since 1970-01-01T00:00:00Z (what geopack.recalc expects)."""
    if t_utc.tzinfo is None:
        t_utc = t_utc.replace(tzinfo=dt.timezone.utc)
    epoch = dt.datetime(1970, 1, 1, tzinfo=dt.timezone.utc)
    return (t_utc - epoch).total_seconds()


def get_omni2_params(t_iso: str):
    """
    Pull OMNI2 hourly values at the requested time.

    SpacePy's OMNI2 naming varies by version; on your install we use:
      Flow_pressure (nPa), Dst_index (nT), By_GSM (nT), Bz_GSM (nT)
    """
    ticks = spt.Ticktock([t_iso], "ISO")

    # Ensure OMNI2 cache exists locally (this is *not* Qin–Denton)
    try:
        omni = om.get_omni(ticks, dbase="OMNI2hourly")
    except Exception as e:
        print("Could not load OMNI2 via SpacePy. Trying to download/update OMNI2 now...")
        print(f"Original error: {e}\n")
        tb.update(omni2=True)
        omni = om.get_omni(ticks, dbase="OMNI2hourly")

    # Helper: find the first existing key among candidates
    def pick_key(candidates):
        for k in candidates:
            if k in omni:
                return k
        raise KeyError(f"None of these keys were found: {candidates}. Keys present include: {list(omni.keys())[:40]} ...")

    k_pdyn = pick_key(["Flow_pressure", "Pdyn", "Pressure", "dynamic_pressure"])
    k_dst  = pick_key(["Dst_index", "Dst", "DST"])
    k_by   = pick_key(["By_GSM", "ByIMF", "ByIMF_GSM"])
    k_bz   = pick_key(["Bz_GSM", "BzIMF", "BzIMF_GSM"])
    k_ma   = pick_key(["Alfven_mach_number"])

    Pdyn = float(omni[k_pdyn][0])
    Dst  = float(omni[k_dst][0])
    By   = float(omni[k_by][0])
    Bz   = float(omni[k_bz][0])
    Ma   = float(omni[k_ma][0])

    return Pdyn, Dst, By, Bz, Ma


def bowshock_xy_from_Ma(r0_mp, Ma, ntheta=361, e=0.9, gamma=5/3,
                        delta_min=1.0, delta_max=8.0):
    """
    Bow shock curve where the *nose standoff* is set by Alfvén Mach number Ma.

    r0_mp : magnetopause subsolar distance (Re) from Shue r0
    Ma    : Alfvén Mach number (must be > 1 for a shock)
    e     : conic shaping parameter (0.8–1.1 is reasonable)
    gamma : ratio of specific heats

    delta_min/max : clamp the sheath thickness at the nose (Re) so it doesn't blow up.
    Returns x, y, Rbs0, delta
    """

    # Guard against unphysical / fill values
    if (not np.isfinite(Ma)) or (Ma <= 1.0):
        # No proper bow shock if Ma <= 1; fall back to a modest offset
        delta = 2.5
    else:
        M2 = Ma * Ma
        # Farris & Russell style standoff scaling (captures Ma dependence)
        frac = ((gamma - 1.0) * M2 + 2.0) / ((gamma + 1.0) * (M2 - 1.0))
        delta = r0_mp * frac

        # Clamp to keep plots sane for Ma near 1
        delta = float(np.clip(delta, delta_min, delta_max))

    Rbs0 = r0_mp + delta

    eps = 1e-6
    theta = np.linspace(0.0, np.pi - eps, ntheta)

    # Scale conic so that r(0) = Rbs0
    r = Rbs0 * (1.0 + e) / (1.0 + e * np.cos(theta))

    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, Rbs0, delta

def simple_bowshock_xy(r0_mp, ntheta=361, e=0.9, scale=1.35):
    """
    Simple bow shock curve as a conic section.

    r0_mp : magnetopause standoff distance (Re)
    scale : bow shock nose / magnetopause nose (typ. 1.25–1.5)
    e     : eccentricity-like shaping parameter (0.8–1.1 is reasonable)

    Returns x, y, Rbs0 where Rbs0 is the bow shock subsolar distance.
    """
    Rbs0 = scale * r0_mp  # desired subsolar bow shock distance

    eps = 1e-6
    theta = np.linspace(0.0, np.pi - eps, ntheta)

    # Scale so that r(theta=0) = Rbs0
    r = Rbs0 * (1.0 + e) / (1.0 + e * np.cos(theta))

    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, Rbs0

def t96_parmod(Pdyn, Dst, By, Bz):
    par = np.zeros(10, dtype=float)
    par[0] = Pdyn
    par[1] = Dst
    par[2] = By
    par[3] = Bz
    return par


def shue98_magnetopause_xy(Pdyn, Bz, ntheta=361):
    """
    Shue et al. (1998) magnetopause model:
      r(θ) = r0 * (2/(1+cosθ))^α
      r0 = (10.22 + 1.29*tanh(0.184*(Bz + 8.14))) * Pdyn^(-1/6.6)
      α  = (0.58 - 0.007*Bz) * (1 + 0.024*ln(Pdyn))
    (Pdyn in nPa, Bz in nT).  Returns x,y in the GSM equatorial plane (Re).
    """
    # Protect against weird/missing values
    if not np.isfinite(Pdyn) or Pdyn <= 0:
        raise ValueError(f"Nonpositive/invalid Pdyn={Pdyn}")
    if not np.isfinite(Bz):
        raise ValueError(f"Invalid Bz={Bz}")

    r0 = (10.22 + 1.29 * np.tanh(0.184 * (Bz + 8.14))) * (Pdyn ** (-1.0 / 6.6))
    alpha = (0.58 - 0.007 * Bz) * (1.0 + 0.024 * np.log(Pdyn))

    eps = 1e-6
    theta = np.linspace(0, np.pi-eps, ntheta)  # 0: +X (subsolar) -> pi: tail
    r = r0 * (2.0 / (1.0 + np.cos(theta))) ** alpha

    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, r0, alpha


def trace_field_line(seed_xyz_gsm, parmod, ut_seconds, rlim=50.0, r0=1.0):
    """
    Trace a field line from seed_xyz_gsm in GSM using external model T96 + internal IGRF.
    geopack.trace returns: xf,yf,zf, xx,yy,zz  (endpoint + full path arrays).
    """
    ps = gp.recalc(ut_seconds)  # updates tilt etc. in geopack's global state

    x0, y0, z0 = seed_xyz_gsm

    # dir=-1 and dir=+1 traces in opposite directions along the field
    xf1, yf1, zf1, xx1, yy1, zz1 = gp.trace(x0, y0, z0, -1, rlim, r0, parmod, "t96", "igrf")
    xf2, yf2, zf2, xx2, yy2, zz2 = gp.trace(x0, y0, z0,  1, rlim, r0, parmod, "t96", "igrf")

    # Stitch into one continuous line through the seed point
    X = np.concatenate([xx1[::-1], xx2[1:]])
    Y = np.concatenate([yy1[::-1], yy2[1:]])
    Z = np.concatenate([zz1[::-1], zz2[1:]])

    return X, Y, Z, ps


def main():
    t_iso = "1999-02-21T00:00:00"
    t_utc = dt.datetime(1999, 2, 21, 0, 0, 0, tzinfo=dt.timezone.utc)

    Pdyn, Dst, By, Bz, Ma = get_omni2_params(t_iso)
    parmod = t96_parmod(Pdyn, Dst, By, Bz)
    ut = unix_seconds(t_utc)

    # Magnetopause curve (equatorial plane)
    xmp, ymp, r0mp, alphamp = shue98_magnetopause_xy(Pdyn, Bz, ntheta=361)
    #xbs, ybs, rbs = simple_bowshock_xy(r0mp)
    xbs, ybs, rbs, delta_bs = bowshock_xy_from_Ma(r0mp, Ma)
    
    # Field lines in meridian plane (Y=0 seeds)
    seeds = [
        (3.0, 0.0,  0.0),
        (3.0, 0.0,  1.0),
        (3.0, 0.0,  2.0),
        (3.0, 0.0,  3.0),
        (3.0, 0.0, -1.0),
        (3.0, 0.0, -2.0),
        (3.0, 0.0, -3.0),
        (-3.0, 0.0,  0.0),
        (-3.0, 0.0,  1.0),
        (-3.0, 0.0,  2.0),
        (-3.0, 0.0,  3.0),
        (-3.0, 0.0, -1.0),
        (-3.0, 0.0, -2.0),
        (-3.0, 0.0, -3.0),
    ]
    field_lines = []
    ps_last = None
    for s in seeds:
        X, Y, Z, ps = trace_field_line(s, parmod, ut, rlim=60.0, r0=1.0)
        field_lines.append((X, Y, Z))
        ps_last = ps

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    ax_xy, ax_xz = axs

    # Left: X–Y (equatorial) magnetopause
    ax_xy.plot(xmp,  ymp, lw=2, label="Magnetopause (Shue 1998)")
    ax_xy.plot(xmp, -ymp, lw=2)
    ax_xy.plot(0, 0, "o", ms=6)
    ax_xy.plot(xbs, ybs, 'k--', lw=2, label="Bow shock")
    ax_xy.plot(xbs, -ybs, 'k--', lw=2)
    ax_xy.set_aspect("equal", adjustable="box")
    ax_xy.set_xlabel("X_GSM [R$_E$]")
    ax_xy.set_ylabel("Y_GSM [R$_E$]")
    ax_xy.set_title("Equatorial plane")
    ax_xy.grid(True)

    # Right: X–Z field lines
    for X, Y, Z in field_lines:
        ax_xz.plot(X, Z, lw=1.5)

    th = np.linspace(0, 2*np.pi, 361)
    ax_xz.plot(np.cos(th), np.sin(th), lw=2)  # Earth disk
    ax_xz.plot(xbs,  ybs, "k--", lw=2)
    ax_xz.plot(xbs, -ybs, "k--", lw=2)
    ax_xz.set_aspect("equal", adjustable="box")
    ax_xz.set_xlabel("X_GSM [R$_E$]")
    ax_xz.set_ylabel("Z_GSM [R$_E$]")
    ax_xz.set_title("Meridian plane (T96 + IGRF)")
    ax_xz.grid(True)

    ax_xy.set_xlim(-30, 30)
    ax_xy.set_ylim(-30, 30)
    ax_xz.set_xlim(-30, 30)
    ax_xz.set_ylim(-30, 30)

    fig.suptitle(
        f"{t_iso}Z | OMNI2: Pdyn={Pdyn:.2f} nPa, Dst={Dst:.0f} nT, "
        f"By={By:.2f} nT, Bz={Bz:.2f} nT, M_A={Ma:.2f}\n"
        f"MP: r0={r0mp:.2f} Re, α={alphamp:.2f} | "
        f"BS: r0={rbs:.2f} Re (Δ={delta_bs:.2f} Re) | "
        f"tilt(ps)={ps_last:.3f} rad",
        fontsize=10,
    )
    
    
    ax_xy.legend(loc="upper right", frameon=False)

    outpng = "figs/magnetosphere_T96_GSM_1999-02-21T00Z.png"
    fig.savefig(outpng, dpi=200)
    print(f"Wrote {outpng}")


if __name__ == "__main__":
    main()
