#!/usr/bin/env python3
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

import spacepy.time as spt
import spacepy.omni as om
import spacepy.toolbox as tb
import spacepy.coordinates as spc

from geopack import geopack as gp
from sscws.sscws import SscWs


# ---------------------------
# Utilities
# ---------------------------

RE_KM = 6371.2  # conventional Earth radius in km


def unix_seconds(t_utc: dt.datetime) -> float:
    if t_utc.tzinfo is None:
        t_utc = t_utc.replace(tzinfo=dt.timezone.utc)
    epoch = dt.datetime(1970, 1, 1, tzinfo=dt.timezone.utc)
    return (t_utc - epoch).total_seconds()


def _parse_times(time_field):
    """Accept list/array of datetime objects, ISO strings, or epoch seconds."""
    if isinstance(time_field, (dt.datetime, str, int, float)):
        time_list = [time_field]
    else:
        time_list = list(time_field)

    out = []
    for t in time_list:
        if isinstance(t, dt.datetime):
            if t.tzinfo is None:
                t = t.replace(tzinfo=dt.timezone.utc)
            out.append(t)
        elif isinstance(t, str):
            out.append(dt.datetime.fromisoformat(t.replace("Z", "+00:00")))
        elif isinstance(t, (int, float)):
            out.append(dt.datetime.fromtimestamp(t, tz=dt.timezone.utc))
        else:
            raise TypeError(f"Unrecognized time element type: {type(t)} value={t}")
    return out


def get_omni2_params(t_iso: str):
    """
    Use SpacePy OMNI2 hourly.
    Your install uses:
      Flow_pressure (nPa), Dst_index (nT), By_GSM (nT), Bz_GSM (nT), Alfven_mach_number
    """
    ticks = spt.Ticktock([t_iso], "ISO")
    try:
        omni = om.get_omni(ticks, dbase="OMNI2hourly")
    except Exception:
        tb.update(omni2=True)
        omni = om.get_omni(ticks, dbase="OMNI2hourly")

    Pdyn = float(omni["Flow_pressure"][0])
    Dst  = float(omni["Dst_index"][0])
    By   = float(omni["By_GSM"][0])
    Bz   = float(omni["Bz_GSM"][0])
    Ma   = float(omni["Alfven_mach_number"][0])
    return Pdyn, Dst, By, Bz, Ma


def t96_parmod(Pdyn, Dst, By, Bz):
    par = np.zeros(10, dtype=float)
    par[0] = Pdyn
    par[1] = Dst
    par[2] = By
    par[3] = Bz
    return par


# ---------------------------
# Magnetopause & bow shock
# ---------------------------

def shue98_magnetopause_rtheta(Pdyn, Bz, ntheta=361):
    """
    Shue et al. (1998):
      r(θ) = r0 * (2/(1+cosθ))^α
      r0 = (10.22 + 1.29*tanh(0.184*(Bz + 8.14))) * Pdyn^(-1/6.6)
      α  = (0.58 - 0.007*Bz) * (1 + 0.024*ln(Pdyn))
    """
    if not np.isfinite(Pdyn) or Pdyn <= 0:
        raise ValueError(f"Nonpositive/invalid Pdyn={Pdyn}")
    if not np.isfinite(Bz):
        raise ValueError(f"Invalid Bz={Bz}")

    r0 = (10.22 + 1.29 * np.tanh(0.184 * (Bz + 8.14))) * (Pdyn ** (-1.0 / 6.6))
    alpha = (0.58 - 0.007 * Bz) * (1.0 + 0.024 * np.log(Pdyn))

    eps = 1e-6
    theta = np.linspace(0, np.pi - eps, ntheta)
    r = r0 * (2.0 / (1.0 + np.cos(theta))) ** alpha

    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, r0, alpha


def bowshock_rtheta_from_Ma(r0_mp, Ma, ntheta=361, e=0.9, gamma=5/3,
                        delta_min=1.0, delta_max=10.0):
    """
    Mach-dependent bow shock nose distance:
      delta = Rbs0 - Rmp0 ~ r0_mp * ((γ-1)Ma^2 + 2)/((γ+1)(Ma^2 - 1)), clamped.
    Then conic shape scaled so r(0)=Rbs0.
    """
    if (not np.isfinite(Ma)) or (Ma <= 1.0):
        delta = 2.5
    else:
        M2 = Ma * Ma
        frac = ((gamma - 1.0) * M2 + 2.0) / ((gamma + 1.0) * (M2 - 1.0))
        delta = float(np.clip(r0_mp * frac, delta_min, delta_max))

    Rbs0 = r0_mp + delta

    eps = 1e-6
    theta = np.linspace(0.0, np.pi - eps, ntheta)
    r = Rbs0 * (1.0 + e) / (1.0 + e * np.cos(theta))

    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, Rbs0, delta


# ---------------------------
# Field line tracing
# ---------------------------

def trace_field_line(seed_xyz_gsm, parmod, ut_seconds, rlim=60.0, r0=1.0):
    ps = gp.recalc(ut_seconds)

    x0, y0, z0 = seed_xyz_gsm
    xf1, yf1, zf1, xx1, yy1, zz1 = gp.trace(x0, y0, z0, -1, rlim, r0, parmod, "t96", "igrf")
    xf2, yf2, zf2, xx2, yy2, zz2 = gp.trace(x0, y0, z0,  1, rlim, r0, parmod, "t96", "igrf")

    X = np.concatenate([np.asarray(xx1)[::-1], np.asarray(xx2)[1:]])
    Y = np.concatenate([np.asarray(yy1)[::-1], np.asarray(yy2)[1:]])
    Z = np.concatenate([np.asarray(zz1)[::-1], np.asarray(zz2)[1:]])
    return X, Y, Z, ps


# ---------------------------
# Wind ephemeris (SSC -> GSE -> GSM)
# ---------------------------

def _extract_xyz_any(sat_data):
    coords = sat_data.get("Coordinates")

    if isinstance(coords, dict):
        entries = [coords]
    else:
        entries = list(coords)  # works for numpy arrays too

    for c in entries:
        if isinstance(c, dict) and all(k in c for k in ("X", "Y", "Z")):
            x = np.asarray(c["X"], float)
            y = np.asarray(c["Y"], float)
            z = np.asarray(c["Z"], float)
            units = c.get("Units") or c.get("units")
            cs = c.get("CoordinateSystem")
            coord_sys = cs.value if hasattr(cs, "value") else (str(cs) if cs is not None else None)
            return x, y, z, units, coord_sys

    raise KeyError(f"Could not find X/Y/Z in Coordinates; type={type(coords)}")


def get_wind_gsm_from_ssc(start_iso, stop_iso, observatory_id="wind"):
    ssc = SscWs()
    result = ssc.get_locations([observatory_id], [start_iso, stop_iso])

    if "Data" not in result or not result.get("Data"):
        msg = "\n".join([f"{k}: {result.get(k)}" for k in ["HttpStatus", "ErrorMessage", "ErrorDescription"]])
        raise RuntimeError(f"SSCWeb returned no Data.\n{msg}")

    sat_data = result["Data"][0]
    times = _parse_times(sat_data["Time"])
    x, y, z, units, coord_sys = _extract_xyz_any(sat_data)

    # Convert to Re if necessary / infer
    if units is None:
        med = np.nanmedian(np.abs(x))
        units = "km" if med > 1000 else "Re"

    if units.lower() in ["km", "kilometer", "kilometers"]:
        x = x / RE_KM
        y = y / RE_KM
        z = z / RE_KM

    xyz_gse = np.vstack([x, y, z]).T

    # Convert GSE -> GSM with SpacePy
    tt = spt.Ticktock(times, "UTC")
    c_gse = spc.Coords(xyz_gse, "GSE", "car", ticks=tt)
    c_gsm = c_gse.convert("GSM", "car")
    xyz_gsm = np.asarray(c_gsm.data)

    return times, xyz_gsm


# ---------------------------
# Main
# ---------------------------

def main():
    t_iso = "1999-02-21T00:00:00"
    t0 = dt.datetime(1999, 2, 21, 0, 0, 0, tzinfo=dt.timezone.utc)

    # OMNI2 drivers
    Pdyn, Dst, By, Bz, Ma = get_omni2_params(t_iso)
    parmod = t96_parmod(Pdyn, Dst, By, Bz)
    ut = unix_seconds(t0)

    # MP + BS
    xmp, ymp, r0mp, alphamp = shue98_magnetopause_rtheta(Pdyn, Bz, ntheta=361)
    xbs, ybs, rbs0, delta_bs = bowshock_rtheta_from_Ma(r0mp, Ma, ntheta=361)

    r_seeds = [2.5]
    theta_max = np.deg2rad(85)
    
    thetas = np.linspace(-theta_max, theta_max, 11)

    # symmetric about z=0
    seeds = [
        (r*np.cos(th), 0.0,  r*np.sin(th))
        for r in r_seeds
        for th in thetas
    ]

    tail_max = np.deg2rad(35)
    tail = np.linspace(0.0, tail_max, 9)

    seeds += [
        (r*np.cos(np.pi - th), 0.0,  sgn*r*np.sin(np.pi - th))
        for r in [3.5, 5.0]         # a bit larger helps in the tail
        for th in tail
        for sgn in (+1.0, -1.0)
    ]

    # include the equator explicitly (optional; avoids duplicate z=0 from the mirror)
    seeds += [(r, 0.0, 0.0) for r in r_seeds]
    
    # Field lines (meridian plane seeds)
    #x_vals = [3.0, 2.0, -2.0, -3.0]
    #z_vals = [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]

    #seeds = [(x, 0.0, z) for x in x_vals for z in z_vals]
    field_lines = []
    ps_last = None
    for s in seeds:
        X, Y, Z, ps = trace_field_line(s, parmod, ut, rlim=60.0, r0=1.0)
        field_lines.append((X, Y, Z))
        ps_last = ps

    # Wind ephemeris (12h window)
    start = "1999-02-20T12:00:00Z"
    stop  = "1999-02-21T12:00:00Z"
    times_w, xyz_w = get_wind_gsm_from_ssc(start, stop, observatory_id="wind")

    # Find closest Wind point to snapshot
    i0 = int(np.argmin([abs((ti - t0).total_seconds()) for ti in times_w]))

    # Plot
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    ax_xy, ax_xz = axs

    # Equatorial: MP + BS
    ax_xy.plot(xmp,  ymp, "-.", lw=2, color='#32a852', label="Magnetopause (Shue 1998)")
    ax_xy.plot(xmp, -ymp, "-.", lw=2, color='#32a852')
    ax_xy.plot(xbs,  ybs, "--", lw=2, color='#326ba8', label="Bow shock ($M_A$)")
    ax_xy.plot(xbs, -ybs, "--", lw=2, color='#326ba8')

    # Wind (XY)
    ax_xy.plot(xyz_w[:, 0], xyz_w[:, 1], lw=2, color='#c01d14', label="Wind")
    ax_xy.plot(xyz_w[i0, 0], xyz_w[i0, 1], "o", color='#fc0394')

    ax_xy.plot(0, 0, "o", ms=6)
    ax_xy.set_aspect("equal", adjustable="box")
    ax_xy.set_xlabel(r"$X_{GSM} [R_E]$")
    ax_xy.set_ylabel(r"$Y_{GSM} [R_E]$")
    ax_xy.set_title("Equatorial plane")
    ax_xy.grid(True)
    
    # Meridian: field lines + BS + Wind
    for X, Y, Z in field_lines:
        ax_xz.plot(X, Z, lw=0.75)

    ax_xz.plot(xyz_w[:, 0], xyz_w[:, 2], lw=2, color='#c01d14')
    ax_xz.plot(xyz_w[i0, 0], xyz_w[i0, 2], "o", color='#fc0394')

    ax_xz.plot(xmp,  ymp, "-.", lw=2, color='#32a852')
    ax_xz.plot(xmp, -ymp, "-.", lw=2, color='#32a852')
    ax_xz.plot(xbs,  ybs, "--", lw=2, color='#326ba8')
    ax_xz.plot(xbs, -ybs, "--", lw=2, color='#326ba8')
    ax_xz.grid(True)

    th = np.linspace(0, 2*np.pi, 361)
    ax_xz.plot(np.cos(th), np.sin(th), lw=2)
    ax_xz.set_aspect("equal", adjustable="box")
    ax_xz.set_xlabel(r"$X_{GSM} [R_E]$")
    ax_xz.set_ylabel(r"$Z_{GSM} [R_E]$")
    ax_xz.set_title("Meridian plane (T96 + IGRF)")

    # Limits (your Wind is near-Earth, so these are fine)
    ax_xy.set_xlim(-40, 40)
    ax_xy.set_ylim(-40, 40)
    ax_xz.set_xlim(-40, 40)
    ax_xz.set_ylim(-40, 40)

    fig.text(
        0.5, 0.02,
        "Boundaries: axisymmetric empirical models (Shue 1998 MP + conic bow shock)",
        ha="center",
        fontsize=9,
        style="italic",
        alpha=0.8,
    )
    
    # Title (two lines)
    fig.suptitle(
        f"{t_iso}Z | OMNI2: Pdyn={Pdyn:.2f} nPa, Dst={Dst:.0f} nT, By={By:.2f} nT, Bz={Bz:.2f} nT, M_A={Ma:.2f}\n"
        f"MP: r0={r0mp:.2f} Re, α={alphamp:.2f} | BS: r0={rbs0:.2f} Re (Δ={delta_bs:.2f} Re) | tilt(ps)={ps_last:.3f} rad",
        fontsize=10,
    )

    ax_xy.legend(loc="upper left", frameon=False)

    outpng = "figs/Simple-Magnetosphere_wind_1999-02-21T00Z.png"
    fig.savefig(outpng, dpi=200)
    print(f"Wrote {outpng}")


if __name__ == "__main__":
    main()
