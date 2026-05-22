
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


def parse_time_utc(t_iso):
    """
    Parse ISO time string and return:
      - timezone-aware datetime (UTC)
      - unix seconds

    Accepts:
      '1999-02-21T00:00:00'
      '1999-02-21T00:00’
      '1999-02-21T00:00:00Z'
    """

    # normalize trailing Z if present
    if t_iso.endswith("Z"):
        t_iso = t_iso[:-1]

    t = dt.datetime.fromisoformat(t_iso)

    # force UTC if not specified
    if t.tzinfo is None:
        t = t.replace(tzinfo=dt.timezone.utc)

    ut = t.timestamp()

    return t, ut

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

def _trace_end_radius(x0, y0, z0, direction, rlim, r0, parmod):
    """
    Trace from (x0,y0,z0) in 'direction' (+1 or -1) and return radius of final point.
    geopack.trace returns: xf, yf, zf, xx, yy, zz
    """
    xf, yf, zf, xx, yy, zz = gp.trace(x0, y0, z0, direction, rlim, r0, parmod, "t96", "igrf")

    # Be defensive: if trace returns NaNs, treat as open
    if not np.isfinite(xf) or not np.isfinite(yf) or not np.isfinite(zf):
        return np.inf

    return float(np.sqrt(xf*xf + yf*yf + zf*zf))

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


def is_closed_fieldline(x0, y0, z0, ut_seconds, parmod, rlim=60.0, r0=1.0, rtol=0.05):
    """
    A field line is 'closed' if both tracing directions return to Earth (radius <= r0+rtol).
    Otherwise it's open (one/both ends reach rlim or fail).

    Assumes gp.recalc(ut_seconds) has been called for this time already, or calls it here.
    """
    gp.recalc(ut_seconds)

    r_end_1 = _trace_end_radius(x0, y0, z0, +1, rlim, r0, parmod)
    r_end_2 = _trace_end_radius(x0, y0, z0, -1, rlim, r0, parmod)

    # "hits Earth" test
    hits_earth_1 = (r_end_1 <= (r0 + rtol))
    hits_earth_2 = (r_end_2 <= (r0 + rtol))

    return bool(hits_earth_1 and hits_earth_2)


def find_boundary_on_ray(theta, plane, ut_seconds, parmod,
                         r_min=1.2, r_max=60.0, dr=0.5,
                         rlim=60.0, r0=1.0, rtol=0.05,
                         bisect_tol=0.05, max_bisect=25,
                         last_r_guess=None):
    """
    Find magnetopause boundary radius r(theta) along a ray in a given plane.

    plane: "xy" (z=0) or "xz" (y=0)
    theta: polar angle in that plane, 0 at +X, pi at -X

    Strategy:
      - start at small r (closed)
      - march outward until open (or hit r_max)
      - bisection to refine transition

    Returns:
      r_boundary (float or np.nan), and a boolean indicating if a transition was found.
    """
    # Ray unit vector in chosen plane
    ct, st = np.cos(theta), np.sin(theta)
    if plane == "xy":        
        def point(r):
            return (r*ct, r*st, 0.0)
    elif plane == "xz":
        def point(r):
            return (r*ct, 0.0, r*st)
    else:
        raise ValueError("plane must be 'xy' or 'xz'")

    # Use last_r_guess to speed things up (boundary varies smoothly with theta)
    if last_r_guess is not None and np.isfinite(last_r_guess):
        r_start = max(r_min, last_r_guess - 2.0)
    else:
        r_start = r_min

    # Must start from closed region; if it's already open at r_start, walk inward a bit
    r = r_start
    x, y, z = point(r)
    closed = is_closed_fieldline(x, y, z, ut_seconds, parmod, rlim=rlim, r0=r0, rtol=rtol)

    if not closed:
        # try stepping inward (but don't go below r_min)
        for _ in range(10):
            r_try = max(r_min, r - dr)
            if r_try == r:
                break
            r = r_try
            x, y, z = point(r)
            closed = is_closed_fieldline(x, y, z, ut_seconds, parmod, rlim=rlim, r0=r0, rtol=rtol)
            if closed:
                break

    if not closed:
        # even near Earth this ray is open in the model (can happen in far tail)
        return np.nan, False

    # march outward to find first open point
    r_lo = r
    r_hi = None

    r = r_lo
    while r < r_max:
        r_next = r + dr
        x, y, z = point(r_next)
        closed_next = is_closed_fieldline(x, y, z, ut_seconds, parmod, rlim=rlim, r0=r0, rtol=rtol)
        if not closed_next:
            r_hi = r_next
            break
        r = r_next
        r_lo = r

    if r_hi is None:
        # never turned open out to r_max
        return np.nan, False

    # bisection between r_lo (closed) and r_hi (open)
    a, b = r_lo, r_hi
    for _ in range(max_bisect):
        m = 0.5*(a + b)
        x, y, z = point(m)
        closed_m = is_closed_fieldline(x, y, z, ut_seconds, parmod, rlim=rlim, r0=r0, rtol=rtol)
        if closed_m:
            a = m
        else:
            b = m
        if (b - a) <= bisect_tol:
            break

    r_boundary = 0.5*(a + b)
    return float(r_boundary), True

def magnetopause_slice_from_last_closed(ut_seconds, parmod, plane,
                                        ntheta=181, r_min=1.2, r_max=60.0,
                                        dr=0.5, rlim=60.0,
                                        bisect_tol=0.05):
    """
    Returns x, transverse (y or z) arrays for the magnetopause slice in the chosen plane.
    plane: 'xy' or 'xz'
    """
    thetas = np.linspace(0.0, np.pi, ntheta)  # +X to -X
    r_prev = None

    x_list = []
    t_list = []  # transverse: y for xy plane, z for xz plane

    if plane == "xy":
        print("Evaluating XY Plane")
    elif plane == "xz":
        print("Evaluating XZ Plane")
    
    for th in thetas:
        print(f"theta={th}")
        r_b, ok = find_boundary_on_ray(
            th, plane, ut_seconds, parmod,
            r_min=r_min, r_max=r_max, dr=dr,
            rlim=rlim, bisect_tol=bisect_tol,
            last_r_guess=r_prev
        )
        if ok and np.isfinite(r_b):
            ct, st = np.cos(th), np.sin(th)
            x = r_b*ct
            trans = r_b*st
            x_list.append(x)
            t_list.append(trans)
            r_prev = r_b
        else:
            # skip points where no boundary found (often deep tail)
            r_prev = r_prev  # keep previous guess

    return np.array(x_list), np.array(t_list)


# ---------------------------
# Main
# ---------------------------

def main():
    t_iso = "1999-02-21T00:00:00"
    t0, ut = parse_time_utc(t_iso)
    
    # OMNI2 drivers
    Pdyn, Dst, By, Bz, Ma = get_omni2_params(t_iso)
    parmod = t96_parmod(Pdyn, Dst, By, Bz)

    # Self-consistent boundary from T96+IGRF last-closed method
    xmp_xy, ymp_xy = magnetopause_slice_from_last_closed(ut, parmod, plane="xy", ntheta=91, dr=0.6, bisect_tol=0.05)
    xmp_xz, zmp_xz = magnetopause_slice_from_last_closed(ut, parmod, plane="xz", ntheta=91, dr=0.6, bisect_tol=0.05)

    # Empirical (optional – for comparison)
    xmp, ymp, r0mp, alphamp = shue98_magnetopause_rtheta(Pdyn, Bz, ntheta=361)
    xbs, ybs, rbs0, delta_bs = bowshock_rtheta_from_Ma(r0mp, Ma, ntheta=361)
    
    # Field lines (meridian plane seeds)
    x_vals = [3.0, 2.0, -2.0, -3.0]
    z_vals = [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]

    seeds = [(x, 0.0, z) for x in x_vals for z in z_vals]
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
    ax_xy.plot(xmp_xy,  ymp_xy, "-.", lw=2, color='#32a852', label="MP (last closed, T96+IGRF)")
    ax_xy.plot(xmp_xy, -ymp_xy, "-.", lw=2, color='#32a852')
    ax_xy.plot(xmp,  ymp, "-.", lw=2, color='#a0a0a0', label="Magnetopause (Shue 1998)")
    ax_xy.plot(xmp, -ymp, "-.", lw=2, color='#a0a0a0')
    ax_xy.plot(xbs,  ybs, "--", lw=2, color='#a0a0a0', label="Bow shock ($M_A$)")
    ax_xy.plot(xbs, -ybs, "--", lw=2, color='#a0a0a0')

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

    ax_xz.plot(xmp_xz,  zmp_xz, "-.", lw=2, color='#32a852')
    ax_xz.plot(xmp_xz, -zmp_xz, "-.", lw=2, color='#32a852')
    ax_xz.plot(xmp,  ymp, "-.", lw=2, color='#a0a0a0')
    ax_xz.plot(xmp, -ymp, "-.", lw=2, color='#a0a0a0')
    ax_xz.plot(xbs,  ybs, "--", lw=2, color='#a0a0a0')
    ax_xz.plot(xbs, -ybs, "--", lw=2, color='#a0a0a0')
    ax_xz.grid(True)

    th = np.linspace(0, 2*np.pi, 361)
    ax_xz.plot(np.cos(th), np.sin(th), lw=2)
    ax_xz.set_aspect("equal", adjustable="box")
    ax_xz.set_xlabel(r"$X_{GSM} [R_E]$")
    ax_xz.set_ylabel(r"$Z_{GSM} [R_E]$")
    ax_xz.set_title("Meridian plane (T96 + IGRF)")

    # Limits (your Wind is near-Earth, so these are fine)
    ax_xy.set_xlim(-35, 35)
    ax_xy.set_ylim(-35, 35)
    ax_xz.set_xlim(-35, 35)
    ax_xz.set_ylim(-35, 35)

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
        f"MP: r0={r0mp:.2f} Re, α={alphamp:.2f} | tilt(ps)={ps_last:.3f} rad",
        fontsize=10,
    )

    ax_xy.legend(loc="upper left", frameon=False)

    outpng = "figs/MP-T96_wind_1999-02-21T00Z.png"
    fig.savefig(outpng, dpi=200)
    print(f"Wrote {outpng}")


if __name__ == "__main__":
    main()
