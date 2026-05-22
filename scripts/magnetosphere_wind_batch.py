#!/usr/bin/env python3
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import os

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


def parse_iso_utc(s: str) -> dt.datetime:
    """
    Robust time parser for collaborator formats.

    Accepts:
      1995-08-19 18:28:55.0
      1995-08-19T18:28:55.0
      1995-08-19 18:28:55
      1995-08-19T18:28:55
      1995-08-19T18:28:55Z

    Returns timezone-aware UTC datetime.
    """
    s = s.strip()

    # remove trailing Z
    if s.endswith("Z"):
        s = s[:-1]

    # allow space separator
    s = s.replace(" ", "T")

    # try fast path first
    try:
        t = dt.datetime.fromisoformat(s)
    except ValueError:
        # fallback formats (handles fractional seconds)
        for fmt in (
            "%Y-%m-%dT%H:%M:%S.%f",
            "%Y-%m-%dT%H:%M:%S",
        ):
            try:
                t = dt.datetime.strptime(s, fmt)
                break
            except ValueError:
                continue
        else:
            raise ValueError(f"Could not parse time string: {s}")

    # ensure UTC
    if t.tzinfo is None:
        t = t.replace(tzinfo=dt.timezone.utc)

    return t


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


# ---------------------------
# OMNI2 drivers
# ---------------------------

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


def convert_xyz(xyz, times, from_sys, to_sys):
    """
    xyz: (N,3) array in Re
    times: list of datetime (UTC)
    from_sys/to_sys: 'GSM' or 'GSE'
    """
    if from_sys.upper() == to_sys.upper():
        return np.asarray(xyz)

    tt = spt.Ticktock(times, "UTC")
    c = spc.Coords(np.asarray(xyz), from_sys.upper(), "car", ticks=tt)
    c2 = c.convert(to_sys.upper(), "car")
    return np.asarray(c2.data)

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

def plane_seeds_in_frame(frame, plane="xz"):
    """
    Return seeds defined IN the requested frame:
      - plane='xy': z=0 (equatorial)
      - plane='xz': y=0 (meridian)
    Seeds returned as (N,3) array IN THAT frame (Re).
    """
    frame = frame.upper()
    if plane == "xz":
        r_seeds = [2.5]
        theta_max = np.deg2rad(85)
        thetas = np.linspace(-theta_max, theta_max, 11)
        seeds = [(r*np.cos(th), 0.0, r*np.sin(th)) for r in r_seeds for th in thetas]

        tail_max = np.deg2rad(35)
        tail = np.linspace(0.0, tail_max, 9)
        seeds += [
            (r*np.cos(np.pi - th), 0.0, sgn*r*np.sin(np.pi - th))
            for r in [3.5, 5.0]
            for th in tail
            for sgn in (+1.0, -1.0)
        ]
        seeds += [(r, 0.0, 0.0) for r in r_seeds]

    elif plane == "xy":
        r_seeds = [2.5, 4.5]
        phi = np.linspace(0, 2*np.pi, 16, endpoint=False)
        seeds = [(r*np.cos(p), r*np.sin(p), 0.0) for r in r_seeds for p in phi]

    else:
        raise ValueError("plane must be 'xy' or 'xz'")

    return np.asarray(seeds, float)

def trace_field_lines_for_plot(frame, plane, parmod, ut_seconds, t_snap_dt, rlim=60.0, r0=1.0):
    """
    Returns list of (X, Y, Z) arrays IN 'frame' coordinates.
    """
    # seeds defined in requested frame
    seeds_frame = plane_seeds_in_frame(frame, plane=plane)

    # convert seeds to GSM for geopack.trace
    seeds_gsm = convert_xyz(seeds_frame, [t_snap_dt]*len(seeds_frame), frame, "GSM")

    out = []
    ps_last = None
    for s_gsm in seeds_gsm:
        Xg, Yg, Zg, ps = trace_field_line(tuple(s_gsm), parmod, ut_seconds, rlim=rlim, r0=r0)

        # convert traced polyline back to requested frame
        xyz_gsm = np.vstack([Xg, Yg, Zg]).T
        xyz_frame = convert_xyz(xyz_gsm, [t_snap_dt]*len(Xg), "GSM", frame)

        out.append((xyz_frame[:, 0], xyz_frame[:, 1], xyz_frame[:, 2]))
        ps_last = ps

    return out, ps_last

def convert_curve_xy_to_frame(x, y, frame, t_snap_dt, plane="xy"):
    """
    Convert a curve from GSM to 'frame' (GSE/GSM) keeping it in the chosen plane
    by embedding into 3D first.
      plane='xy' -> z=0
      plane='xz' -> y=0, use y-array as z-array
    Returns x_plot, t_plot where t is y for xy plane or z for xz plane.
    """
    frame = frame.upper()
    if plane == "xy":
        xyz_gsm = np.vstack([x, y, np.zeros_like(x)]).T
        xyz_f = convert_xyz(xyz_gsm, [t_snap_dt]*len(x), "GSM", frame)
        return xyz_f[:, 0], xyz_f[:, 1]
    elif plane == "xz":
        xyz_gsm = np.vstack([x, np.zeros_like(x), y]).T
        xyz_f = convert_xyz(xyz_gsm, [t_snap_dt]*len(x), "GSM", frame)
        return xyz_f[:, 0], xyz_f[:, 2]
    else:
        raise ValueError("plane must be 'xy' or 'xz'")

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


def get_wind_gse_and_gsm_from_ssc(start_iso, stop_iso, observatory_id="wind"):
    ssc = SscWs()
    result = ssc.get_locations([observatory_id], [start_iso, stop_iso])

    data = result.get("Data", None)
    empty = (
        data is None or
        (hasattr(data, "size") and data.size == 0) or
        (hasattr(data, "__len__") and len(data) == 0)
    )
    if empty:
        msg = "\n".join([f"{k}: {result.get(k)}" for k in ["HttpStatus", "ErrorMessage", "ErrorDescription"]])
        raise RuntimeError(f"SSCWeb returned no Data.\n{msg}\nKeys: {list(result.keys())}")

    sat_data = result["Data"][0]
    times = _parse_times(sat_data["Time"])
    x, y, z, units, coord_sys = _extract_xyz_any(sat_data)

    if units is None:
        med = np.nanmedian(np.abs(x))
        units = "km" if med > 1000 else "Re"
    if units.lower() in ["km", "kilometer", "kilometers"]:
        x = x / RE_KM; y = y / RE_KM; z = z / RE_KM

    xyz_gse = np.vstack([x, y, z]).T  # SSC delivers GSE in your current workflow
    xyz_gsm = convert_xyz(xyz_gse, times, "GSE", "GSM")

    return times, xyz_gse, xyz_gsm


# ---------------------------
# Seeds
# ---------------------------

def default_meridian_seeds():
    """
    Angular seeds in XZ (y=0), plus a few tail seeds.
    Returns list of (x,y,z) in Re.
    """
    r_seeds = [2.5]
    theta_max = np.deg2rad(85)
    thetas = np.linspace(-theta_max, theta_max, 11)

    seeds = [(r*np.cos(th), 0.0, r*np.sin(th)) for r in r_seeds for th in thetas]

    tail_max = np.deg2rad(35)
    tail = np.linspace(0.0, tail_max, 9)
    seeds += [
        (r*np.cos(np.pi - th), 0.0, sgn*r*np.sin(np.pi - th))
        for r in [3.5, 5.0]
        for th in tail
        for sgn in (+1.0, -1.0)
    ]

    seeds += [(r, 0.0, 0.0) for r in r_seeds]
    return seeds

def find_valid_omni_within_interval(start_dt, stop_dt, prefer="midpoint", step_hours=1, max_steps=2000):
    """
    Try to find a snapshot time within [start_dt, stop_dt] where OMNI2 drivers are valid.
    Steps forward by step_hours until stop_dt is exceeded.

    Returns:
      t_snap (datetime UTC), t_iso (str), Pdyn, Dst, By, Bz, Ma

    Raises:
      RuntimeError if no valid OMNI point found in the interval.
    """
    if start_dt.tzinfo is None:
        start_dt = start_dt.replace(tzinfo=dt.timezone.utc)
    if stop_dt.tzinfo is None:
        stop_dt = stop_dt.replace(tzinfo=dt.timezone.utc)

    if prefer == "midpoint":
        t0 = start_dt + (stop_dt - start_dt) / 2
    elif prefer == "start":
        t0 = start_dt
    else:
        raise ValueError("prefer must be 'midpoint' or 'start'")

    t = t0
    n = 0

    while t <= stop_dt and n < max_steps:
        t_iso = t.strftime("%Y-%m-%dT%H:%M:%S")
        Pdyn, Dst, By, Bz, Ma = get_omni2_params(t_iso)

        ok = (
            np.isfinite(Pdyn) and (Pdyn > 0) and
            np.isfinite(Dst) and
            np.isfinite(By) and
            np.isfinite(Bz) and
            np.isfinite(Ma)
        )

        print(
            f"[OMNI2 try] idx? {t_iso} | "
            f"Pdyn={Pdyn}, Dst={Dst}, By={By}, Bz={Bz}, Ma={Ma} | ok={ok}"
        )

        if ok:
            return t, t_iso, Pdyn, Dst, By, Bz, Ma

        t = t + dt.timedelta(hours=step_hours)
        n += 1

    raise RuntimeError("No valid OMNI2 drivers found within interval")

# ---------------------------
# Plot one interval
# ---------------------------
def plot_interval(idx: int, start_dt: dt.datetime, stop_dt: dt.datetime,
                  outdir="figs", observatory_id="wind"):

    print(f"Interval {idx:02d} | start={start_dt} stop={stop_dt}")
    os.makedirs(outdir, exist_ok=True)

    # -------------------------
    # 1) Choose a valid OMNI snapshot (once)
    # -------------------------
    try:
        t_snap, t_iso, Pdyn, Dst, By, Bz, Ma = find_valid_omni_within_interval(
            start_dt, stop_dt, prefer="midpoint", step_hours=1
        )
    except RuntimeError as e:
        print(f"⚠ Skipping interval {idx:02d}: {e}")
        return

    parmod = t96_parmod(Pdyn, Dst, By, Bz)
    ut = unix_seconds(t_snap)

    # MP + BS curves computed once (treat as defined in GSM, then transform for plotting)
    xmp_gsm, ymp_gsm, r0mp, alphamp = shue98_magnetopause_rtheta(Pdyn, Bz, ntheta=361)
    xbs_gsm, ybs_gsm, rbs0, delta_bs = bowshock_rtheta_from_Ma(r0mp, Ma, ntheta=361)

    # -------------------------
    # 2) Wind ephemeris once (returns both frames)
    # -------------------------
    start_iso = start_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    stop_iso  = stop_dt.strftime("%Y-%m-%dT%H:%M:%SZ")

    try:
        times_w, xyz_w_gse, xyz_w_gsm = get_wind_gse_and_gsm_from_ssc(
            start_iso, stop_iso, observatory_id=observatory_id
        )
    except Exception as e:
        print(f"⚠ Skipping interval {idx:02d}: SSC error: {e}")
        return

    # closest Wind sample to snapshot
    i0 = int(np.argmin([abs((ti - t_snap).total_seconds()) for ti in times_w]))

    # -------------------------
    # 3) Plot in each frame (GSM and GSE) without re-downloading anything
    # -------------------------
    for frame in ("GSM", "GSE"):

        # choose Wind coords for this frame
        xyz_w = xyz_w_gsm if frame == "GSM" else xyz_w_gse

        # --- field lines for this frame
        # XY plane (equatorial) and XZ plane (meridian) IN THIS FRAME
        try:
            fl_xy, ps_last_xy = trace_field_lines_for_plot(frame, "xy", parmod, ut, t_snap, rlim=60.0, r0=1.0)
            fl_xz, ps_last_xz = trace_field_lines_for_plot(frame, "xz", parmod, ut, t_snap, rlim=60.0, r0=1.0)
            ps_last = ps_last_xz if ps_last_xz is not None else ps_last_xy
        except Exception as e:
            print(f"⚠ Skipping interval {idx:02d} ({frame}): trace error: {e}")
            continue

        # --- transform boundary curves into this frame & plane
        xmp_xy, ymp_xy = convert_curve_xy_to_frame(xmp_gsm, ymp_gsm, frame, t_snap, plane="xy")
        xbs_xy, ybs_xy = convert_curve_xy_to_frame(xbs_gsm, ybs_gsm, frame, t_snap, plane="xy")

        xmp_xz, zmp_xz = convert_curve_xy_to_frame(xmp_gsm, ymp_gsm, frame, t_snap, plane="xz")
        xbs_xz, zbs_xz = convert_curve_xy_to_frame(xbs_gsm, ybs_gsm, frame, t_snap, plane="xz")

        # -------------------------
        # 4) Render figure
        # -------------------------
        fig, axs = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        ax_xy, ax_xz = axs

        # ---- XY panel (equatorial)
        ax_xy.plot(xmp_xy,  ymp_xy, "-.", lw=2, color="#32a852", label="Magnetopause (Shue 1998)")
        ax_xy.plot(xmp_xy, -ymp_xy, "-.", lw=2, color="#32a852")
        ax_xy.plot(xbs_xy,  ybs_xy, "--", lw=2, color="#326ba8", label="Bow shock ($M_A$)")
        ax_xy.plot(xbs_xy, -ybs_xy, "--", lw=2, color="#326ba8")

        # Wind trajectory (XY)
        ax_xy.plot(xyz_w[:, 0], xyz_w[:, 1], lw=2, color="#c01d14", label="Wind (interval)")
        ax_xy.plot(xyz_w[i0, 0], xyz_w[i0, 1], "o", color="#fc0394")

        # ---- XY panel (meridian)
        # field lines (projected into XY for this frame)
        for X, Y, Z in fl_xy:
            ax_xy.plot(X, Z, lw=0.75)
        
        ax_xy.plot(0, 0, "o", ms=6)
        ax_xy.set_aspect("equal", adjustable="box")
        ax_xy.set_xlabel(rf"$X_{{{frame}}}\ [R_E]$")
        ax_xy.set_ylabel(rf"$Y_{{{frame}}}\ [R_E]$")
        ax_xy.set_title(f"Equatorial plane ({frame})")
        ax_xy.grid(True)

        # ---- XZ panel (meridian)
        # field lines (projected into XZ for this frame)
        for X, Y, Z in fl_xz:
            ax_xz.plot(X, Z, lw=0.75)

        # Wind trajectory (XZ)
        ax_xz.plot(xyz_w[:, 0], xyz_w[:, 2], lw=2, color="#c01d14")
        ax_xz.plot(xyz_w[i0, 0], xyz_w[i0, 2], "o", color="#fc0394")

        # boundaries in XZ plane
        ax_xz.plot(xmp_xz,  zmp_xz, "-.", lw=2, color="#32a852")
        ax_xz.plot(xmp_xz, -zmp_xz, "-.", lw=2, color="#32a852")
        ax_xz.plot(xbs_xz,  zbs_xz, "--", lw=2, color="#326ba8")
        ax_xz.plot(xbs_xz, -zbs_xz, "--", lw=2, color="#326ba8")

        # Earth
        th = np.linspace(0, 2*np.pi, 361)
        ax_xz.plot(np.cos(th), np.sin(th), lw=2)

        ax_xz.set_aspect("equal", adjustable="box")
        ax_xz.set_xlabel(rf"$X_{{{frame}}}\ [R_E]$")
        ax_xz.set_ylabel(rf"$Z_{{{frame}}}\ [R_E]$")
        ax_xz.set_title(f"Meridian plane ({frame})")
        ax_xz.grid(True)

        # Limits
        ax_xy.set_xlim(-70, 40)
        ax_xy.set_ylim(-40, 40)
        ax_xz.set_xlim(-70, 40)
        ax_xz.set_ylim(-40, 40)

        fig.text(
            0.5, 0.02,
            "Boundaries: axisymmetric empirical models.\n"
            "axisymmetric empirical models (Shue 1998 MP + conic bow shock).\n"
            "Field: T96 + IGRF traced in GSM, plotted in selected frame.",
            ha="center",
            fontsize=8.5,
            style="italic",
            alpha=0.85,
        )        
        fig.suptitle(
            f"{t_snap.strftime('%Y-%m-%dT%H:%M:%SZ')} | Wind {start_iso} → {stop_iso} | Frame={frame}\n"
            f"OMNI2@used: Pdyn={Pdyn:.2f} nPa, Dst={Dst:.0f} nT, By={By:.2f} nT, Bz={Bz:.2f} nT, M_A={Ma:.2f} | "
            f"MP r0={r0mp:.2f} Re α={alphamp:.2f} | BS r0={rbs0:.2f} Re (Δ={delta_bs:.2f} Re) | tilt(ps)={ps_last:.3f} rad",
            fontsize=10,
        )

        ax_xy.legend(loc="upper right", frameon=False)

        outpng = os.path.join(outdir, f"Crossing_{idx:02d}_{frame}.png")
        fig.savefig(outpng, dpi=200)
        plt.close(fig)
        print(f"Wrote {outpng}")


# ---------------------------
# Parse wind-bowshock.txt and run
# ---------------------------

def parse_intervals_file(path):
    """
    Expected format (whitespace-delimited):
      idx  yyyy-mm-dd  hh:mm:ss.s  yyyy-mm-dd  hh:mm:ss.s
    Example:
      2 1995-08-19 18:28:55.0 1995-08-20 23:46:15.5
    """
    intervals = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            idx = int(parts[0])
            start = parse_iso_utc(parts[1] + " " + parts[2])
            stop  = parse_iso_utc(parts[3] + " " + parts[4])
            intervals.append((idx, start, stop))
    return intervals


def main():
    #intervals_path = "data/wind-bowshock.txt"
    intervals_path = "data/wind-bowshock-short.txt"
    intervals = parse_intervals_file(intervals_path)

    print(f"Found {len(intervals)} intervals in {intervals_path}")
    for idx, start_dt, stop_dt in intervals:
        plot_interval(idx, start_dt, stop_dt, outdir="figs", observatory_id="wind")


if __name__ == "__main__":
    main()
