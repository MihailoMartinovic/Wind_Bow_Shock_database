#!/usr/bin/env python3
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator
from pathlib import Path
import os

from spacepy import pycdf
import spacepy.time as spt
import spacepy.omni as om
import spacepy.toolbox as tb
import spacepy.coordinates as spc

from geopack import geopack as gp

import spiceypy as sp
from spacepy.coordinates import Coords


def _ensure_utc(dt_in: dt.datetime) -> dt.datetime:
    if dt_in.tzinfo is None:
        return dt_in.replace(tzinfo=dt.timezone.utc)
    return dt_in.astimezone(dt.timezone.utc)


def load_spice_kernels(kernel_dir="data/spice",
                       lsk="naif0012.tls",
                       spk="de440s.bsp",
                       pck="pck00010.tpc"):
    """
    Load SPICE kernels needed for Moon position.
    Call once at startup.
    """
    kernel_dir = os.path.expanduser(kernel_dir)
    paths = [
        os.path.join(kernel_dir, lsk),
        os.path.join(kernel_dir, spk),
        os.path.join(kernel_dir, pck),
    ]
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        msg = "Missing SPICE kernel(s):\n" + "\n".join("  - " + m for m in missing)
        msg += f"\nExpected in: {os.path.abspath(kernel_dir)}"
        raise FileNotFoundError(msg)

    # Clear any previously loaded kernels (optional but nice in dev)
    sp.kclear()
    for p in paths:
        sp.furnsh(p)


def moon_pos_geocentric_j2000_km(t_utc: dt.datetime) -> np.ndarray:
    """
    Moon position relative to Earth center in J2000 (ECI2000), km.
    Returns np.array([x,y,z]) in km.
    """
    t_utc = _ensure_utc(t_utc)

    # SPICE wants ET seconds past J2000
    et = sp.str2et(t_utc.strftime("%Y-%m-%dT%H:%M:%S"))

    # Position of MOON relative to EARTH, in J2000, no aberration corrections
    # Returns (pos, lt) where pos is km
    pos_km, _ = sp.spkpos("MOON", et, "J2000", "NONE", "EARTH")
    return np.array(pos_km, dtype=float)


def convert_j2000_to_gse_gsm(pos_km: np.ndarray, t_utc: dt.datetime):
    """
    Convert a single J2000 position vector (km) to GSE and GSM using SpacePy.
    Returns (pos_gse_km, pos_gsm_km).
    """
    t_utc = _ensure_utc(t_utc)

    # SpacePy tick
    ticks = spt.Ticktock([t_utc], "UTC")

    # Coords expects shape (N,3)
    c = Coords(pos_km.reshape(1, 3), "ECI2000", "car")  # J2000-ish inertial
    c.ticks = ticks
    c.units = "km"

    c_gse = c.convert("GSE", "car")
    c_gsm = c.convert("GSM", "car")

    return np.array(c_gse.data[0], float), np.array(c_gsm.data[0], float)


def get_moon_position(t_utc: dt.datetime, frame="GSE", units="Re") -> np.ndarray:
    """
    Convenience wrapper.
    frame: "J2000", "GSE", or "GSM"
    units: "km" or "Re"
    """
    t_utc = _ensure_utc(t_utc)
    pos_j2000_km = moon_pos_geocentric_j2000_km(t_utc)

    if frame.upper() == "J2000":
        pos_km = pos_j2000_km
    else:
        pos_gse_km, pos_gsm_km = convert_j2000_to_gse_gsm(pos_j2000_km, t_utc)
        if frame.upper() == "GSE":
            pos_km = pos_gse_km
        elif frame.upper() == "GSM":
            pos_km = pos_gsm_km
        else:
            raise ValueError(f"Unknown frame={frame}. Use 'J2000','GSE','GSM'.")

    if units.lower() == "km":
        return pos_km
    elif units.lower() == "re":
        return pos_km / RE_KM
    else:
        raise ValueError("units must be 'km' or 'Re'")

# ---------------------------
# Utilities
# ---------------------------

RE_KM = 6371.2  # Earth radius in km

FS_color='#00aa00'
MP_color="#a0a0a0"
BS_color="#326ba8"

HS_COLORS = [
    "#56B4E9", "#E69F00", "#8E4D4D", "#888888",
    "#F0E442", "#D55E00", "#009E73", "#8D00B2", "#CC79A7"
]  # hub+8 nodes, matching hs-single-hour-baselines.py
HUB_COLOR = HS_COLORS[0]
NODE_COLORS = HS_COLORS[1:9]

def get_omni2_params(t_iso: str, t_drv: dt.datetime | None = None):
    """
    Robust OMNI2hourly driver getter.

    Returns:
      Pdyn, Dst,
      Bx_gse, By_gse, Bz_gse,
      Bx_gsm, By_gsm, Bz_gsm,
      Ma

    Notes:
    - SpacePy OMNI2hourly often lacks Bx_GSM. We reconstruct GSM from GSE when needed.
    - If reconstructing, we use t_drv for the coordinate transform time if provided;
      otherwise we use the same time represented by t_iso.
    """
    ticks = spt.Ticktock([t_iso], "ISO")
    try:
        omni = om.get_omni(ticks, dbase="OMNI2hourly")
    except Exception:
        tb.update(omni2=True)
        omni = om.get_omni(ticks, dbase="OMNI2hourly")

    # Core drivers
    Pdyn = float(omni["Flow_pressure"][0])
    Dst  = float(omni["Dst_index"][0])
    Ma   = float(omni["Alfven_mach_number"][0])

    # ---- GSE IMF (these are usually present)
    bx_gse = float(omni["Bx_GSE"][0]) if "Bx_GSE" in omni else np.nan
    by_gse = float(omni["By_GSE"][0]) if "By_GSE" in omni else np.nan
    bz_gse = float(omni["Bz_GSE"][0]) if "Bz_GSE" in omni else np.nan

    # ---- GSM IMF (By/Bz often present; Bx often missing)
    #by_gsm = float(omni["By_GSM"][0]) if "By_GSM" in omni else np.nan
    #bz_gsm = float(omni["Bz_GSM"][0]) if "Bz_GSM" in omni else np.nan
    #bx_gsm = float(omni["Bx_GSM"][0]) if "Bx_GSM" in omni else np.nan
    
    # Reconstruct GSM components if needed / desired
    #need_gsm = (not np.isfinite(bx_gsm)) or (not np.isfinite(by_gsm)) or (not np.isfinite(bz_gsm))
    have_gse = np.isfinite(bx_gse) and np.isfinite(by_gse) and np.isfinite(bz_gse)

    #if need_gsm and have_gse:
    if have_gse:
        # Use driver time for transform if provided
        if t_drv is None:
            # derive from t_iso
            t_drv = dt.datetime.fromisoformat(t_iso.replace("Z", "+00:00"))
        if t_drv.tzinfo is None:
            t_drv = t_drv.replace(tzinfo=dt.timezone.utc)

        v_gse = np.array([[bx_gse, by_gse, bz_gse]], dtype=float)
        v_gsm = convert_xyz(v_gse, [t_drv], "GSE", "GSM")[0]
        bx_gsm, by_gsm, bz_gsm = float(v_gsm[0]), float(v_gsm[1]), float(v_gsm[2])

    return Pdyn, Dst, bx_gse, by_gse, bz_gse, bx_gsm, by_gsm, bz_gsm, Ma

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

def bowshock_tangent_index_signed(xbs, ybs, bhat, side_mask=None):
    """
    Returns index i that maximizes dot(t_hat, bhat) (signed, not abs).
    If side_mask is provided, restrict search to that boolean mask.
    """
    xbs = np.asarray(xbs, float)
    ybs = np.asarray(ybs, float)

    dx = np.gradient(xbs)
    dy = np.gradient(ybs)
    t = np.vstack([dx, dy]).T
    nt = np.linalg.norm(t, axis=1)
    nt[nt == 0] = np.nan
    that = t / nt[:, None]

    score = that @ bhat  # signed cos(angle)

    if side_mask is not None:
        score = np.where(side_mask, score, -np.inf)

    i = int(np.nanargmax(score))
    return i, score

def bowshock_tangent_indices_to_imf(xbs, ybs, bx, by, ysplit=0.0):
    """
    Find bow shock point(s) where the curve tangent is most parallel to IMF in XY.

    Returns:
      i_pos, i_neg, bhat, that, score
    where score ~ |that·bhat| (closer to 1 is more tangent).
    """
    xbs = np.asarray(xbs, float)
    ybs = np.asarray(ybs, float)

    b = np.array([bx, by], float)
    nb = np.linalg.norm(b)
    if not np.isfinite(nb) or nb == 0:
        raise ValueError(f"Bad IMF in-plane vector: bx={bx}, by={by}")
    bhat = b / nb

    # tangent along the curve
    dx = np.gradient(xbs)
    dy = np.gradient(ybs)
    t = np.vstack([dx, dy]).T
    nt = np.linalg.norm(t, axis=1)
    nt[nt == 0] = np.nan
    that = t / nt[:, None]

    score = np.abs(that @ bhat)  # |cos(angle)|

    # pick best match separately on +y and -y sides
    pos = np.where(ybs >= ysplit)[0]
    neg = np.where(ybs <= ysplit)[0]

    i_pos = pos[np.nanargmax(score[pos])] if pos.size else None
    i_neg = neg[np.nanargmax(score[neg])] if neg.size else None

    return i_pos, i_neg, bhat, that, score

def local_tangent_normal(x, y, i):
    """
    Estimate tangent and outward normal at index i on a 2D curve.
    """
    dx = np.gradient(x)
    dy = np.gradient(y)
    t = np.array([dx[i], dy[i]], float)
    nt = np.linalg.norm(t)
    if nt == 0 or not np.isfinite(nt):
        return None, None
    that = t / nt

    # normal is tangent rotated by +90 deg (one of the two normals)
    nhat = np.array([-that[1], that[0]], float)

    # choose the normal that points roughly sunward (toward +x) OR outward.
    # For a bow shock, "outward" is generally away from Earth: ~+x near nose.
    if nhat[0] < 0:
        nhat = -nhat

    return that, nhat


def choose_foreshock_tangent(xbs, ybs, bhat, i_candidates):
    """
    From candidate tangency indices, choose the one with smallest thetaBn.
    Returns: i_best, theta_best_deg
    """
    best = None
    best_theta = None
    for i in i_candidates:
        if i is None:
            continue
        that, nhat = local_tangent_normal(xbs, ybs, i)
        if nhat is None:
            continue
        # thetaBn between IMF and normal
        cosang = np.clip(np.abs(np.dot(bhat, nhat)), -1.0, 1.0)
        theta = np.degrees(np.arccos(cosang))
        if (best_theta is None) or (theta < best_theta):
            best_theta = theta
            best = i
    return best, best_theta

def tangent_line_xy(x0, y0, bhat, smin=-80.0, smax=80.0, ns=200):
    """
    Parametric line through (x0,y0) along direction bhat in XY.
    s is in Re units.
    """
    s = np.linspace(smin, smax, ns)
    x = x0 + s * bhat[0]
    y = y0 + s * bhat[1]
    return x, y

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

def foreshock_mask_from_imf_xy(
    xbs, ybs,
    bx, by,
    theta_c_deg=45.0,
    require_dayside=True,
):
    """
    Identify quasi-parallel bow shock segment in the XY plane.

    Parameters
    ----------
    xbs, ybs : 1D arrays
        Bow shock curve in chosen frame, in Re.
    bx, by : float
        IMF components in the same frame as xbs/ybs (nT is fine; only direction matters).
    theta_c_deg : float
        Threshold for quasi-parallel: theta_Bn < theta_c.
    require_dayside : bool
        If True, only consider the dayside-ish region (x > 0) as foreshock source.

    Returns
    -------
    mask : bool array
        True where quasi-parallel.
    theta_bn_deg : 1D array
        Angle at each curve point (degrees).
    bhat : (2,) array
        Unit IMF direction in plane.
    nhat : (N,2) array
        Unit normals in plane (chosen to point outward-ish).
    """
    xbs = np.asarray(xbs, float)
    ybs = np.asarray(ybs, float)

    # IMF direction in plane
    bvec = np.array([bx, by], float)
    bnorm = np.linalg.norm(bvec)
    if not np.isfinite(bnorm) or bnorm == 0:
        raise ValueError(f"Invalid IMF plane vector: bx={bx}, by={by}")
    bhat = bvec / bnorm

    # Tangent along curve via gradient
    dx = np.gradient(xbs)
    dy = np.gradient(ybs)
    t = np.vstack([dx, dy]).T
    tnorm = np.linalg.norm(t, axis=1)
    tnorm[tnorm == 0] = np.nan
    that = t / tnorm[:, None]

    # Two candidate normals: rotate tangent by +/-90deg
    n1 = np.vstack([-that[:, 1], that[:, 0]]).T
    n2 = -n1

    # Pick the normal that points roughly outward (away from origin)
    r = np.vstack([xbs, ybs]).T
    pick_n1 = np.sum(n1 * r, axis=1) > np.sum(n2 * r, axis=1)
    nhat = np.where(pick_n1[:, None], n1, n2)

    # theta_Bn = acos(|n·b|)  (use abs because parallel vs anti-parallel both quasi-parallel)
    cosang = np.abs(np.sum(nhat * bhat[None, :], axis=1))
    cosang = np.clip(cosang, -1.0, 1.0)
    theta_bn = np.degrees(np.arccos(cosang))

    mask = theta_bn < float(theta_c_deg)

    if require_dayside:
        mask &= (xbs > 0)

    return mask, theta_bn, bhat, nhat

def read_hs_cdf_directory(directory: str, pos_var="Position", epoch_var="Epoch"):
    """
    Read and concatenate HelioSwarm representative trajectory CDF files.

    Returns:
      times_dt: list[datetime] (UTC-aware if possible)
      pos_re:   np.ndarray shape (Nt, Nsc, 3) in Earth radii (Re)
               frame: GSE (as stored in the CDF set)
    """
    files = sorted([fn for fn in os.listdir(directory) if fn.endswith(".cdf")])
    if not files:
        raise FileNotFoundError(f"No .cdf files found in {directory}")

    pos_list = []
    t_list = []

    for fn in files:
        fp = os.path.join(directory, fn)
        with pycdf.CDF(fp) as cdf:
            pos = np.array(cdf[pos_var][:])
            tt  = np.array(cdf[epoch_var][:])

        pos_list.append(pos)
        t_list.append(tt)

    pos = np.concatenate(pos_list, axis=0)
    tt  = np.concatenate(t_list, axis=0)

    # ---- Normalize time objects to python datetimes
    # SpacePy usually returns python datetime objects already.
    # We’ll coerce to list for consistency.
    times = list(tt)

    # ---- Normalize positions to shape (Nt, Nsc, 3)
    pos_norm = normalize_hs_position_array(pos)

    # ---- km -> Re
    # (HS file units are km; convert to Earth Radius)
    pos_re = pos_norm / RE_KM

    print(f"CDF Files Joined")
    
    return times, pos_re

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
        r_seeds = [2.5, 5.0, 10.]
        theta_max = np.deg2rad(85)
        thetas = np.linspace(-theta_max, theta_max, 7)
        seeds = [(r*np.cos(th), 0.0, r*np.sin(th)) for r in r_seeds for th in thetas]

        tail_max = np.deg2rad(35)
        tail = np.linspace(0.0, tail_max, 7)
        seeds += [
            (r*np.cos(np.pi - th), 0.0, sgn*r*np.sin(np.pi - th))
            for r in [3.5, 7.5]
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

def normalize_hs_position_array(pos):
    """
    Normalize HS 'Position' array to shape (Nt, Nsc, 3).

    Accepts a few common shapes:
      (Nt, Nsc, 3)
      (Nt, Nsc, 1, 3)
      (Nt, Nsc, 3, 1)
      (Nt, 1, Nsc, 3)

    Raises if it can’t interpret.
    """

    print(f"Standardizing Position Array")
    
    pos = np.asarray(pos)
    if pos.ndim == 3 and pos.shape[-1] == 3:
        # (Nt, Nsc, 3)
        return pos

    if pos.ndim == 4:
        # common oddball: (Nt, Nsc, 1, 3)
        if pos.shape[-1] == 3 and pos.shape[-2] == 1:
            return pos[..., 0, :]  # -> (Nt, Nsc, 3)

        # (Nt, Nsc, 3, 1)
        if pos.shape[-2] == 3 and pos.shape[-1] == 1:
            return pos[..., 0]     # -> (Nt, Nsc, 3)

        # (Nt, 1, Nsc, 3)
        if pos.shape[1] == 1 and pos.shape[-1] == 3:
            return pos[:, 0, :, :] # -> (Nt, Nsc, 3)

    raise ValueError(f"Unrecognized HS Position array shape: {pos.shape}")


def add_zoom_box(ax, center_xy, halfwidth_km, *, edgecolor="0.25", lw=1.2, ls="--", alpha=0.9, zorder=5):
    """
    Draw a square box centered at center_xy with halfwidth_km on a parent axis.
    center_xy: (2,) in km
    """
    cx, cy = float(center_xy[0]), float(center_xy[1])
    hw = float(halfwidth_km)

    rect = Rectangle(
        (cx - hw, cy - hw),  # lower-left
        2 * hw, 2 * hw,
        fill=False,
        edgecolor=edgecolor,
        linewidth=lw,
        linestyle=ls,
        alpha=alpha,
        zorder=zorder,
    )
    ax.add_patch(rect)
    return rect

def add_hub_zoom_inset(ax, hub_xy_re, nodes_xy_re, colors, *,
                       hub_color=HUB_COLOR,
                       halfwidth_km=200.0,
                       loc="upper right",
                       inset_size="25%",   # tweak: "35%".."45%"
                       title=None):
    """
    ax: parent axis
    hub_xy_re: (2,) hub position in that plane (Re)
    nodes_xy_re: (N,2) node positions in that plane (Re)
    colors: list of N matplotlib colors for nodes
    halfwidth_km: zoom half-width around hub in km
    """
    #hw = halfwidth_km / RE_KM  # in Re
    hw = halfwidth_km

    axins = inset_axes(ax, width=inset_size, height=inset_size, loc=loc, borderpad=1.0)

    nodes_xy_re = np.asarray(nodes_xy_re)

    # Force to shape (N, 2)
    if nodes_xy_re.ndim == 1 and nodes_xy_re.size == 2:
        nodes_xy_re = nodes_xy_re.reshape(1, 2)
    elif nodes_xy_re.ndim != 2 or nodes_xy_re.shape[1] != 2:
        raise ValueError(f"nodes_xy_re must be (N,2); got shape {nodes_xy_re.shape}")
    
    # plot nodes relative structure (same style as your relative panels)
    axins.plot(hub_xy_re[0], hub_xy_re[1], "o", ms=5, color=hub_color)  # hub marker
    for i in range(nodes_xy_re.shape[0]):
        axins.plot(nodes_xy_re[i, 0], nodes_xy_re[i, 1], "o", ms=4, color=colors[i])

    # nice tight limits around hub
    axins.set_xlim(hub_xy_re[0] - hw, hub_xy_re[0] + hw)
    axins.set_ylim(hub_xy_re[1] - hw, hub_xy_re[1] + hw)
    axins.set_aspect("equal", adjustable="box")

    # choose a nice spacing
    major_step = halfwidth_km / 2     # e.g. 250 if halfwidth=500
    minor_step = major_step / 2       # e.g. 125
    
    axins.xaxis.set_major_locator(MultipleLocator(major_step))
    axins.yaxis.set_major_locator(MultipleLocator(major_step))
    
    axins.xaxis.set_minor_locator(MultipleLocator(minor_step))
    axins.yaxis.set_minor_locator(MultipleLocator(minor_step))

    # draw the same zoom box on the parent axis
    add_zoom_box(ax, hub_xy_re, halfwidth_km, edgecolor="0.25", lw=1.2, ls="--", alpha=0.9)
    
    # optional cosmetics
    axins.grid(True, alpha=0.3)
    axins.tick_params(labelsize=7)
    if title:
        axins.set_title(title, fontsize=8)

    # make inset spines a touch stronger so it reads clearly
    for sp in axins.spines.values():
        sp.set_linewidth(1.0)

    return axins

def plot_hs_hub(
    t_drv: dt.datetime,
    hs_dir: str,
    t_hs_center: dt.datetime,
    window_days: float = 2.0,
    cadence_minutes: int = 60,
    outdir: str = "figs",
    xylim: float = 40.0,
):
    """
    Plot MP/BS + T96 field lines + HelioSwarm HUB trajectory window
    in two separate figures: GSM and GSE.
    """

    os.makedirs(outdir, exist_ok=True)

    # -------------------------
    # 1) Model snapshot (OMNI -> T96)
    # -------------------------
    if t_drv.tzinfo is None:
        t_drv = t_drv.replace(tzinfo=dt.timezone.utc)
    t_iso = t_drv.strftime("%Y-%m-%dT%H:%M:%S")

    Pdyn, Dst, bx_gse, by_gse, bz_gse, bx_gsm, by_gsm, bz_gsm, Ma = \
        get_omni2_params(t_iso, t_drv=t_drv)

    # T96 wants GSM By/Bz
    parmod = t96_parmod(Pdyn, Dst, by_gsm, bz_gsm)
    ut = unix_seconds(t_drv)

    # Boundaries in GSM (then we will transform for plotting)
    xmp_gsm, ymp_gsm, r0mp, alphamp = shue98_magnetopause_rtheta(Pdyn, bz_gsm, ntheta=361)
    xbs_gsm, ybs_gsm, rbs0, delta_bs = bowshock_rtheta_from_Ma(r0mp, Ma, ntheta=361)

    # -------------------------
    # 2) HelioSwarm ephemeris window (GSE base + GSM conversion)
    # -------------------------
    print(f"Extracting HS Ephemeris")
    times_hs, pos_gse, pos_gsm = get_hs_ephemeris_window(
        hs_dir,
        t_hs_center,
        t_drv,
        window_days=window_days,
        cadence_minutes=cadence_minutes,
    )

    moon_gse_re = get_moon_position(t_hs_center, frame="GSE", units="Re")
    moon_gsm_re = get_moon_position(t_hs_center, frame="GSM", units="Re")    

    # HUB only (spacecraft index 0)
    hub_gse = pos_gse[:, 0, :]  # (Nt, 3)
    nodes_gse = pos_gse[:, 1:9, :]      # (Nt, 8, 3)
    rel_gse = (nodes_gse - hub_gse[:, None, :])*RE_KM     # (Nt, 8, 3)
    
    hub_gsm = pos_gsm[:, 0, :]  # (Nt, 3)
    nodes_gsm = pos_gsm[:, 1:9, :]      # (Nt, 8, 3)
    rel_gsm = (nodes_gsm - hub_gsm[:, None, :])*RE_KM     # (Nt, 8, 3)

    # Choose a “marker” index near the HS center time (not the driver time)
    if t_hs_center.tzinfo is None:
        t_hs_center = t_hs_center.replace(tzinfo=dt.timezone.utc)
    i_center = int(np.argmin([abs((ti - t_hs_center).total_seconds()) for ti in times_hs]))

    # -------------------------
    # 3) Plot loop: GSM then GSE (separate figures)
    # -------------------------
    for frame in ("GSM", "GSE"):
        print(f"Calculating {frame}:")
        hub = hub_gsm if frame == "GSM" else hub_gse
        nodes = nodes_gsm if frame == "GSM" else nodes_gse
        rel = rel_gsm if frame == "GSM" else rel_gse
        moon_re = moon_gsm_re if frame == "GSM" else moon_gse_re

        # Field lines in this frame (trace in GSM internally)
        field_xy, ps_xy = trace_field_lines_for_plot(frame, "xy", parmod, ut, t_drv, rlim=60.0, r0=1.0)
        field_xz, ps_xz = trace_field_lines_for_plot(frame, "xz", parmod, ut, t_drv, rlim=60.0, r0=1.0)
        ps_last = ps_xz if ps_xz is not None else ps_xy

        # Transform boundaries into this frame for each plane
        xmp_xy, ymp_xy = convert_curve_xy_to_frame(xmp_gsm, ymp_gsm, frame, t_drv, plane="xy")
        xbs_xy, ybs_xy = convert_curve_xy_to_frame(xbs_gsm, ybs_gsm, frame, t_drv, plane="xy")

        xmp_xz, zmp_xz = convert_curve_xy_to_frame(xmp_gsm, ymp_gsm, frame, t_drv, plane="xz")
        xbs_xz, zbs_xz = convert_curve_xy_to_frame(xbs_gsm, ybs_gsm, frame, t_drv, plane="xz")
        
        # -------------------------
        # 4) Render: XY and XZ
        # -------------------------
        fig, axs = plt.subplots(2, 2, figsize=(12.5, 9.0), constrained_layout=True)
        ax_xy   = axs[0, 0]
        ax_xz   = axs[0, 1]
        ax_rxy  = axs[1, 0]
        ax_rxz  = axs[1, 1]

        # Boundaries
        ax_xy.plot(xmp_xy,  ymp_xy, "-.", lw=1.5, color=MP_color, label="Magnetopause")
        ax_xy.plot(xmp_xy, -ymp_xy, "-.", lw=1.5, color=MP_color)
        ax_xy.plot(xbs_xy,  ybs_xy, "--", lw=1.5, color=BS_color, label="Bow shock")
        ax_xy.plot(xbs_xy, -ybs_xy, "--", lw=1.5, color=BS_color)

        ax_xz.plot(xmp_xz,  zmp_xz, "-.", lw=1.5, color=MP_color)
        ax_xz.plot(xmp_xz, -zmp_xz, "-.", lw=1.5, color=MP_color)
        ax_xz.plot(xbs_xz,  zbs_xz, "--", lw=1.5, color=BS_color)
        ax_xz.plot(xbs_xz, -zbs_xz, "--", lw=1.5, color=BS_color)

        
        # IMF components in this frame (you provide these)
        if frame == "GSE":
            bx, by = bx_gse, by_gse
        else:
            bx, by = bx_gsm, by_gsm
        

        i_pos, i_neg, bhat, that, score = bowshock_tangent_indices_to_imf(xbs_xy, ybs_xy, bx, by)

        print(f"[{frame}] two points: {i_pos}, {i_neg}")
        
        # ensure IMF points upstream (toward -X) for drawing rays
        if bhat[0] > 0:
            bhat = -bhat

        #i_best, theta_best = choose_foreshock_tangent(xbs_xy, ybs_xy, bhat, [i_pos, i_neg])
        i_best, theta_best = choose_foreshock_tangent(xbs_xy, ybs_xy, bhat, [i_pos])
        print(f"[{frame}] foreshock tangency: i={i_best}, thetaBn~{theta_best:.1f} deg")

        if i_best is not None:
            xt, yt = tangent_line_xy(xbs_xy[i_best], ybs_xy[i_best], bhat, smin=0, smax=140)
            ax_xy.plot(xt, yt, ":", lw=1., color=FS_color, label="Foreshock")
            xt, yt = tangent_line_xy(xbs_xy[i_best], ybs_xy[i_best], -bhat, smin=0, smax=140)
            ax_xy.plot(xt, yt, ":", lw=1., color=FS_color)
        
        #print(f"Frame: {frame}; Bx {bx}, By {by}")
        #mask, theta_bn, bhat, nhat = foreshock_mask_from_imf_xy(
        #    xbs_xy, ybs_xy,
        #    bx, by,
        #    theta_c_deg=45.0,
        #    require_dayside=False,
        #)

        #print(
        #    f"[{frame}] IMF_xy=({bx:.2f},{by:.2f}) | "
        #    f"|B_xy|={np.hypot(bx,by):.2f} | "
        #    f"thetaBn min/med/max = {np.nanmin(theta_bn):.1f}/{np.nanmedian(theta_bn):.1f}/{np.nanmax(theta_bn):.1f} deg | "
      #      f"mask={np.count_nonzero(mask)}/{mask.size}"
      #  )

        # Foreshock proxy (put ON TOP of the bow shock)
        #ax_xy.plot(
        #    xbs_xy[mask],  ybs_xy[mask],
        #    lw=4.0, alpha=0.9, color="#1f77b4",
        #    zorder=10,
        #    label=r"Foreshock (proxy: $\theta_{Bn}<45^\circ$)"
        #)
        #ax_xy.plot(
        #    xbs_xy[mask], -ybs_xy[mask],
        #    lw=4.0, alpha=0.9, color="#1f77b4",
        #    zorder=10,
        #)

        ax_xy.plot(moon_re[0], moon_re[1], marker="o", fillstyle='left',
                   color='#96b5ad', ms=6, label="Moon")
        ax_xz.plot(moon_re[0], moon_re[2], marker="o", fillstyle='left',
                   color='#96b5ad', ms=6)

        
        # Optional: IMF direction arrow near the top-left of the plot
        #arrow_origin = np.array([0.85*xylim, 0.85*xylim]) * 0.0  # or choose a nicer anchored location
        # Better: place it at e.g. (-30, 30) if that’s in view
        #arrow_origin = np.array([-0.8*xylim, 0.8*xylim])
        arrow_origin = np.array([0.7*xylim, 0.8*xylim])
        arrow_len = 8.0  # Re, purely visual
        ax_xy.arrow(
            arrow_origin[0], arrow_origin[1],
            arrow_len*bhat[0], arrow_len*bhat[1],
            head_width=1.2, head_length=1.6,
            length_includes_head=True,
        )
        ax_xy.text(
            arrow_origin[0] + arrow_len*bhat[0],
            arrow_origin[1] + arrow_len*bhat[1],
            rf"$\mathbf{{B}}_{{IMF}}$ ({frame})",
            fontsize=9,
            ha="left", va="center"
        )
        
        # Field lines
        for X, Y, Z in field_xy:
            ax_xy.plot(X, Y, lw=0.6, alpha=0.8)
        for X, Y, Z in field_xz:
            ax_xz.plot(X, Z, lw=0.6, alpha=0.8)

        # HUB trajectory
        ax_xy.plot(hub[:, 0], hub[:, 1], lw=2, color=HUB_COLOR, label="HelioSwarm")
        ax_xy.plot(hub[i_center, 0], hub[i_center, 1], "o", color=HUB_COLOR)

        ax_xz.plot(hub[:, 0], hub[:, 2], lw=2, color=HUB_COLOR)
        ax_xz.plot(hub[i_center, 0], hub[i_center, 2], "o", color=HUB_COLOR)

        # Earth
        th = np.linspace(0, 2*np.pi, 361)
        ax_xy.plot(np.cos(th), np.sin(th), lw=2, color="k")
        ax_xz.plot(np.cos(th), np.sin(th), lw=2, color="k")

        # Axes cosmetics
        ax_xy.set_aspect("equal", adjustable="box")
        ax_xz.set_aspect("equal", adjustable="box")
        ax_xy.set_xlabel(rf"$X_{{{frame}}}\ [R_E]$")
        ax_xy.set_ylabel(rf"$Y_{{{frame}}}\ [R_E]$")
        ax_xz.set_xlabel(rf"$X_{{{frame}}}\ [R_E]$")
        ax_xz.set_ylabel(rf"$Z_{{{frame}}}\ [R_E]$")
        ax_xy.set_title(f"Equatorial plane ({frame})")
        ax_xz.set_title(f"Meridian plane ({frame})")
        ax_xy.grid(True); ax_xz.grid(True)

        ax_xy.set_xlim(-xylim, xylim); ax_xy.set_ylim(-xylim, xylim)
        ax_xz.set_xlim(-xylim, xylim); ax_xz.set_ylim(-xylim, xylim)

        fig.text(
            0.5, 0.45,
            #"Boundaries: axisymmetric empirical models.\n"
            "Axisymmetric Boundary Models: Shue 1998 MP, conic bow shock.\n"
            "Field: T96 + IGRF traced in GSM, plotted in selected frame.",
            ha="center",
            fontsize=8.5,
            style="italic",
            alpha=0.85,
        )        
        

        fig.suptitle(
            f"Driver: {t_drv.strftime('%Y-%m-%dT%H:%M:%SZ')} | HS window centered: {t_hs_center.strftime('%Y-%m-%dT%H:%M:%SZ')} | Frame={frame}\n"
            f"OMNI2: Pdyn={Pdyn:.2f} nPa, Dst={Dst:.0f} nT, By={by_gsm:.2f} nT, bz={bz_gsm:.2f} nT, M_A={Ma:.2f} | "
            f"MP r0={r0mp:.2f} Re α={alphamp:.2f} | BS r0={rbs0:.2f} Re (Δ={delta_bs:.2f} Re) | tilt(ps)={ps_last:.3f} rad",
            fontsize=10,
        )
        #ax_xy.legend(loc="lower right", frameon=False)
        #fix this location!
        ax_xy.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
            frameon=True,
        )
        
        # -------------------------
        # Relative swarm panels: nodes wrt hub
        # -------------------------
        # rel has shape (Nt, 8, 3) in Re
        node_colors = NODE_COLORS

        # Plot each node track with its assigned color
        for j in range(8):
            c = node_colors[j]
            ax_rxy.plot(rel[:, j, 0], rel[:, j, 1], lw=1.8, color=c)
            ax_rxz.plot(rel[:, j, 0], rel[:, j, 2], lw=1.8, color=c)
            
            # mark the center time point
            ax_rxy.plot(rel[i_center, j, 0], rel[i_center, j, 1], "o", ms=5, color=c)
            ax_rxz.plot(rel[i_center, j, 0], rel[i_center, j, 2], "o", ms=5, color=c)
            
        # hub is origin in relative plots
        ax_rxy.plot(0, 0, "o", ms=6, color=HUB_COLOR)
        ax_rxz.plot(0, 0, "o", ms=6, color=HUB_COLOR)

        ax_rxy.set_aspect("equal", adjustable="box")
        ax_rxz.set_aspect("equal", adjustable="box")
        ax_rxy.grid(True); ax_rxz.grid(True)

        ax_rxy.set_xlabel(rf"$\Delta X_{{{frame}}}\ [km]$")
        ax_rxy.set_ylabel(rf"$\Delta Y_{{{frame}}}\ [km]$")
        ax_rxz.set_xlabel(rf"$\Delta X_{{{frame}}}\ [km]$")
        ax_rxz.set_ylabel(rf"$\Delta Z_{{{frame}}}\ [km]$")
        
        ax_rxy.set_title("Nodes relative to Hub")
        ax_rxz.set_title("Nodes relative to Hub")

        nodes_xyz = rel[i_center, :, :]
        nodes_xy  = nodes_xyz[:, :2]
        nodes_xz  = nodes_xyz[:, [0, 2]]

        hw = 500.0
        inside_xy = np.sum((np.abs(nodes_xy[:,0]) <= hw) & (np.abs(nodes_xy[:,1]) <= hw))
        inside_xz = np.sum((np.abs(nodes_xz[:,0]) <= hw) & (np.abs(nodes_xz[:,1]) <= hw))
        print(f"Inset counts inside ±{hw:.0f} km: XY {inside_xy}/8, XZ {inside_xz}/8")
        
        add_hub_zoom_inset(
            ax_rxy,
            np.array([0.0, 0.0]),               # hub at origin
            nodes_xy,           # current time XY positions (km)
            node_colors,
            hub_color=HUB_COLOR,
            halfwidth_km=200.0,
            loc="upper right",
            title=""
        )
        
        add_hub_zoom_inset(
            ax_rxz,
            np.array([0.0, 0.0]),               # hub at origin
            nodes_xz,           # current time XZ positions (km)
            node_colors,
            hub_color=HUB_COLOR,
            halfwidth_km=200.0,
            loc="upper right",
            title=""
        )
        
        # plotting in km
        rel_lim = 3000.0
        ax_rxy.set_xlim(-rel_lim, rel_lim); ax_rxy.set_ylim(-rel_lim, rel_lim)
        ax_rxz.set_xlim(-rel_lim, rel_lim); ax_rxz.set_ylim(-rel_lim, rel_lim)

        outpng = os.path.join(outdir, f"HS_{t_hs_center.strftime('%Y%m%dT%H%M%SZ')}_{frame}.png")
        fig.savefig(outpng, dpi=200)
        plt.close(fig)
        print(f"Wrote {outpng}")

def imf_gse_to_gsm(bx_gse, by_gse, bz_gse, t_drv):
    v = np.array([[bx_gse, by_gse, bz_gse]], dtype=float)
    v_gsm = convert_xyz(v, [t_drv], "GSE", "GSM")[0]
    return v_gsm[0], v_gsm[1], v_gsm[2]
        
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

def trace_field_line(seed_xyz_gsm, parmod, ut_seconds, rlim=60.0, r0=1.0):
    ps = gp.recalc(ut_seconds)

    x0, y0, z0 = seed_xyz_gsm
    xf1, yf1, zf1, xx1, yy1, zz1 = gp.trace(x0, y0, z0, -1, rlim, r0, parmod, "t96", "igrf")
    xf2, yf2, zf2, xx2, yy2, zz2 = gp.trace(x0, y0, z0,  1, rlim, r0, parmod, "t96", "igrf")

    X = np.concatenate([np.asarray(xx1)[::-1], np.asarray(xx2)[1:]])
    Y = np.concatenate([np.asarray(yy1)[::-1], np.asarray(yy2)[1:]])
    Z = np.concatenate([np.asarray(zz1)[::-1], np.asarray(zz2)[1:]])
    return X, Y, Z, ps

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

def select_time_window_old(times, center_dt, window, cadence=None):
    """
    Select indices where |t - center| <= window.

    times: list[datetime]
    center_dt: datetime
    window: timedelta
    cadence: optional timedelta for downsampling (e.g., timedelta(hours=1))

    Returns: np.ndarray of selected indices (sorted)

    """
    if center_dt.tzinfo is None:
        center_dt = center_dt.replace(tzinfo=dt.timezone.utc)

    times = [
        (t.replace(tzinfo=dt.timezone.utc) if isinstance(t, dt.datetime) and t.tzinfo is None else t)
        for t in times
    ]

    diffs = np.array([abs((t - center_dt).total_seconds()) for t in times], float)
    mask = diffs <= window.total_seconds()
    idx = np.where(mask)[0]

    if idx.size == 0:
        raise RuntimeError("No HS points found in requested time window.")

    if cadence is None:
        return idx

    # cadence downsample: keep points spaced by >= cadence
    cad_s = cadence.total_seconds()
    keep = [idx[0]]
    last_t = times[idx[0]]
    for j in idx[1:]:
        if (times[j] - last_t).total_seconds() >= cad_s:
            keep.append(j)
            last_t = times[j]
    return np.array(keep, dtype=int)

def select_time_window(times, center_dt, window, cadence=None):
    """
    Select indices where center_dt - window <= t <= center_dt.

    times: list[datetime]
    center_dt: datetime
    window: timedelta
    cadence: optional timedelta for downsampling (e.g., timedelta(hours=1))

    Returns: np.ndarray of selected indices (sorted)
    """
    if center_dt.tzinfo is None:
        center_dt = center_dt.replace(tzinfo=dt.timezone.utc)

    times = [
        (t.replace(tzinfo=dt.timezone.utc)
         if isinstance(t, dt.datetime) and t.tzinfo is None else t)
        for t in times
    ]

    t0 = center_dt - window
    t1 = center_dt

    # boolean mask for [t0, t1]
    mask = np.array([(t0 <= t <= t1) for t in times], dtype=bool)
    idx = np.where(mask)[0]

    if idx.size == 0:
        raise RuntimeError("No HS points found in requested time window.")

    if cadence is None:
        return idx

    # cadence downsample: keep points spaced by >= cadence
    cad_s = cadence.total_seconds()
    keep = [idx[0]]
    last_t = times[idx[0]]
    for j in idx[1:]:
        if (times[j] - last_t).total_seconds() >= cad_s:
            keep.append(j)
            last_t = times[j]

    return np.array(keep, dtype=int)

def convert_hs_positions_frame(pos_re,
                               times,
                               from_sys="GSE",
                               to_sys="GSM",
                               t_drv=None,):
    """
    pos_re: (Nt, Nsc, 3) in Re
    times:  list[datetime] length Nt
    t_drv:  datetime to use for coordinate transform (optional)
            If provided, overrides HS times for transform.

    returns: (Nt, Nsc, 3) in Re
    """
    Nt, Nsc, _ = pos_re.shape
    flat = pos_re.reshape(Nt * Nsc, 3)

    if t_drv is not None:
        # Use driver time for all transforms (avoids IGRF 2029 warning)
        if t_drv.tzinfo is None:
            t_drv = t_drv.replace(tzinfo=dt.timezone.utc)

        times_rep = [t_drv] * (Nt * Nsc)
    else:
        # Original behavior: use HS timestamps
        times_rep = [t for t in times for _ in range(Nsc)]

    # repeat each time Nsc times to match flat rows
    #times_rep = [t for t in times for _ in range(Nsc)]

    flat2 = convert_xyz(flat, times_rep, from_sys, to_sys)  # uses your function
    return flat2.reshape(Nt, Nsc, 3)

def get_hs_ephemeris_window(hs_dir, t_hs_center, t_drv, window_days=2.0, cadence_minutes=60):
    """
    Returns:
      times_w: list[datetime] (selected window, downsampled)
      xyz_gse: (Nt, Nsc, 3) in Re
      xyz_gsm: (Nt, Nsc, 3) in Re
    """
    times_all, pos_all_gse = read_hs_cdf_directory(hs_dir)

    # Make CDF times UTC-aware (fixes naive/aware subtraction)
    times_all = [
        (t.replace(tzinfo=dt.timezone.utc) if isinstance(t, dt.datetime) and t.tzinfo is None else t)
        for t in times_all
    ]

    if isinstance(t_hs_center, str):
        center_dt = parse_iso_utc(t_hs_center)
    else:
        center_dt = t_hs_center
    if center_dt.tzinfo is None:
        center_dt = center_dt.replace(tzinfo=dt.timezone.utc)

    print("HS time sample:", times_all[0], "tzinfo:", times_all[0].tzinfo)
    print("Center time:", center_dt, "tzinfo:", center_dt.tzinfo)

    window = dt.timedelta(days=float(window_days))
    cadence = None if cadence_minutes is None else dt.timedelta(minutes=int(cadence_minutes))

    idx = select_time_window(times_all, center_dt, window, cadence=cadence)

    times = [times_all[i] for i in idx]
    pos_gse = pos_all_gse[idx, :, :]

    pos_gsm = convert_hs_positions_frame(pos_gse, times, "GSE", "GSM",t_drv=t_drv)

    return times, pos_gse, pos_gsm

if __name__ == "__main__":
    hs_dir = "/home/kgklein/Codes/HS-RT/PhB_SRD5B_0x75b"

    load_spice_kernels(kernel_dir="data/spice")
    # User-chosen solar wind / magnetosphere driver time (must be within OMNI2 coverage)
    t_drv = parse_iso_utc("1999-02-21T00:00:00Z")

    # Independent HS time (must fall within the HS trajectory CDF span)
    t_hs_center = parse_iso_utc("2029-08-07T06:00:00Z")

    plot_hs_hub(
        t_drv=t_drv,
        hs_dir=hs_dir,
        t_hs_center=t_hs_center,
        window_days=1.0,
        cadence_minutes=60,
        outdir="figs",
        xylim=65.0,
    )
