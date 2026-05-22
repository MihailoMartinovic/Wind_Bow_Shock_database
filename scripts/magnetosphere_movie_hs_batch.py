#!/usr/bin/env python3
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
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

import pyvista as pv

STATE_COLORS = {
    0: "#1f77b4",  # inside MP: blue
    1: "#1aa6a6",  # MP to BS: teal
    2: "#0af27e",  # foreshock-connected: green
    3: "#d62728",  # pristine SW: red
}

STATE_LABELS = {
    0: "Inside magnetopause",
    1: "Magnetosheath",
    2: "Foreshock-connected",
    3: "Pristine solar wind",
}

RE_KM = 6371.2  # Earth radius in km

FS_color='#00aa00'
MP_color="#a0a0a0"
BS_color="#326ba8"

#mp_color = (0.48, 0.56, 0.68)
#bs_color = (0.92, 0.78, 0.55)

mp_color = (0.55, 0.60, 0.65)   # soft cool gray-blue
mp_opacity = 0.28

bs_color = (0.85, 0.70, 0.45)   # warm tan
bs_opacity = 0.18

HS_COLORS = [
    "#56B4E9", "#E69F00", "#8E4D4D", "#888888",
    "#F0E442", "#D55E00", "#009E73", "#8D00B2", "#CC79A7"
]  # hub+8 nodes, matching hs-single-hour-baselines.py
HUB_COLOR = HS_COLORS[0]
NODE_COLORS = HS_COLORS[1:9]

def world_to_viewport(plotter, point):
    renderer = plotter.renderer
    coord = renderer.GetActiveCamera()

    # Use VTK coordinate transform
    import vtk
    coordinate = vtk.vtkCoordinate()
    coordinate.SetCoordinateSystemToWorld()
    coordinate.SetValue(point[0], point[1], point[2])

    display = coordinate.GetComputedDisplayValue(renderer)

    w, h = plotter.window_size
    return (display[0] / w, display[1] / h)

def boundary_curve_to_r_theta(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    theta = np.arctan2(y, x)
    theta = np.where(theta < 0.0, theta + np.pi, theta)
    r = np.sqrt(x * x + y * y)

    order = np.argsort(theta)
    return theta[order], r[order]

def revolve_r_theta_to_surface(theta, r_theta, n_phi=96, r_clip=None, x_min=None):
    theta = np.asarray(theta, float)
    r_theta = np.asarray(r_theta, float)

    if r_clip is not None:
        r_theta = np.clip(r_theta, None, r_clip)

    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=True)
    th, ph = np.meshgrid(theta, phi, indexing="ij")

    r = r_theta[:, None] * np.ones_like(ph)

    x = r * np.cos(th)
    rho = r * np.sin(th)
    y = rho * np.cos(ph)
    z = rho * np.sin(ph)

    grid = pv.StructuredGrid(x, y, z)
    surf = grid.extract_surface(algorithm="dataset_surface").triangulate()

    if x_min is not None:
        surf = surf.clip(origin=(x_min, 0.0, 0.0), normal=(1.0, 0.0, 0.0), invert=False)

    return surf



def transform_polydata_points(mesh, t_snap_dt, from_sys, to_sys):
    pts = np.asarray(mesh.points)
    pts2 = convert_xyz(pts, [t_snap_dt] * len(pts), from_sys, to_sys)
    out = mesh.copy(deep=True)
    out.points = pts2
    return out


def polyline_tube(points_xyz, radius=0.12, n_sides=16):
    poly = pv.lines_from_points(np.asarray(points_xyz, float))
    return poly.tube(radius=radius, n_sides=n_sides)

def make_colored_linecollection(x, y, states, state_colors, lw=2.5, alpha=1.0, zorder=6):
    """
    Build a LineCollection from x,y arrays, coloring each segment by state of endpoint i+1.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    states = np.asarray(states, int)

    pts = np.column_stack([x, y])
    segs = np.stack([pts[:-1], pts[1:]], axis=1)

    seg_colors = [state_colors[s] for s in states[1:]]

    lc = LineCollection(
        segs,
        colors=seg_colors,
        linewidths=lw,
        alpha=alpha,
        zorder=zorder,
        capstyle="round",
        joinstyle="round",
    )
    return lc

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

def mp_radius_shue1998(theta, r0mp, alphamp):
    """
    Magnetopause radius at polar angle theta from +X axis.
    theta in radians.
    """
    return r0mp * (2.0 / (1.0 + np.cos(theta))) ** alphamp


def bs_radius_conic(theta, rbs0, e=0.9):
    """
    Bow shock radius for the same conic form used in bowshock_rtheta_from_Ma.
    theta in radians, measured from +X axis.
    """
    return rbs0 * (1.0 + e) / (1.0 + e * np.cos(theta))

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

def cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]

def line_intersects_segment_xy(p0, d, q0, q1, tol=1e-10):
    """
    Does the infinite line p0 + u d intersect the segment q0->q1 ?
    p0, d, q0, q1 are 2-vectors.
    """
    p0 = np.asarray(p0, float)
    d  = np.asarray(d,  float)
    q0 = np.asarray(q0, float)
    q1 = np.asarray(q1, float)

    e = q1 - q0
    denom = cross2(d, e)

    if abs(denom) < tol:
        return False  # parallel or nearly so

    w = q0 - p0
    u = cross2(w, e) / denom
    t = cross2(w, d) / denom

    # infinite line in u, finite segment in t
    return (0.0 <= t <= 1.0)

def line_intersects_polyline_xy(x0, y0, bhat, xcurve, ycurve):
    """
    Test whether the infinite line through (x0,y0) parallel to bhat
    intersects the polyline (xcurve, ycurve).
    """
    p0 = np.array([x0, y0], float)
    d = normalize_vec2(bhat)

    for i in range(len(xcurve) - 1):
        q0 = np.array([xcurve[i],   ycurve[i]], float)
        q1 = np.array([xcurve[i+1], ycurve[i+1]], float)
        if line_intersects_segment_xy(p0, d, q0, q1):
            return True
    return False

def bowshock_connection_mask_xy(x, y, bhat, xbs, ybs):
    """
    For each point (x,y), return True if the infinite line through that point
    parallel to bhat intersects either branch of the bow shock.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)

    mask = np.zeros(len(x), dtype=bool)

    for i in range(len(x)):
        hit_upper = line_intersects_polyline_xy(x[i], y[i], bhat, xbs,  ybs)
        hit_lower = line_intersects_polyline_xy(x[i], y[i], bhat, xbs, -ybs)
        mask[i] = hit_upper or hit_lower

    return mask

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

def normalize_vec2(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    if n == 0 or not np.isfinite(n):
        raise ValueError(f"Cannot normalize vector {v}")
    return v / n


def point_side_of_line_xy(x, y, x0, y0, bhat):
    """
    Signed side of point(s) relative to an oriented line through (x0,y0)
    pointing along bhat = (bx,by).

    Positive/negative sign depends on orientation convention.
    Zero means on the line.
    """
    bhat = normalize_vec2(bhat)
    dx = np.asarray(x) - x0
    dy = np.asarray(y) - y0

    # 2D cross product z-component: b x r
    return bhat[0] * dy - bhat[1] * dx

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

def foreshock_side_sign(x0, y0, bhat, offset=5.0):
    """
    Determine which side of the tangent line is the foreshock side
    by probing a point a little upstream along bhat.
    """
    bhat = normalize_vec2(bhat)

    # choose upstream-pointing IMF direction (toward -X)
    bup = bhat.copy()
    if bup[0] > 0:
        bup = -bup

    x_ref = x0 + offset * bup[0]
    y_ref = y0 + offset * bup[1]

    sref = point_side_of_line_xy(x_ref, y_ref, x0, y0, bup)
    return np.sign(sref), bup

def foreshock_connected_mask_xy(x, y, xtan, ytan, bhat, tol=0.0):
    """
    Return boolean mask for points on the foreshock-connected side
    of the IMF tangent line in the XY plane.
    """
    side_sign, bup = foreshock_side_sign(xtan, ytan, bhat)

    s = point_side_of_line_xy(x, y, xtan, ytan, bup)
    return side_sign * s >= -tol

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

def classify_hub_positions_gsm(
    hub_gsm,
    r0mp,
    alphamp,
    rbs0,
    bhat_xy,
    xbs_gsm,
    ybs_gsm,
    bowshock_e=0.9,
):
    """
    Classify each hub position in GSM.

    States:
      0 = inside MP
      1 = between MP and BS
      2 = outside BS but magnetically connected to BS in XY (foreshock proxy)
      3 = pristine solar wind
    """
    hub_gsm = np.asarray(hub_gsm, float)
    x = hub_gsm[:, 0]
    y = hub_gsm[:, 1]
    z = hub_gsm[:, 2]

    r = np.sqrt(x*x + y*y + z*z)
    mu = np.clip(x / np.where(r > 0, r, np.nan), -1.0, 1.0)
    theta = np.arccos(mu)

    r_mp = mp_radius_shue1998(theta, r0mp, alphamp)
    r_bs = bs_radius_conic(theta, rbs0, e=bowshock_e)

    inside_mp = r <= r_mp
    inside_bs = r <= r_bs
    outside_bs = ~inside_bs

    # New: actual magnetic-connection proxy in XY
    connected_xy = bowshock_connection_mask_xy(x, y, bhat_xy, xbs_gsm, ybs_gsm)

    states = np.full(len(r), 3, dtype=int)     # default pristine solar wind
    states[outside_bs & connected_xy] = 2      # foreshock-connected
    states[inside_bs] = 1                      # magnetosheath / between MP and BS
    states[inside_mp] = 0                      # inside MP overrides all

    return states

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
        phi = np.linspace(0, 2*np.pi, 16, endpoint=True)
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
    tag=None,
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
    # 2b) GSM IMF direction for trajectory classification
    # -------------------------
    bhat_gsm_xy = normalize_vec2([bx_gsm, by_gsm])

    # -------------------------
    # 2c) Classify hub states once in GSM
    # -------------------------
    states = classify_hub_positions_gsm(
        hub_gsm,
        r0mp=r0mp,
        alphamp=alphamp,
        rbs0=rbs0,
        bhat_xy=bhat_gsm_xy,
        xbs_gsm=xbs_gsm,
        ybs_gsm=ybs_gsm,
        bowshock_e=0.9,
    )

    unique, counts = np.unique(states, return_counts=True)
    print("Hub state counts:", dict(zip(unique, counts)))
    
    #print(f"[GSM classification] tangent index={i_best_gsm}, thetaBn={theta_best_gsm:.1f} deg, "
    #      f"point=({xtan_gsm:.2f}, {ytan_gsm:.2f})")
    
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
        fig, axs = plt.subplots(2, 2, figsize=(12., 9.0), constrained_layout=True)
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
            ax_xy.plot(xt, yt, ":", lw=1., color=FS_color, label="Foreshock boundary")
            xt, yt = tangent_line_xy(xbs_xy[i_best], ybs_xy[i_best], -bhat, smin=0, smax=140)
            ax_xy.plot(xt, yt, ":", lw=1., color=FS_color)
        
        ax_xy.plot(moon_re[0], moon_re[1], marker="o", fillstyle='left',
                   color='#96b5ad', ms=6, label="Moon")
        ax_xz.plot(moon_re[0], moon_re[2], marker="o", fillstyle='left',
                   color='#96b5ad', ms=6)
        
        # Optional: IMF direction arrow near the top-right of the plot
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
        # HUB trajectory, colored by GSM-based region classification
        lc_xy = make_colored_linecollection(
            hub[:, 0], hub[:, 1],
            states,
            STATE_COLORS,
            lw=2.5,
            alpha=1.0,
            zorder=7,
        )
        ax_xy.add_collection(lc_xy)

        lc_xz = make_colored_linecollection(
            hub[:, 0], hub[:, 2],
            states,
            STATE_COLORS,
            lw=2.5,
            alpha=1.0,
            zorder=7,
        )
        ax_xz.add_collection(lc_xz)

        ax_xy.plot(hub[i_center, 0], hub[i_center, 1], "o", color=HUB_COLOR, ms=6, zorder=9)
        ax_xz.plot(hub[i_center, 0], hub[i_center, 2], "o", color=HUB_COLOR, ms=6, zorder=9)

        #ax_xy.plot(hub[:, 0], hub[:, 1], lw=2, color="#c01d14", label="HelioSwarm")
        #ax_xy.plot(hub[i_center, 0], hub[i_center, 1], "o", color="#fc0394")

        #ax_xz.plot(hub[:, 0], hub[:, 2], lw=2, color="#c01d14")
        #ax_xz.plot(hub[i_center, 0], hub[i_center, 2], "o", color="#fc0394")

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
            0.5, 0.475,
            #"Boundaries: axisymmetric empirical models.\n"
            "Axisymmetric Boundary Models:\n"
            "Shue 1998 MP + conic bow shock $(M_A)$.\n"
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
        #ax_xy.legend(loc="upper right", frameon=False)
        state_handles = [
            Line2D([0], [0], color=STATE_COLORS[0], lw=2.5, label="Inside MP"),
            Line2D([0], [0], color=STATE_COLORS[1], lw=2.5, label="MP to BS"),
            Line2D([0], [0], color=STATE_COLORS[2], lw=2.5, label="Foreshock-connected"),
            Line2D([0], [0], color=STATE_COLORS[3], lw=2.5, label="Pristine SW"),
        ]
        #ax_xy.legend(
        #    loc="upper left",
        #    bbox_to_anchor=(1.02, 1.0),
        #    borderaxespad=0.0,
        #    frameon=True,
        #)
        handles, labels = ax_xy.get_legend_handles_labels()
        handles.extend(state_handles)

        ax_xy.legend(handles=handles,
                     loc="upper right",
                     bbox_to_anchor=(1.75, 1.0),
                     borderaxespad=0.0,
                     frameon=True)
        
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

        outpng = os.path.join(outdir, f"HS_{t_hs_center.strftime('%Y%m%dT%H%M%SZ')}_{frame}_{tag}.png")
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

def make_relative_inset_png(
    rel_xyz_km,
    outfile,
    plane="xy",
    halfwidth_km=3000.0,
    node_colors=None,
    hub_color=HUB_COLOR,
    title=None,
    show_zoom_inset=True,
    zoom_halfwidth_km=200.0,
    zoom_box_on_parent=True,
):
    """
    Create a transparent PNG inset showing node positions relative to the hub,
    with an optional inset-within-inset zoom panel.

    Parameters
    ----------
    rel_xyz_km : ndarray, shape (N, 3)
        Relative spacecraft positions in km, typically rel_gse[i_center].
    outfile : str
        Output PNG path.
    plane : {"xy", "xz", "yz"}
        Projection plane for the inset.
    halfwidth_km : float
        Half-width of the main inset axes in km.
    node_colors : list of colors, optional
        One color per node.
    title : str, optional
        Title for the main inset.
    show_zoom_inset : bool
        Whether to add the smaller inset inside the larger inset.
    zoom_halfwidth_km : float
        Half-width of the nested zoom inset in km.
    zoom_box_on_parent : bool
        Whether to draw a dashed rectangle on the main inset marking the zoom region.
    """
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    rel_xyz_km = np.asarray(rel_xyz_km, float)
    if rel_xyz_km.ndim != 2 or rel_xyz_km.shape[1] != 3:
        raise ValueError("rel_xyz_km must have shape (N, 3)")

    plane = plane.lower()
    if plane == "xy":
        ii, jj = 0, 1
        xlabel, ylabel = r"$\Delta X$ [km]", r"$\Delta Y$ [km]"
    elif plane == "xz":
        ii, jj = 0, 2
        xlabel, ylabel = r"$\Delta X$ [km]", r"$\Delta Z$ [km]"
    elif plane == "yz":
        ii, jj = 1, 2
        xlabel, ylabel = r"$\Delta Y$ [km]", r"$\Delta Z$ [km]"
    else:
        raise ValueError("plane must be one of: 'xy', 'xz', 'yz'")

    if node_colors is None:
        node_colors = NODE_COLORS

    n_nodes = rel_xyz_km.shape[0]
    if len(node_colors) < n_nodes:
        raise ValueError("node_colors must have at least as many entries as nodes")

    xy = rel_xyz_km[:, [ii, jj]]

    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)

    fig = plt.figure(figsize=(3.3, 3.3), dpi=200, facecolor=(1, 1, 1, 0))
    ax = fig.add_axes([0.18, 0.12, 0.76, 0.80])
    ax.set_facecolor((1, 1, 1, 0.80))

    # -------------------------
    # Main inset: large scale
    # -------------------------
    ax.scatter(
        [0.0], [0.0],
        s=46, c=hub_color, marker="o",
        edgecolors="white", linewidths=0.6, zorder=5
    )

    for i in range(n_nodes):
        ax.scatter(
            xy[i, 0], xy[i, 1],
            s=28, c=node_colors[i], marker="o",
            edgecolors="black", linewidths=0.35, zorder=4
        )

    ax.set_xlim(-halfwidth_km, halfwidth_km)
    ax.set_ylim(-halfwidth_km, halfwidth_km)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(xlabel, fontsize=7,labelpad=2)
    ax.set_ylabel(ylabel, fontsize=7,labelpad=2)
    if title:
        ax.set_title(title, fontsize=8, pad=2)

    ax.tick_params(axis="both", which="major", labelsize=6, length=3)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="minor", length=2)
    ax.grid(True, which="major", alpha=0.25, linewidth=0.5)
    ax.grid(True, which="minor", alpha=0.12, linewidth=0.35)

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    # Dashed zoom box on parent
    if show_zoom_inset and zoom_box_on_parent:
        rect = Rectangle(
            (-zoom_halfwidth_km, -zoom_halfwidth_km),
            2.0 * zoom_halfwidth_km,
            2.0 * zoom_halfwidth_km,
            fill=False,
            linestyle="--",
            linewidth=0.9,
            edgecolor="0.35",
            zorder=6,
        )
        ax.add_patch(rect)

    # -------------------------
    # Nested zoom inset
    # -------------------------
    if show_zoom_inset:
        ax_zoom = ax.inset_axes([0.66, 0.60, 0.30, 0.30])
        ax_zoom.set_facecolor((1, 1, 1, 0.88))

        ax_zoom.scatter(
            [0.0], [0.0],
            s=30, c=hub_color, marker="o",
            edgecolors="white", linewidths=0.45, zorder=5
        )

        for i in range(n_nodes):
            ax_zoom.scatter(
                xy[i, 0], xy[i, 1],
                s=18, c=node_colors[i], marker="o",
                edgecolors="black", linewidths=0.25, zorder=4
            )

        ax_zoom.set_xlim(-zoom_halfwidth_km, zoom_halfwidth_km)
        ax_zoom.set_ylim(-zoom_halfwidth_km, zoom_halfwidth_km)
        ax_zoom.set_aspect("equal", adjustable="box")

        ax_zoom.tick_params(axis="both", which="major", labelsize=4.5, length=2)
        ax_zoom.minorticks_on()
        ax_zoom.tick_params(axis="both", which="minor", length=1.5)

        ax_zoom.grid(True, which="major", alpha=0.20, linewidth=0.4)
        ax_zoom.grid(True, which="minor", alpha=0.10, linewidth=0.25)

        for spine in ax_zoom.spines.values():
            spine.set_linewidth(0.6)

    fig.savefig(outfile, transparent=True)
    plt.close(fig)

    return outfile

def composite_multiple_insets_on_png(
    base_png,
    inset_specs,
    outfile=None,
):
    """
    Composite multiple transparent inset PNGs onto a base PNG.

    Parameters
    ----------
    base_png : str
        Path to main rendered frame.
    inset_specs : list of dict
        Each dict should contain:
            {
                "inset_png": <path>,
                "corner": "upper right" | "upper left" | "lower right" | "lower left",
                "pad_px": 30,          # optional
                "scale": 1.0,          # optional
            }
    outfile : str, optional
        Output file. If omitted, overwrite base_png.

    Returns
    -------
    outfile : str
        Output PNG path.
    """
    from PIL import Image

    if outfile is None:
        outfile = base_png

    base = Image.open(base_png).convert("RGBA")
    bw, bh = base.size

    composed = base.copy()

    for spec in inset_specs:
        inset_png = spec["inset_png"]
        corner = spec.get("corner", "upper right").lower().strip()
        pad_px = spec.get("pad_px", 30)
        scale = spec.get("scale", 1.0)

        inset = Image.open(inset_png).convert("RGBA")

        if scale is not None and scale != 1.0:
            iw, ih = inset.size
            new_size = (int(iw * scale), int(ih * scale))
            inset = inset.resize(new_size, resample=Image.LANCZOS)

        iw, ih = inset.size

        if corner == "upper right":
            x = bw - iw - pad_px
            y = pad_px
        elif corner == "upper left":
            x = pad_px
            y = pad_px
        elif corner == "lower right":
            x = bw - iw - pad_px
            y = bh - ih - pad_px
        elif corner == "lower left":
            x = pad_px
            y = bh - ih - pad_px
        else:
            raise ValueError(
                "corner must be one of: upper right, upper left, lower right, lower left"
            )

        composed.alpha_composite(inset, dest=(x, y))

    composed.save(outfile)
    return outfile

def composite_inset_on_png(
    base_png,
    inset_png,
    outfile=None,
    corner="upper right",
    pad_px=30,
    scale=0.35,   # NEW
):
    from PIL import Image

    if outfile is None:
        outfile = base_png

    base = Image.open(base_png).convert("RGBA")
    inset = Image.open(inset_png).convert("RGBA")

    # ---------------------------------
    # NEW: resize inset
    # ---------------------------------
    if scale is not None and scale != 1.0:
        iw, ih = inset.size
        new_size = (int(iw * scale), int(ih * scale))
        inset = inset.resize(new_size, resample=Image.LANCZOS)

    bw, bh = base.size
    iw, ih = inset.size  # updated after scaling

    corner = corner.lower().strip()
    if corner == "upper right":
        x = bw - iw - pad_px
        y = pad_px
    elif corner == "upper left":
        x = pad_px
        y = pad_px
    elif corner == "lower right":
        x = bw - iw - pad_px
        y = bh - ih - pad_px
    elif corner == "lower left":
        x = pad_px
        y = bh - ih - pad_px
    else:
        raise ValueError("corner must be one of: upper right, upper left, lower right, lower left")

    composed = base.copy()
    composed.alpha_composite(inset, dest=(x, y))
    composed.save(outfile)

    return outfile


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

def build_hs_3d_frame_data(
    t_drv: dt.datetime,
    hs_dir: str,
    t_hs_center: dt.datetime,
    window_days: float = 1.0,
    cadence_minutes: int = 60,
):
    """
    Build all data needed for one 3D GSE frame.

    Returns
    -------
    frame_data : dict
        Dictionary containing:
          - driver / IMF / boundary metadata
          - HS ephemeris window in GSE and GSM
          - HUB trajectory and current marker index
          - Moon position
          - MP/BS boundary curves in GSM
    """
    # -------------------------
    # 1) Model snapshot (OMNI -> T96)
    # -------------------------
    if t_drv.tzinfo is None:
        t_drv = t_drv.replace(tzinfo=dt.timezone.utc)
    if t_hs_center.tzinfo is None:
        t_hs_center = t_hs_center.replace(tzinfo=dt.timezone.utc)

    t_iso = t_drv.strftime("%Y-%m-%dT%H:%M:%S")

    Pdyn, Dst, bx_gse, by_gse, bz_gse, bx_gsm, by_gsm, bz_gsm, Ma = \
        get_omni2_params(t_iso, t_drv=t_drv)

    parmod = t96_parmod(Pdyn, Dst, by_gsm, bz_gsm)
    ut = unix_seconds(t_drv)

    # Boundaries in GSM
    xmp_gsm, ymp_gsm, r0mp, alphamp = shue98_magnetopause_rtheta(
        Pdyn, bz_gsm, ntheta=361
    )
    xbs_gsm, ybs_gsm, rbs0, delta_bs = bowshock_rtheta_from_Ma(
        r0mp, Ma, ntheta=361
    )

    # -------------------------
    # 2) HelioSwarm ephemeris window
    # -------------------------
    print("Extracting HS Ephemeris")
    times_hs, pos_gse, pos_gsm = get_hs_ephemeris_window(
        hs_dir,
        t_hs_center,
        t_drv,
        window_days=window_days,
        cadence_minutes=cadence_minutes,
    )

    moon_gse_re = get_moon_position(t_hs_center, frame="GSE", units="Re")
    moon_gsm_re = get_moon_position(t_hs_center, frame="GSM", units="Re")

    # HUB only
    hub_gse = pos_gse[:, 0, :]   # (Nt, 3)
    nodes_gse = pos_gse[:, 1:9, :]      # (Nt, 8, 3)
    rel_gse = (nodes_gse - hub_gse[:, None, :])*RE_KM     # (Nt, 8, 3)
    
    hub_gsm = pos_gsm[:, 0, :]   # (Nt, 3)
    nodes_gsm = pos_gsm[:, 1:9, :]      # (Nt, 8, 3)
    rel_gsm = (nodes_gsm - hub_gsm[:, None, :])*RE_KM     # (Nt, 8, 3)

    # Closest index to requested HS center time
    i_center = int(np.argmin([abs((ti - t_hs_center).total_seconds()) for ti in times_hs]))

    return {
        "t_drv": t_drv,
        "t_hs_center": t_hs_center,
        "t_iso": t_iso,
        "ut": ut,
        "parmod": parmod,
        "Pdyn": Pdyn,
        "Dst": Dst,
        "Ma": Ma,
        "bx_gse": bx_gse,
        "by_gse": by_gse,
        "bz_gse": bz_gse,
        "bx_gsm": bx_gsm,
        "by_gsm": by_gsm,
        "bz_gsm": bz_gsm,
        "xmp_gsm": xmp_gsm,
        "ymp_gsm": ymp_gsm,
        "xbs_gsm": xbs_gsm,
        "ybs_gsm": ybs_gsm,
        "r0mp": r0mp,
        "alphamp": alphamp,
        "rbs0": rbs0,
        "delta_bs": delta_bs,
        "times_hs": times_hs,
        "pos_gse": pos_gse,
        "pos_gsm": pos_gsm,
        "hub_gse": hub_gse,
        "hub_gsm": hub_gsm,
        "nodes_gse": nodes_gse,
        "nodes_gsm": nodes_gsm,
        "rel_gse": rel_gse,
        "rel_gsm": rel_gsm,
        "i_center": i_center,
        "moon_gse_re": moon_gse_re,
        "moon_gsm_re": moon_gsm_re,
    }


def render_hs_hub_3d_gse(
    frame_data,
    outfile,
    window_size=(1600, 1200),
    mp_r_clip=80.0,
    bs_r_clip=120.0,
    hub_radius=0.18,
    marker_radius=0.45,
    earth_radius=1.0,
    moon_radius=0.27,
    show_moon=True,
    show_axes=False,
    show_relative_inset=True,
    inset_planes=("xy", "xz"),
    inset_halfwidth_km=3000.0,
    inset_zoom_halfwidth_km=200.0,
    inset_corner="upper right",
    inset_pad_px=28,
    node_colors=None,
    hub_color=HUB_COLOR,
):

    """
    Render one 3D frame in GSE:
      - Earth
      - magnetopause surface
      - bow shock surface
      - hub trajectory (GSE)
      - current hub marker
      - optional Moon marker

    Writes a PNG to `outfile`.
    """

    if node_colors is None:
        node_colors = NODE_COLORS

    camera_position = (200.0, 125.0, 150.0)
    camera_focal_point = (-25.0, 0.0, 0.0)
    camera_viewup = (0.0, 0.0, 1.0)
    
    # -------------------------
    # 1) Build MP / BS surfaces in GSM
    # -------------------------
    theta_mp, r_mp = boundary_curve_to_r_theta(frame_data["xmp_gsm"], frame_data["ymp_gsm"])
    theta_bs, r_bs = boundary_curve_to_r_theta(frame_data["xbs_gsm"], frame_data["ybs_gsm"])

    #mp_surf_gsm = revolve_r_theta_to_surface(theta_mp, r_mp, n_phi=96, r_clip=mp_r_clip)
    #bs_surf_gsm = revolve_r_theta_to_surface(theta_bs, r_bs, n_phi=96, r_clip=bs_r_clip)

    mp_surf_gsm = revolve_r_theta_to_surface(
        theta_mp, r_mp, n_phi=96, r_clip=mp_r_clip, x_min=-120.0
    )
    bs_surf_gsm = revolve_r_theta_to_surface(
        theta_bs, r_bs, n_phi=96, r_clip=bs_r_clip, x_min=-160.0
    )
    
    # Transform surface vertices GSM -> GSE using the driver time
    t_drv = frame_data["t_drv"]
    mp_surf_gse = transform_polydata_points(mp_surf_gsm, t_drv, "GSM", "GSE")
    bs_surf_gse = transform_polydata_points(bs_surf_gsm, t_drv, "GSM", "GSE")

    # -------------------------
    # 2) Trajectory objects in GSE
    # -------------------------
    hub_gse = np.asarray(frame_data["hub_gse"], float)
    i_center = int(frame_data["i_center"])
    hub_now = hub_gse[i_center]

    traj_tube = polyline_tube(hub_gse, radius=hub_radius, n_sides=16)

    # -------------------------
    # 3) Scene setup
    # -------------------------
    pl = pv.Plotter(off_screen=True, window_size=list(window_size))
    pl.set_background("white")

    # Earth
    earth = pv.Sphere(radius=earth_radius, center=(0.0, 0.0, 0.0),
                      theta_resolution=64, phi_resolution=64)
    pl.add_mesh(earth, color="lightblue", smooth_shading=True)

    # Boundaries
    #pl.add_mesh(mp_surf_gse, color="slategray", opacity=0.23, smooth_shading=True)
    #pl.add_mesh(bs_surf_gse, color="tan", opacity=0.18, smooth_shading=True)
    pl.add_mesh(
        mp_surf_gse,
        color=mp_color,
        opacity=mp_opacity,
        smooth_shading=True,
        specular=0.15,
        specular_power=20,
    )
    
    pl.add_mesh(
        bs_surf_gse,
        color=bs_color,
        opacity=bs_opacity,
        smooth_shading=True,
        specular=0.10,
        specular_power=15,
    )
        
    # HUB trajectory and marker
    pl.add_mesh(traj_tube, color=hub_color, smooth_shading=True)

    hub_marker = pv.Sphere(radius=marker_radius, center=tuple(hub_now),
                           theta_resolution=48, phi_resolution=48)
    pl.add_mesh(hub_marker, color=hub_color, smooth_shading=True)

    # Optional Moon
    if show_moon and frame_data["moon_gse_re"] is not None:
        moon_xyz = np.asarray(frame_data["moon_gse_re"], float)

        print("Moon GSE (Re):", moon_xyz)
        print("Moon distance (Re):", np.linalg.norm(moon_xyz))
        
        moon = pv.Sphere(radius=moon_radius, center=tuple(moon_xyz),
                         theta_resolution=40, phi_resolution=40)
        pl.add_mesh(moon, color="#D9D9D9", smooth_shading=True)
        #pl.add_mesh(moon, color="0.85", smooth_shading=True)

    if show_axes:
        pl.add_axes()
    
    # -------------------------
    # 4) Camera
    # -------------------------
    # A stable, slightly sunward/above view in GSE
    pl.camera_position = [
        camera_position,
        camera_focal_point,
        camera_viewup,
    ]
    
    # Initialize render state before screenshot
    pl.show(auto_close=False)

    # Title text
    txt1 = f"GSE: {frame_data['t_hs_center'].strftime('%Y-%m-%d %H:%M UTC')}"
    txt2 = f"Drivers: {t_drv.strftime('%Y-%m-%dT%H:%M:%SZ')}"
    txt3 = f"Pdyn={frame_data['Pdyn']:.2f} nPa  Dst={frame_data['Dst']:.1f} nT  Ma={frame_data['Ma']:.2f}"    
    txt4 = f"Shue 1998 Magnetopause + conic Bowshock.\n"
    pl.add_text(
        txt1,
        position=(0.8, 0.87),   # just below logo
        viewport=True,
        font_size=12,
        color="black",
    )
    
    pl.add_text(
        txt2,
        position=(0.65, 0.07),
        viewport=True,
        font_size=10,
        color="black",
    )

    pl.add_text(
        txt3,
        position=(0.65, 0.05),
        viewport=True,
        font_size=10,
        color="black",
    )

    pl.add_text(
        txt4,
        position=(0.65, 0.01),
        viewport=True,
        font_size=10,
        color="black",
    )
    
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    pl.screenshot(outfile)
    pl.close()

    # ---------------------------------
    # Optional relative-position inset
    # ---------------------------------
    if show_relative_inset:
        i_center = int(frame_data["i_center"])

        # Prefer rel_gse if present; otherwise derive it from nodes/hub
        if "rel_gse" in frame_data:
            rel_now_km = np.asarray(frame_data["rel_gse"][i_center], float)
        else:
            hub_now_gse = np.asarray(frame_data["hub_gse"][i_center], float)
            nodes_now_gse = np.asarray(frame_data["nodes_gse"][i_center], float)
            rel_now_km = (nodes_now_gse - hub_now_gse[None, :]) * RE_KM

        inset_png = outfile.replace(".png", "_inset.png")

        xy_inset_png = outfile.replace(".png", "_inset_xy.png")
        xz_inset_png = outfile.replace(".png", "_inset_xz.png")
        
        make_relative_inset_png(
            rel_xyz_km=rel_now_km,
            outfile=xy_inset_png,
            plane="xy",
            halfwidth_km=3000.0,
            node_colors=node_colors,
            hub_color=hub_color,
            title="Relative Hub-Node Separations: XY",
            show_zoom_inset=True,
            zoom_halfwidth_km=200.0,
            zoom_box_on_parent=True,
        )
        
        make_relative_inset_png(
            rel_xyz_km=rel_now_km,
            outfile=xz_inset_png,
            plane="xz",
            halfwidth_km=3000.0,
            node_colors=node_colors,
            hub_color=hub_color,
            title="Relative Hub-Node Separations: XZ",
            show_zoom_inset=True,
            zoom_halfwidth_km=200.0,
            zoom_box_on_parent=True,
        )

        composite_multiple_insets_on_png(
            base_png=outfile,
            outfile=outfile,
            inset_specs=[
                {
                    "inset_png": xy_inset_png,
                    "corner": "upper left",
                    "pad_px": 28,
                    "scale": 0.8,
                },
                {
                    "inset_png": xz_inset_png,
                    "corner": "lower left",
                    "pad_px": 28,
                    "scale": 0.8,
                },
            ],
        )

        logo_png = "figs/HS001-logo_FinalBlackA.png"

        #composite_overlay_on_png(
        #    base_png=outfile,
        #    overlay_png=logo_png,
        #    outfile=outfile,
        #    corner="upper left",
        #    pad_px=18,
        #    scale=0.16,
        #)
        
        composite_inset_on_png(
            base_png=outfile,
            inset_png=logo_png,
            outfile=outfile,
            corner="upper right",
            pad_px=5,
            scale=0.12,   # <-- try 0.25–0.4
        )
        

    print(f"Wrote {outfile}")


def run_time_sequence_3d_hub_gse(
    t_drv_start,
    t_hs_start,
    hs_dir,
    n_steps=12,
    step_hours=1,
    window_days=1.0,
    cadence_minutes=60,
    outdir="movies",
):
    """
    March forward in time and render a PNG sequence of 3D GSE frames
    showing:
      - Earth
      - MP / BS
      - HUB orbit over the window
      - current HUB location

    Both:
      - OMNI driver time (t_drv)
      - HelioSwarm time (t_hs)

    advance together.
    """
    os.makedirs(outdir, exist_ok=True)

    dt_step = dt.timedelta(hours=float(step_hours))

    for i in range(n_steps):
        t_drv = t_drv_start + i * dt_step
        t_hs = t_hs_start + i * dt_step

        print(f"\n--- 3D Frame {i:03d} ---")
        print("Driver:", t_drv.isoformat())
        print("HS:", t_hs.isoformat())

        try:
            frame_data = build_hs_3d_frame_data(
                t_drv=t_drv,
                hs_dir=hs_dir,
                t_hs_center=t_hs,
                window_days=window_days,
                cadence_minutes=cadence_minutes,
            )

            outpng = os.path.join(
                outdir,
                f"HS3D_GSE_{t_hs.strftime('%Y%m%dT%H%M%SZ')}_j{i:03d}.png"
            )

            render_hs_hub_3d_gse(
                frame_data,
                outfile=outpng,
            )

        except Exception as e:
            print("Skipping frame:", e)

def run_time_sequence(
    t_drv_start,
    t_hs_start,
    hs_dir,
    n_steps=12,
    step_hours=1,
    window_days=1.0,
    cadence_minutes=60,
    outdir="figs",
    xylim=65.0,
):
    """
    March forward in time and generate a sequence of plots.

    Both:
        - OMNI driver time (t_drv)
        - HelioSwarm time (t_hs)

    advance together.

    Parameters
    ----------
    t_drv_start : datetime
    t_hs_start : datetime
    n_steps : int
        number of frames
    step_hours : float
        time step between frames
    """

    dt_step = dt.timedelta(hours=float(step_hours))

    for i in range(n_steps):
        t_drv = t_drv_start + i * dt_step
        t_hs  = t_hs_start + i * dt_step

        print(f"\n--- Frame {i:03d} ---")
        print("Driver:", t_drv.isoformat())
        print("HS:", t_hs.isoformat())

        try:
            plot_hs_hub(
                t_drv=t_drv,
                hs_dir=hs_dir,
                t_hs_center=t_hs,
                window_days=window_days,
                cadence_minutes=cadence_minutes,
                outdir=outdir,
                xylim=xylim,
                tag=f"j{i:03d}",   # optional filename tag
            )
        except Exception as e:
            print("Skipping frame:", e)

            
if __name__ == "__main__":

    hs_dir = "/home/kgklein/Codes/HS-RT/PhB_SRD5B_0x75b"

    load_spice_kernels(kernel_dir="data/spice")
    
    # starting times
    #t_drv_start = parse_iso_utc("1999-02-21T00:00:00Z")
    t_drv_start = parse_iso_utc("2020-01-01T00:00:00Z")
    t_hs_start  = parse_iso_utc("2029-08-07T06:00:00Z")

    frame_data = build_hs_3d_frame_data(
        t_drv=parse_iso_utc("2020-01-01T00:00:00Z"),
        hs_dir=hs_dir,
        t_hs_center=parse_iso_utc("2029-08-07T06:00:00Z"),
        window_days=1.0,
        cadence_minutes=60,
    )

    #render_hs_hub_3d_gse(frame_data, "test_hs3d.png")

    #render_hs_hub_3d_gse(
    #    frame_data,
    #    "test_hs3d.png",
    #    show_relative_inset=True,
    #    inset_halfwidth_km=3000.0,
    #)

    #run_time_sequence(
    #    t_drv_start=t_drv_start,
    #    t_hs_start=t_hs_start,
    #    hs_dir=hs_dir,
    #    n_steps=330,        # number of frames
    #    step_hours=2,      # cadence
    #    window_days=1.0,
    #    cadence_minutes=60,
    #    outdir="figs",
    #    xylim=65.0,
    #)


    run_time_sequence_3d_hub_gse(
        t_drv_start=t_drv_start,
        t_hs_start=t_hs_start,
        hs_dir=hs_dir,
        n_steps=4,
        step_hours=2,
        window_days=1.0,
        cadence_minutes=60,
        outdir="movies",
    )
