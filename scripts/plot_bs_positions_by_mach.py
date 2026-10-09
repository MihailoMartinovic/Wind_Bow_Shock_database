"""
Plot Wind bow shock crossing positions binned by Alfvén Mach number.

Produces three figures: GSE, GSM, and GPE coordinates.
Each figure has three panels:
  Left   — X-Y plane
  Centre — X-Z plane
  Right  — X-R plane  (R = sqrt(Y² + Z²))

GPE (Geocentric Plasma Ecliptic): GSE rotated around Z by the per-crossing
aberration angle α = arctan(V_Earth / V_sw), where V_Earth = 29.78 km/s.
V_sw is taken from OMNI2 (column vsw_km_s in mach_per_crossing.csv).

Crossing positions are marked with +.
Bow shock model curves from Merka et al. (2005), Table 2, GPE coordinates:
  3-D second-order surface: a1*x² + y² + a3*z² + 2*a4*xy + 2*a7*x + 2*a8*y + a10 = 0
  (north-south symmetry assumed: a2=1, a5=a6=a9=0)
  Pressure-normalized to n_a=7.0 cm⁻³, v_a=457.5 km/s.

Gray envelope = span between lowest (MA 2-5) and highest (MA 13-20) Mach bin curves.

Mach bins: MA = 2-5, 5-8, 8-13, 13-20  (Jelínek+ / Peredo+ / Merka+ convention)

Data sources:
  mach_per_crossing.csv   — crossing times + alfven_mach  (from add_mach_number.py)
  Wind_Ephemerids.h5      — GSE (key='GSE2', fallback '/GSE') and GSM (key='GSM2') positions

Usage (from bow_shock/ notebook):
    %run "Bow_Shock_Stats/plot_bs_positions_by_mach.py"
"""

import sys
import warnings
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR    = Path(__file__).resolve().parent
BOW_SHOCK_DIR = SCRIPT_DIR.parent

MACH_CSV     = SCRIPT_DIR / "mach_per_crossing.csv"
EPHEM_H5     = BOW_SHOCK_DIR / "Wind_Ephemerids.h5"

OUT_GSE        = BOW_SHOCK_DIR / "Figures" / "bs_positions_by_mach_GSE.png"
OUT_GSM        = BOW_SHOCK_DIR / "Figures" / "bs_positions_by_mach_GSM.png"
OUT_GPE        = BOW_SHOCK_DIR / "Figures" / "bs_positions_by_mach_GPE.png"
OUT_ABER       = BOW_SHOCK_DIR / "Figures" / "bs_aberration_comparison.png"
OUT_GSE_GEOTAIL = BOW_SHOCK_DIR / "Figures" / "bs_positions_geotail_GSE.png"
OUT_GSM_GEOTAIL = BOW_SHOCK_DIR / "Figures" / "bs_positions_geotail_GSM.png"
OUT_GPE_GEOTAIL = BOW_SHOCK_DIR / "Figures" / "bs_positions_geotail_GPE.png"
# Equal-count Mach bin figures
OUT_GSE_EQBIN         = BOW_SHOCK_DIR / "Figures" / "bs_positions_eqbins_GSE.png"
OUT_GSM_EQBIN         = BOW_SHOCK_DIR / "Figures" / "bs_positions_eqbins_GSM.png"
OUT_GPE_EQBIN         = BOW_SHOCK_DIR / "Figures" / "bs_positions_eqbins_GPE.png"
OUT_GSE_GEOTAIL_EQBIN = BOW_SHOCK_DIR / "Figures" / "bs_positions_geotail_eqbins_GSE.png"
OUT_GSM_GEOTAIL_EQBIN = BOW_SHOCK_DIR / "Figures" / "bs_positions_geotail_eqbins_GSM.png"
OUT_GPE_GEOTAIL_EQBIN = BOW_SHOCK_DIR / "Figures" / "bs_positions_geotail_eqbins_GPE.png"
OUT_GPE_GEOTAIL_EQBIN_4PANEL = BOW_SHOCK_DIR / "Figures" / "bs_positions_geotail_eqbins_GPE_4panel.png"
OUT_GSE_LINES    = BOW_SHOCK_DIR / "Figures" / "bs_positions_by_mach_GSE_lines.png"
OUT_GSE_STATS    = BOW_SHOCK_DIR / "Figures" / "bs_positions_by_mach_GSE_stats.png"
OUT_GSE_ALL      = BOW_SHOCK_DIR / "Figures" / "bs_positions_by_mach_GSE_all_styles.png"
OUT_GSE_BINTEST  = BOW_SHOCK_DIR / "Figures" / "bs_positions_by_mach_GSE_binsize_test.png"

V_EARTH = 29.78   # Earth's orbital speed [km/s]
RE_KM   = 6371.2  # Earth radius in km

# ---------------------------------------------------------------------------
# Mach bins  (edges, labels, line styles — mirrors the reference figure)
# ---------------------------------------------------------------------------

MACH_BINS   = [(2, 5), (5, 8), (8, 13), (13, 20)]
BIN_LABELS  = [r"$M_\mathrm{A} < 5$",
               r"$M_\mathrm{A} = 5$–$8$",
               r"$M_\mathrm{A} = 8$–$13$",
               r"$M_\mathrm{A} > 13$"]
BIN_COLORS  = ["#000000", "#d62728", "#2ca02c", "#1f77b4"]   # black, red, green, blue
BIN_MARKERS = ["o", "s", "^", "D"]                            # circle, square, triangle, diamond
BIN_STYLES  = [
    dict(lw=2.0, ls="-",  color=c) for c in BIN_COLORS
]

# ---------------------------------------------------------------------------
# Bow shock model — Merka et al. (2005), Table 2, GPE coordinates
# ---------------------------------------------------------------------------
# 3-D second-order surface:
#   a1*x^2 + y^2 + a3*z^2 + 2*a4*x*y + 2*a7*x + 2*a8*y + a10 = 0
# (a2=1 fixed; a5=a6=a9=0 for north-south symmetry)
# All coordinates in R_E, positions pressure-normalized to
#   n_a = 7.0 cm^-3,  v_a = 457.5 km/s
#
# Reference: Merka, J., A. Szabo, J. A. Slavin, and M. Peredo (2005),
#   J. Geophys. Res., 110, A04202, doi:10.1029/2004JA010944.

MERKA_GPE = {
    (2,  5):  dict(a1=-0.12292, a3=0.81092, a4=-0.20902, a7=28.270,  a8=-0.75985, a10=-649.29),
    (5,  8):  dict(a1=-0.02996, a3=0.93743, a4= 0.00227, a7=21.196,  a8=-0.57751, a10=-536.81),
    (8,  13): dict(a1=-0.21644, a3=0.91729, a4= 0.01680, a7=19.570,  a8= 0.02222, a10=-471.21),
    (13, 20): dict(a1= 0.04615, a3=1.00468, a4= 0.01286, a7=20.280,  a8=-0.03024, a10=-512.18),
}


# ---------------------------------------------------------------------------
# Merka coefficient interpolation — for arbitrary MA values / equal-count bins
# ---------------------------------------------------------------------------

# Representative MA for each Table 2 bin (arithmetic midpoints of 2-5,5-8,8-13,13-20)
_MERKA_MA_REFS   = np.array([3.5, 6.5, 10.5, 16.5])
_MERKA_COEFF_ARR = {
    name: np.array([MERKA_GPE[bk][name]
                    for bk in [(2,5),(5,8),(8,13),(13,20)]])
    for name in ("a1","a3","a4","a7","a8","a10")
}


def _interp_merka_coeffs(ma_value):
    """
    Linearly interpolate Merka 2005 Table 2 surface coefficients to an
    arbitrary Alfvén Mach number.  Clipped to the Table 2 range (3.5–16.5).
    """
    ma_c = float(np.clip(ma_value, _MERKA_MA_REFS[0], _MERKA_MA_REFS[-1]))
    return {name: float(np.interp(ma_c, _MERKA_MA_REFS, arr))
            for name, arr in _MERKA_COEFF_ARR.items()}


def build_equal_count_bins(mach_values, n_bins=4):
    """
    Split mach_values into n_bins equal-count bins.

    For each bin the median MA is used to interpolate Merka 2005 Table 2
    coefficients, so each bin gets its own model curve matched to the
    typical Mach number of that bin's population.

    Returns
    -------
    mach_bins  : list of (lo, hi) float tuples
    bin_labels : list of str label strings
    merka_gpe  : dict {(lo,hi): coeff_dict}  — interpolated coefficients
    """
    mach_clean = mach_values[np.isfinite(mach_values)]
    edges      = np.unique(np.nanpercentile(mach_clean, np.linspace(0, 100, n_bins + 1)))
    mach_bins  = [(float(edges[i]), float(edges[i+1])) for i in range(len(edges) - 1)]

    bin_labels = []
    for i, (lo, hi) in enumerate(mach_bins):
        if i == 0:
            bin_labels.append(rf"$M_\mathrm{{A}} < {hi:.1f}$")
        elif i == len(mach_bins) - 1:
            bin_labels.append(rf"$M_\mathrm{{A}} > {lo:.1f}$")
        else:
            bin_labels.append(rf"$M_\mathrm{{A}} = {lo:.1f}$–${hi:.1f}$")

    merka_gpe = {}
    for lo, hi in mach_bins:
        in_bin = (mach_clean >= lo) & (mach_clean <= hi)
        ma_rep = float(np.nanmedian(mach_clean[in_bin]))
        merka_gpe[(lo, hi)] = _interp_merka_coeffs(ma_rep)

    return mach_bins, bin_labels, merka_gpe


@contextmanager
def _override_bin_config(mach_bins, bin_labels, merka_gpe):
    """
    Temporarily replace the global bin configuration (MACH_BINS, BIN_LABELS,
    MERKA_GPE) so all figure functions pick up the new bins without needing
    parameter changes.  Restores the originals on exit.
    """
    global MACH_BINS, BIN_LABELS, MERKA_GPE
    old = MACH_BINS, BIN_LABELS, MERKA_GPE
    MACH_BINS, BIN_LABELS, MERKA_GPE = mach_bins, bin_labels, merka_gpe
    try:
        yield
    finally:
        MACH_BINS, BIN_LABELS, MERKA_GPE = old


def merka_xy(ma_bin_key, x_min=-40.0, x_max=22.0, n_points=3000):
    """
    X-Y plane (z=0) cross-section of the Merka 2005 bow shock.

    The sampling range defaults to the plotting window.  The geotail criterion asks for a
    wider one, because a crossing beyond the end of the sampled curve cannot be tested.

    Solving  a1*x^2 + y^2 + 2*a4*x*y + 2*a7*x + 2*a8*y + a10 = 0
    as a quadratic in y:
        y^2 + 2*(a4*x + a8)*y + (a1*x^2 + 2*a7*x + a10) = 0

    Returns (x_arr, y_dusk, y_dawn)  where y_dusk > 0 and y_dawn < 0.
    """
    p = MERKA_GPE[ma_bin_key]
    a1, a4, a7, a8, a10 = p["a1"], p["a4"], p["a7"], p["a8"], p["a10"]
    x    = np.linspace(x_min, x_max, n_points)
    b    = a4 * x + a8               # (1/2) coefficient of y  [quadratic half-coeff]
    c    = a1 * x**2 + 2*a7*x + a10  # constant term
    disc = b**2 - c
    ok   = disc >= 0
    sq   = np.sqrt(np.maximum(disc[ok], 0.0))
    y_p  = -b[ok] + sq               # dusk (positive y) branch
    y_m  = -b[ok] - sq               # dawn (negative y) branch
    return x[ok], y_p, y_m


def merka_xz(ma_bin_key):
    """
    X-Z plane (y=0) cross-section of the Merka 2005 bow shock.

    Solving  a1*x^2 + a3*z^2 + 2*a7*x + a10 = 0
    gives z = ±sqrt(-(a1*x^2 + 2*a7*x + a10) / a3)

    Returns (x_arr, z_north)  where z_north >= 0 (northern hemisphere).
    """
    p = MERKA_GPE[ma_bin_key]
    a1, a3, a7, a10 = p["a1"], p["a3"], p["a7"], p["a10"]
    x   = np.linspace(-40, 22, 3000)
    rhs = -(a1 * x**2 + 2*a7*x + a10) / a3
    ok  = rhs >= 0
    return x[ok], np.sqrt(rhs[ok])


def _merka_nose_x(ma_bin_key):
    """
    X coordinate of the Merka 2005 bow shock nose (subsolar point).

    The nose is where the XR curve touches R=0, i.e. where the discriminant
    of the quadratic in Y equals zero:
        (a4*x + a8)^2 - (a1*x^2 + 2*a7*x + a10) = 0
    which is a quadratic in x.  Returns the physically meaningful (smaller) root.
    """
    p = MERKA_GPE[ma_bin_key]
    a1, a4, a7, a8, a10 = p["a1"], p["a4"], p["a7"], p["a8"], p["a10"]
    A =  a4 ** 2 - a1
    B =  2.0 * (a4 * a8 - a7)
    C =  a8 ** 2 - a10
    disc = B ** 2 - 4 * A * C
    if disc < 0 or abs(A) < 1e-12:
        return np.nan
    sq = np.sqrt(disc)
    x1 = (-B + sq) / (2 * A)
    x2 = (-B - sq) / (2 * A)
    # Return the smaller positive root (the actual dayside nose; the other is
    # a large spurious root on the far side of the conic)
    candidates = [v for v in (x1, x2) if 0 < v < 100]
    return min(candidates) if candidates else np.nan


# Sampling range for the geotail reference curve.  It has to reach past the most distant
# crossing in the database (x_GPE = -94.7 R_E, on pass 69 outbound); -140 leaves room and the
# flag count is unchanged by going further out.
GEOTAIL_CURVE_X_MIN = -140.0

# Style of the optional extra curve a caller can overlay on each panel of the four-panel
# figure (used for a surface fitted to the crossings of that bin).
EXTRA_CURVE_COLOR, EXTRA_CURVE_STYLE = "black", "--"
EXTRA_CURVE_WIDTH, EXTRA_CURVE_ALPHA = 2.0, 0.85


def _merka_shifted_xr(ma_bin_key, x_min=GEOTAIL_CURVE_X_MIN):
    """
    Merka 2005 XR curve shifted so the nose is at X=0.

    Returns (x_shifted, R) arrays.  The nose lands at (0, 0) and the curve
    opens to negative x (the tail side).  Use this as the geotail boundary:
    crossings at X < 0 in the original frame are on the geotail side of the
    shifted curve.
    """
    xb, yp, _ = merka_xy(ma_bin_key, x_min=x_min)
    x_nose     = _merka_nose_x(ma_bin_key)
    if not np.isfinite(x_nose):
        x_nose = xb.max()           # fallback: use the max x in the curve
    return xb - x_nose, yp          # shift so nose → 0


def _envelope_fill(ax, bk_lo, bk_hi, plane, symmetric):
    """Fill gray envelope between lowest and highest MA bin curves."""
    x_grid = np.linspace(-30, 20, 1500)

    if plane in ("xy", "xr"):
        x_lo, yp_lo, _ = merka_xy(bk_lo)
        x_hi, yp_hi, _ = merka_xy(bk_hi)
        ya = np.interp(x_grid, x_lo, yp_lo, left=np.nan, right=np.nan)
        yb = np.interp(x_grid, x_hi, yp_hi, left=np.nan, right=np.nan)
    else:  # xz
        x_lo, zp_lo = merka_xz(bk_lo)
        x_hi, zp_hi = merka_xz(bk_hi)
        ya = np.interp(x_grid, x_lo, zp_lo, left=np.nan, right=np.nan)
        yb = np.interp(x_grid, x_hi, zp_hi, left=np.nan, right=np.nan)

    ok = np.isfinite(ya) & np.isfinite(yb)
    ymin = np.minimum(ya, yb)
    ymax = np.maximum(ya, yb)
    ax.fill_between(x_grid[ok], ymin[ok], ymax[ok],
                    alpha=0.20, color="gray", linewidth=0, zorder=0)
    if symmetric:
        ax.fill_between(x_grid[ok], -ymax[ok], -ymin[ok],
                        alpha=0.20, color="gray", linewidth=0, zorder=0)


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_mach_crossings():
    if not MACH_CSV.exists():
        print(f"ERROR: {MACH_CSV} not found.")
        print("Run add_mach_number.py first.")
        sys.exit(1)
    df = pd.read_csv(MACH_CSV, parse_dates=["time"])
    df = df.dropna(subset=["alfven_mach"]).reset_index(drop=True)
    print(f"✓ {len(df)} crossings with Alfvén Mach data loaded")
    return df


def load_positions_gse():
    """
    Load GSE spacecraft positions from Wind_Ephemerids.h5.
    Prefers the independently-fetched 'GSE2' key (from fetch_gse_sscweb.py);
    falls back to the original '/GSE' key if GSE2 is not yet available.
    """
    for key in ("GSE2", "/GSE"):
        try:
            df = pd.read_hdf(EPHEM_H5, key=key)
            if isinstance(df.columns, pd.MultiIndex):
                df.columns = ["_".join(str(s) for s in c).strip("_") for c in df.columns]
            print(f"  Using GSE key: '{key}'")
            return df
        except KeyError:
            continue
    raise KeyError("Neither 'GSE2' nor '/GSE' found in Wind_Ephemerids.h5")


def load_positions_gsm():
    """Load GSM spacecraft positions from Wind_Ephemerids.h5."""
    try:
        df = pd.read_hdf(EPHEM_H5, key="GSM2")
        if isinstance(df.columns, pd.MultiIndex):
            df.columns = ["_".join(str(s) for s in c).strip("_") for c in df.columns]
        return df
    except KeyError:
        print("⚠  GSM2 key not found in Wind_Ephemerids.h5 — GSM figure will be skipped.")
        return None


def find_pos_cols(df, coord):
    """Find X, Y, Z column names for a given coordinate system label."""
    coord = coord.upper()
    candidates = {
        "X": [f"{coord}_X [km]", f"{coord}_X[km]", f"X_{coord}", f"{coord}X"],
        "Y": [f"{coord}_Y [km]", f"{coord}_Y[km]", f"Y_{coord}", f"{coord}Y"],
        "Z": [f"{coord}_Z [km]", f"{coord}_Z[km]", f"Z_{coord}", f"{coord}Z"],
    }
    result = {}
    for axis, names in candidates.items():
        for n in names:
            if n in df.columns:
                result[axis] = n
                break
        if axis not in result:
            for col in df.columns:
                if axis in col.upper() and coord in col.upper():
                    result[axis] = col
                    break
    return result.get("X"), result.get("Y"), result.get("Z")


def match_positions(df_cross, df_pos, coord):
    """
    Match ephemeris positions to crossing times via time-index interpolation.

    Mirrors the notebook approach (00_Bow_Shock_Database_3D.ipynb):
      outer-merge crossings + ephemeris on time index, then
      interpolate(method='index') to get positions at exact crossing times.
    Returns df_cross with added columns X_RE, Y_RE, Z_RE.
    """
    cx, cy, cz = find_pos_cols(df_pos, coord)
    if cx is None:
        print(f"  ⚠  Could not find {coord} X/Y/Z columns in ephemeris.")
        print(f"     Available columns: {list(df_pos.columns[:10])} ...")
        return None

    # Slim ephemeris to just the 3 position columns
    pos_slim = df_pos[[cx, cy, cz]].loc[df_pos.index.notna()]

    # Outer merge: ephemeris timestamps + crossing timestamps in one frame
    cross_indexed = df_cross.set_index("time")
    merged = pd.merge(
        cross_indexed,
        pos_slim,
        left_index=True,
        right_index=True,
        how="outer"
    ).sort_index()

    # Linear interpolation of positions in time.
    # method='time' is correct for DatetimeIndex (method='index' breaks in pandas 2+)
    merged[[cx, cy, cz]] = merged[[cx, cy, cz]].interpolate(method="time")

    # Keep only crossing rows (identified by 'pass' column from df_cross)
    result = merged[merged["pass"].notna()].copy()
    result.index.name = "time"
    result = result.reset_index()

    result["X_RE"] = result[cx] / RE_KM
    result["Y_RE"] = result[cy] / RE_KM
    result["Z_RE"] = result[cz] / RE_KM
    return result


# ---------------------------------------------------------------------------
# GPE coordinate transformation
# ---------------------------------------------------------------------------

def gse_to_gpe(df):
    """
    Rotate GSE positions to GPE (Geocentric Plasma Ecliptic) per crossing.

    GPE X-axis points into the aberrated solar wind flow direction.
    Rotation around Z by α = arctan(V_Earth / V_sw), where:
      V_Earth = 29.78 km/s (Earth's orbital speed)
      V_sw    = OMNI solar wind speed [km/s]

    Rotation (GSE → GPE):
      X_GPE =  X_GSE * cos(α) + Y_GSE * sin(α)
      Y_GPE = -X_GSE * sin(α) + Y_GSE * cos(α)
      Z_GPE =  Z_GSE

    Falls back to α = 0 (no rotation) for missing V_sw.
    """
    if "vsw_km_s" not in df.columns:
        print("⚠  vsw_km_s not in data — re-run add_mach_number.py, using α=0 fallback")
        df = df.copy()
        df["X_RE_GPE"] = df["X_RE"]
        df["Y_RE_GPE"] = df["Y_RE"]
        df["Z_RE_GPE"] = df["Z_RE"]
        return df

    vsw   = df["vsw_km_s"].fillna(400.0)  # 400 km/s fallback
    alpha = np.arctan(V_EARTH / vsw)      # per-crossing aberration angle [rad]

    df = df.copy()
    df["alpha_deg"]  = np.degrees(alpha)
    df["X_RE_GPE"]   =  df["X_RE"] * np.cos(alpha) - df["Y_RE"] * np.sin(alpha)
    df["Y_RE_GPE"]   =  df["X_RE"] * np.sin(alpha) + df["Y_RE"] * np.cos(alpha)
    df["Z_RE_GPE"]   =  df["Z_RE"]

    print(f"  Aberration angle: {np.degrees(alpha).mean():.2f}° mean  "
          f"(range {np.degrees(alpha).min():.1f}°–{np.degrees(alpha).max():.1f}°)")
    return df


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def make_figure(df, coord_label, out_path):
    """
    Single-panel figure: X-R plane only (R = sqrt(Y²+Z²)).
    Bow shock curves: Merka et al. (2005), Table 2, GPE coordinates,
    plotted semi-transparently so data points read clearly.

    XY and XZ panel definitions are preserved below for future use —
    change `active_panels` to restore them.
    """
    df = df.copy()
    df["R_RE"] = np.sqrt(df["Y_RE"] ** 2 + df["Z_RE"] ** 2)

    # All panels defined; only XR is active in the output figure.
    # To restore multi-panel output set active_planes = {"xy", "xz", "xr"}
    # and change subplots to plt.subplots(1, 3, figsize=(15, 7)).
    _all_panels = [
        ("X_RE", "Y_RE", r"$Y\ [\mathrm{R}_\oplus]$",                          True,  "xy"),
        ("X_RE", "Z_RE", r"$Z\ [\mathrm{R}_\oplus]$",                          True,  "xz"),
        ("X_RE", "R_RE", r"$R = \sqrt{Y^2+Z^2}\ [\mathrm{R}_\oplus]$",     False, "xr"),
    ]
    panels = [p for p in _all_panels if p[4] == "xr"]

    fig, ax = plt.subplots(1, 1, figsize=(7, 7))
    axes = [ax]

    legend_handles = []

    for ax, (h_col, v_col, v_label, symmetric, plane) in zip(axes, panels):

        _envelope_fill(ax, MACH_BINS[0], MACH_BINS[-1], plane, symmetric)

        legend_handles = []
        for idx, ((lo, hi), label, style, color, marker) in enumerate(zip(
                MACH_BINS, BIN_LABELS, BIN_STYLES, BIN_COLORS, BIN_MARKERS)):
            bk = (lo, hi)
            if bk not in MERKA_GPE:
                continue

            # Model curves — semi-transparent
            curve_style = dict(style, alpha=0.35)
            if plane == "xy":
                xb, yp, ym = merka_xy(bk)
                ax.plot(xb, yp, **curve_style, zorder=2)
                ax.plot(xb, ym, **curve_style, zorder=2)
            elif plane == "xz":
                xb, zp = merka_xz(bk)
                ax.plot(xb,  zp, **curve_style, zorder=2)
                ax.plot(xb, -zp, **curve_style, zorder=2)
            else:  # xr
                xb, yp, _ = merka_xy(bk)
                ax.plot(xb, yp, **curve_style, zorder=2)

            # Use <= hi for the last bin to include the maximum value
            if idx == len(MACH_BINS) - 1:
                mask = (df["alfven_mach"] >= lo) & (df["alfven_mach"] <= hi)
            else:
                mask = (df["alfven_mach"] >= lo) & (df["alfven_mach"] < hi)
            sub  = df[mask]
            ax.scatter(
                sub[h_col], sub[v_col],
                marker=marker, s=30, linewidths=0.8,
                color=color, edgecolors=color, alpha=0.75, zorder=3
            )

            legend_handles.append(
                mlines.Line2D([], [], marker=marker, color=color,
                              markerfacecolor=color, markersize=6,
                              ls=style["ls"], lw=style["lw"],
                              label=f"{label}  (n={mask.sum()})")
            )

        _ax_format(ax, v_label, symmetric)

    axes[0].legend(handles=legend_handles, fontsize=9, loc="upper right",
                   framealpha=0.9, frameon=True)

    fig.text(0.99, 0.005,
             "Bow shock model: Merka et al. (2005) GPE, JGR 110, A04202",
             ha="right", va="bottom", fontsize=7.5, color="gray",
             style="italic")

    plt.suptitle(
        f"Wind Bow Shock Crossings by Alfvén Mach Number — {coord_label}",
        fontsize=12, fontweight="bold"
    )
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"✓ Saved: {out_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Aberration comparison figure
# ---------------------------------------------------------------------------

def plot_aberration_comparison(df_gse, df_gpe, out_path):
    """
    Three-panel figure showing the effect of the GSE → GPE aberration correction:
      Left   — Y distribution: GSE vs GPE (histogram, both overlaid)
      Centre — Z distribution: GSE vs GPE (should be identical — sanity check)
      Right  — Per-crossing aberration angle α distribution
    """
    fig, axes = plt.subplots(1, 3, figsize=(13, 5))

    bins_yz = np.arange(-50, 51, 3)

    # ----- Left: Y distribution -----
    ax = axes[0]
    ax.hist(df_gse["Y_RE"],     bins=bins_yz, alpha=0.6,
            color="#1f77b4", label="GSE", edgecolor="none")
    ax.hist(df_gpe["Y_RE_GPE"], bins=bins_yz, alpha=0.6,
            color="#ff7f0e", label="GPE", edgecolor="none")

    med_gse = df_gse["Y_RE"].median()
    med_gpe = df_gpe["Y_RE_GPE"].median()
    ax.axvline(med_gse, color="#1f77b4", lw=1.5, ls="--",
               label=f"GSE median = {med_gse:.1f} R$_E$")
    ax.axvline(med_gpe, color="#ff7f0e", lw=1.5, ls="--",
               label=f"GPE median = {med_gpe:.1f} R$_E$")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel(r"$Y\ [\mathrm{R}_\oplus]$", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Y distribution\n(dawn/dusk asymmetry)", fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()

    # ----- Centre: Z distribution (should not change) -----
    ax = axes[1]
    ax.hist(df_gse["Z_RE"],     bins=bins_yz, alpha=0.6,
            color="#1f77b4", label="GSE", edgecolor="none")
    ax.hist(df_gpe["Z_RE_GPE"], bins=bins_yz, alpha=0.6,
            color="#ff7f0e", label="GPE", edgecolor="none")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel(r"$Z\ [\mathrm{R}_\oplus]$", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Z distribution\n(sanity check — rotation around Z)", fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()

    # ----- Right: aberration angle distribution -----
    ax = axes[2]
    if "alpha_deg" in df_gpe.columns:
        alpha_vals = df_gpe["alpha_deg"].dropna()
        bins_a = np.arange(0, 12, 0.5)
        ax.hist(alpha_vals, bins=bins_a, color="#2ca02c", edgecolor="none", alpha=0.8)
        ax.axvline(alpha_vals.mean(),   color="k",   lw=1.5, ls="--",
                   label=f"Mean = {alpha_vals.mean():.2f}°")
        ax.axvline(alpha_vals.median(), color="gray", lw=1.5, ls=":",
                   label=f"Median = {alpha_vals.median():.2f}°")
        ax.set_xlabel("Aberration angle α [°]", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)
        ax.set_title("Per-crossing aberration angle\nα = arctan(V$_{Earth}$ / V$_{sw}$)",
                     fontsize=11)
        ax.legend(fontsize=9)
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.minorticks_on()
    else:
        ax.text(0.5, 0.5, "alpha_deg not available\n(re-run add_mach_number.py)",
                ha="center", va="center", transform=ax.transAxes)

    plt.suptitle(
        "GSE → GPE Aberration Correction\n"
        f"Δ(Y median) = {med_gpe - med_gse:+.2f} R$_E$",
        fontsize=12, fontweight="bold"
    )
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"✓ Saved: {out_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Helper: draw model curves onto an axis
# ---------------------------------------------------------------------------

def _draw_model_curves(ax, plane, alpha=0.35):
    """Overlay Merka 2005 curves semi-transparently."""
    for (lo, hi), style in zip(MACH_BINS, BIN_STYLES):
        bk = (lo, hi)
        if bk not in MERKA_GPE:
            continue
        s = dict(style, alpha=alpha)
        if plane == "xy":
            xb, yp, ym = merka_xy(bk)
            ax.plot(xb, yp, **s, zorder=2)
            ax.plot(xb, ym, **s, zorder=2)
        elif plane == "xz":
            xb, zp = merka_xz(bk)
            ax.plot(xb,  zp, **s, zorder=2)
            ax.plot(xb, -zp, **s, zorder=2)
        else:
            xb, yp, _ = merka_xy(bk)
            ax.plot(xb, yp, **s, zorder=2)


def _ax_format(ax, v_label, symmetric):
    ax.set_xlabel(r"$X\ [\mathrm{R}_\oplus]$", fontsize=14)
    ax.set_ylabel(v_label, fontsize=12)
    ax.set_xlim(-30, 20)
    ax.set_ylim(-45, 45) if symmetric else ax.set_ylim(0, 45)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()
    ax.set_aspect("equal", adjustable="box")
    ax.axhline(0, color="gray", lw=0.5, ls=":")
    ax.axvline(0, color="gray", lw=0.5, ls=":")


# ---------------------------------------------------------------------------
# Style B — connected-dot figure
# ---------------------------------------------------------------------------

def make_figure_lines(df, coord_label, out_path):
    """
    Like make_figure but data points per Mach bin are connected by a line
    (sorted by X position) with smaller markers — shows the spatial envelope
    of each bin as a ribbon rather than a cloud.
    """
    df = df.copy()
    df["R_RE"] = np.sqrt(df["Y_RE"] ** 2 + df["Z_RE"] ** 2)

    panels = [
        ("X_RE", "Y_RE", r"$Y\ [\mathrm{R}_\oplus]$",                          True,  "xy"),
        ("X_RE", "Z_RE", r"$Z\ [\mathrm{R}_\oplus]$",                          True,  "xz"),
        ("X_RE", "R_RE", r"$R = \sqrt{Y^2+Z^2}\ [\mathrm{R}_\oplus]$",     False, "xr"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 7))

    for ax, (h_col, v_col, v_label, symmetric, plane) in zip(axes, panels):
        _envelope_fill(ax, (2, 5), (13, 20), plane, symmetric)
        _draw_model_curves(ax, plane)

        legend_handles = []
        for (lo, hi), label, style, color, marker in zip(
                MACH_BINS, BIN_LABELS, BIN_STYLES, BIN_COLORS, BIN_MARKERS):
            mask = (df["alfven_mach"] >= lo) & (df["alfven_mach"] < hi)
            sub  = df[mask].sort_values(h_col)   # sort by X for coherent line
            ax.plot(
                sub[h_col].to_numpy(), sub[v_col].to_numpy(),
                color=color, lw=0.8, alpha=0.6,
                marker=marker, markersize=3, markerfacecolor=color,
                markeredgewidth=0.4, zorder=3
            )
            legend_handles.append(
                mlines.Line2D([], [], marker=marker, color=color,
                              markerfacecolor=color, markersize=5,
                              lw=0.8, label=f"{label}  (n={mask.sum()})")
            )

        _ax_format(ax, v_label, symmetric)

    axes[0].legend(handles=legend_handles, fontsize=9, loc="lower left", framealpha=0.9)
    fig.text(0.99, 0.005,
             "Bow shock model: Merka et al. (2005) GPE, JGR 110, A04202",
             ha="right", va="bottom", fontsize=7.5, color="gray", style="italic")
    plt.suptitle(
        f"Wind Bow Shock — {coord_label} — Style: connected dots (sorted by X)",
        fontsize=12, fontweight="bold"
    )
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Style C — binned mean / median figure
# ---------------------------------------------------------------------------

def _binned_stats(x, y, n_bins=12):
    """
    Bin x into n_bins equal-COUNT (quantile) intervals; return (centers, means, medians).

    Equal-count binning ensures every bin has roughly the same number of points,
    so more bins always yields more plotted points (until individual bins drop
    below min_count).  Equal-width bins behave the opposite way for sparse groups:
    narrow bins have 0–1 points and most fall below threshold.

    Bins with fewer than min_count points are dropped.
    """
    min_count = 3
    # Quantile edges → equal number of points per bin
    quantiles  = np.linspace(0, 100, n_bins + 1)
    bin_edges  = np.nanpercentile(x, quantiles)
    bin_edges  = np.unique(bin_edges)   # collapse duplicates (small datasets)

    centers, means, medians = [], [], []
    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        mask = (x >= lo) & (x <= hi)   # <= on right edge to include maximum
        if mask.sum() >= min_count:
            centers.append(0.5 * (lo + hi))
            means.append(np.nanmean(y[mask]))
            medians.append(np.nanmedian(y[mask]))
    return np.array(centers), np.array(means), np.array(medians)


def make_figure_stats(df, coord_label, out_path, n_bins=12):
    """
    Scatter plot with small semi-transparent markers (background) and overlaid
    binned-mean (solid) and binned-median (dashed) lines per Mach bin.
    """
    df = df.copy()
    df["R_RE"] = np.sqrt(df["Y_RE"] ** 2 + df["Z_RE"] ** 2)

    panels = [
        ("X_RE", "Y_RE", r"$Y\ [\mathrm{R}_\oplus]$",                          True,  "xy"),
        ("X_RE", "Z_RE", r"$Z\ [\mathrm{R}_\oplus]$",                          True,  "xz"),
        ("X_RE", "R_RE", r"$R = \sqrt{Y^2+Z^2}\ [\mathrm{R}_\oplus]$",     False, "xr"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 7))

    for ax, (h_col, v_col, v_label, symmetric, plane) in zip(axes, panels):
        _envelope_fill(ax, (2, 5), (13, 20), plane, symmetric)
        _draw_model_curves(ax, plane)

        legend_handles = []
        for (lo, hi), label, style, color, marker in zip(
                MACH_BINS, BIN_LABELS, BIN_STYLES, BIN_COLORS, BIN_MARKERS):
            mask = (df["alfven_mach"] >= lo) & (df["alfven_mach"] < hi)
            sub  = df[mask]

            # Background scatter — small, semi-transparent
            ax.scatter(sub[h_col], sub[v_col],
                       marker=marker, s=12, color=color,
                       alpha=0.30, linewidths=0, zorder=2)

            # Binned mean line (median kept in _binned_stats for future use)
            xc, mn, _md = _binned_stats(
                sub[h_col].values, sub[v_col].values, n_bins=n_bins
            )
            if len(xc) > 1:
                ax.plot(xc, mn, color=color, lw=3.0, ls="-", zorder=4,
                        label="_nolegend_")

            legend_handles.append(
                mlines.Line2D([], [], marker=marker, color=color,
                              markerfacecolor=color, markersize=5,
                              ls="-", lw=2.0,
                              label=f"{label}  (n={mask.sum()})  —mean  --median")
            )

        _ax_format(ax, v_label, symmetric)

    axes[0].legend(handles=legend_handles, fontsize=8.5, loc="lower left", framealpha=0.9)
    fig.text(0.99, 0.005,
             "Bow shock model: Merka et al. (2005) GPE, JGR 110, A04202",
             ha="right", va="bottom", fontsize=7.5, color="gray", style="italic")
    plt.suptitle(
        f"Wind Bow Shock — {coord_label} — Style: scatter + binned mean (—) / median (--)",
        fontsize=12, fontweight="bold"
    )
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# All-styles comparison — 3 rows × 3 panels (XY / XZ / XR)
# ---------------------------------------------------------------------------

def make_figure_all_styles(df, coord_label, out_path, n_bins=12):
    """
    3-row × 3-column figure comparing all display styles for one coordinate system:
      Row 0 — colored scatter (current default)
      Row 1 — connected dots sorted by X
      Row 2 — scatter + binned mean/median lines
    """
    df = df.copy()
    df["R_RE"] = np.sqrt(df["Y_RE"] ** 2 + df["Z_RE"] ** 2)

    panels = [
        ("X_RE", "Y_RE", r"$Y\ [\mathrm{R}_\oplus]$",                      True,  "xy"),
        ("X_RE", "Z_RE", r"$Z\ [\mathrm{R}_\oplus]$",                      True,  "xz"),
        ("X_RE", "R_RE", r"$R = \sqrt{Y^2+Z^2}\ [\mathrm{R}_\oplus]$", False, "xr"),
    ]
    row_titles = [
        "Style A — colored scatter",
        "Style B — connected dots (sorted by X)",
        "Style C — scatter + binned mean (—) / median (--)",
    ]

    fig, axes = plt.subplots(3, 3, figsize=(15, 19))

    for row, row_title in enumerate(row_titles):
        for col, (h_col, v_col, v_label, symmetric, plane) in enumerate(panels):
            ax = axes[row, col]
            _envelope_fill(ax, (2, 5), (13, 20), plane, symmetric)
            _draw_model_curves(ax, plane)

            legend_handles = []
            for (lo, hi), label, style, color, marker in zip(
                    MACH_BINS, BIN_LABELS, BIN_STYLES, BIN_COLORS, BIN_MARKERS):
                mask = (df["alfven_mach"] >= lo) & (df["alfven_mach"] < hi)
                sub  = df[mask]

                if row == 0:        # ---- Style A: scatter ----
                    ax.scatter(sub[h_col], sub[v_col],
                               marker=marker, s=25, color=color,
                               alpha=0.75, linewidths=0, zorder=3)
                    lh = mlines.Line2D([], [], marker=marker, color=color,
                                       markerfacecolor=color, markersize=5,
                                       ls="none",
                                       label=f"{label} (n={mask.sum()})")

                elif row == 1:      # ---- Style B: connected dots ----
                    sub_s = sub.sort_values(h_col)
                    ax.plot(sub_s[h_col].to_numpy(), sub_s[v_col].to_numpy(),
                            color=color, lw=0.8, alpha=0.6,
                            marker=marker, markersize=3,
                            markerfacecolor=color, markeredgewidth=0.3,
                            zorder=3)
                    lh = mlines.Line2D([], [], marker=marker, color=color,
                                       markerfacecolor=color, markersize=5,
                                       lw=0.8,
                                       label=f"{label} (n={mask.sum()})")

                else:               # ---- Style C: stats lines ----
                    ax.scatter(sub[h_col], sub[v_col],
                               marker=marker, s=10, color=color,
                               alpha=0.25, linewidths=0, zorder=2)
                    xc, mn, _md = _binned_stats(
                        sub[h_col].values, sub[v_col].values, n_bins=n_bins)
                    if len(xc) > 1:
                        ax.plot(xc, mn, color=color, lw=3.0, ls="-", zorder=4)
                    lh = mlines.Line2D([], [], marker=marker, color=color,
                                       markerfacecolor=color, markersize=5,
                                       ls="-", lw=2.0,
                                       label=f"{label} (n={mask.sum()})")

                legend_handles.append(lh)

            _ax_format(ax, v_label, symmetric)
            if col == 0:
                ax.set_title(row_title, fontsize=9, loc="left", pad=4)

    fig.text(0.99, 0.002,
             "Bow shock model: Merka et al. (2005) GPE, JGR 110, A04202",
             ha="right", va="bottom", fontsize=7.5, color="gray", style="italic")
    plt.suptitle(
        f"Wind Bow Shock — {coord_label} — All display styles",
        fontsize=13, fontweight="bold", y=1.002
    )

    # One shared legend from the last row's first panel
    axes[2, 0].legend(handles=legend_handles, fontsize=8.5,
                      loc="lower left", framealpha=0.9)

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Bin-size sensitivity — 2×3 grid, XY panel only, varying n_bins
# ---------------------------------------------------------------------------

BIN_SIZES_TEST = [4, 6, 8, 12, 16, 20]


def make_figure_binsize_test(df, coord_label, bin_sizes=BIN_SIZES_TEST):
    """
    2×3 grid of X-R (cylindrical) panels, one per bin count.
    Each panel: faint scatter background + binned mean (solid) / median (dashed).
    Plots inline only — no file saved.
    """
    df = df.copy()
    df["R_RE"] = np.sqrt(df["Y_RE"] ** 2 + df["Z_RE"] ** 2)

    n_rows, n_cols = 2, 3
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 10))
    axes_flat = axes.flatten()

    legend_handles = []
    for ax, n_bins in zip(axes_flat, bin_sizes):
        _envelope_fill(ax, (2, 5), (13, 20), "xr", symmetric=False)
        _draw_model_curves(ax, "xr")

        legend_handles = []
        for (lo, hi), label, style, color, marker in zip(
                MACH_BINS, BIN_LABELS, BIN_STYLES, BIN_COLORS, BIN_MARKERS):
            mask = (df["alfven_mach"] >= lo) & (df["alfven_mach"] < hi)
            sub  = df[mask]

            ax.scatter(sub["X_RE"].to_numpy(), sub["R_RE"].to_numpy(),
                       marker=marker, s=10, color=color,
                       alpha=0.20, linewidths=0, zorder=2)

            xc, mn, _md = _binned_stats(
                sub["X_RE"].to_numpy(), sub["R_RE"].to_numpy(), n_bins=n_bins)
            if len(xc) > 1:
                ax.plot(xc, mn, color=color, lw=3.0, ls="-", zorder=4)

            legend_handles.append(
                mlines.Line2D([], [], marker=marker, color=color,
                              markerfacecolor=color, markersize=5,
                              ls="-", lw=2.0,
                              label=f"{label} (n={mask.sum()})"))

        _ax_format(ax, r"$R = \sqrt{Y^2+Z^2}\ [\mathrm{R}_\oplus]$", symmetric=False)
        ax.set_title(f"n_bins = {n_bins}", fontsize=10, fontweight="bold")

    for ax in axes_flat[len(bin_sizes):]:
        ax.set_visible(False)

    axes_flat[0].legend(handles=legend_handles, fontsize=8.5,
                        loc="upper right", framealpha=0.9)
    fig.text(0.50, 0.01,
             "Solid = binned mean   |   Dashed = binned median   |   "
             "Bow shock model: Merka et al. (2005) GPE",
             ha="center", va="bottom", fontsize=8.5, color="gray", style="italic")
    plt.suptitle(
        f"Wind Bow Shock — {coord_label} X-R — Bin-size sensitivity",
        fontsize=13, fontweight="bold"
    )
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.show()


# ---------------------------------------------------------------------------
# Geotail-aware figure
# ---------------------------------------------------------------------------

def make_figure_geotail(df, coord_label, out_path, n_bins=12):
    """
    X-R figure distinguishing geotail (X<0) from dayside crossings.

    Geotail definition: X_RE < 0 (behind the terminator plane).

    Geotail boundary reference:
      The Merka 2005 MA 2-5 curve is the outermost expected bow shock.
      Its XR cross-section at X=0 gives the terminator radius R₀; for X<0
      the curve traces the expected flank/tail boundary.  This X≤0 portion
      is drawn as a bold dashed line — "the MA 2-5 curve shifted to start at
      X=0" — so it is visually clear which crossings sit inside this envelope.

    Scatter:
      filled markers  — dayside crossings (X ≥ 0)
      open markers    — geotail crossings  (X < 0)

    Binned mean (thick solid line): dayside crossings only.
    """
    df = df.copy()
    df["R_RE"] = np.sqrt(df["Y_RE"] ** 2 + df["Z_RE"] ** 2)

    # Geotail criterion: crossing lies INSIDE the shifted lowest-MA-bin boundary.
    # The shifted curve has its nose at (0, 0); it opens to negative X.
    # A crossing at (X, R) is geotail iff:
    #   X < 0  AND  R < R_shifted_merka_lowest(X)
    # (crossings outside the curve at X<0 are in the flank beyond the bow shock)
    xs, Rs = _merka_shifted_xr(MACH_BINS[0])
    x_vals = df["X_RE"].to_numpy()
    r_vals = df["R_RE"].to_numpy()
    # Interpolate the shifted curve R at each crossing's X position.
    # left/right=nan so points outside the curve's x-range are excluded.
    r_boundary = np.interp(x_vals, xs, Rs, left=np.nan, right=np.nan)
    df["_geotail"] = (x_vals < 0) & (r_vals < r_boundary) & np.isfinite(r_boundary)

    fig, ax = plt.subplots(1, 1, figsize=(7, 7))

    # Gray envelope (same as main figure)
    _envelope_fill(ax, MACH_BINS[0], MACH_BINS[-1], "xr", symmetric=False)

    # All Merka curves semi-transparent (background reference)
    _draw_model_curves(ax, "xr", alpha=0.20)

    # ---- Lowest-MA-bin Merka curve shifted so nose is at X=0 ---------------
    # Translates the lowest-MA (outermost) bow shock curve left by X_nose so
    # the subsolar point lands at X=0.  Crossings that fall inside this shifted
    # envelope at X<0 are geotail crossings.
    lo0, hi0 = MACH_BINS[0]
    x_nose   = _merka_nose_x(MACH_BINS[0])
    xs, Rs   = _merka_shifted_xr(MACH_BINS[0])
    ax.plot(xs, Rs,
            color="black", lw=1.5, ls=":", alpha=0.90, zorder=3,
            label=(rf"Merka 2005 $M_A$={lo0:.1f}–{hi0:.1f}  (nose shifted to X=0"
                   + (f", X$_{{nose}}$={x_nose:.1f} R$_E$)" if np.isfinite(x_nose) else ")")))

    # Vertical terminator line and nose marker
    ax.axvline(0, color="gray", lw=1.0, ls=":", alpha=0.8, zorder=1)
    ax.plot(0, 0, marker="*", color=BIN_COLORS[0], markersize=11,
            zorder=5, label="Shifted bow shock nose (X=0)")

    # ---- Per-Mach-bin scatter + dayside mean --------------------------------
    legend_handles = []
    for idx, ((lo, hi), label, style, color, marker) in enumerate(zip(
            MACH_BINS, BIN_LABELS, BIN_STYLES, BIN_COLORS, BIN_MARKERS)):

        if idx == len(MACH_BINS) - 1:
            mask = (df["alfven_mach"] >= lo) & (df["alfven_mach"] <= hi)
        else:
            mask = (df["alfven_mach"] >= lo) & (df["alfven_mach"] < hi)
        sub  = df[mask]
        day  = sub[~sub["_geotail"]]   # dayside
        tail = sub[sub["_geotail"]]    # inside shifted curve at X < 0

        # Dayside: filled markers
        if len(day) > 0:
            ax.scatter(
                day["X_RE"].to_numpy(), day["R_RE"].to_numpy(),
                marker=marker, s=30, color=color, edgecolors=color,
                linewidths=0.8, alpha=0.75, zorder=3
            )

        # Geotail: open (unfilled) markers — same color, no fill
        if len(tail) > 0:
            ax.scatter(
                tail["X_RE"].to_numpy(), tail["R_RE"].to_numpy(),
                marker=marker, s=30, facecolors="none", edgecolors=color,
                linewidths=0.9, alpha=0.75, zorder=3
            )

        # Binned mean from dayside only
        if len(day) >= 4:
            xc, mn, _md = _binned_stats(
                day["X_RE"].to_numpy(), day["R_RE"].to_numpy(), n_bins=n_bins
            )
            if len(xc) > 1:
                ax.plot(xc, mn, color=color, lw=3.0, ls="-", zorder=4)

        # Print counts to console; omit from legend label
        print(f"    {label}: day n={len(day)}, tail n={len(tail)}")
        legend_handles.append(
            mlines.Line2D(
                [], [], marker=marker, color=color,
                markerfacecolor=color, markersize=8,
                ls="-", lw=style["lw"],
                label=label
            )
        )

    _ax_format(ax, r"$R = \sqrt{Y^2+Z^2}\ [\mathrm{R}_\oplus]$", symmetric=False)

    # Larger axis labels and ticks for publication
    ax.set_xlabel(r"$X\ [\mathrm{R}_\oplus]$", fontsize=18)
    ax.set_ylabel(r"$R = \sqrt{Y^2+Z^2}\ [\mathrm{R}_\oplus]$", fontsize=18)
    ax.tick_params(axis="both", which="major", labelsize=15)

    # Legend: Mach bins only (no Merka / geotail / dayside annotation lines)
    ax.legend(
        handles=legend_handles,
        fontsize=15, loc="upper right", framealpha=0.9
    )

    fig.text(0.99, 0.005,
             "Bow shock model: Merka et al. (2005) GPE, JGR 110, A04202",
             ha="right", va="bottom", fontsize=7.5, color="gray", style="italic")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"✓ Saved: {out_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Geotail-aware 4-panel figure — one panel per Mach bin
# ---------------------------------------------------------------------------

def make_figure_geotail_4panel(df, coord_label, out_path, n_bins=12,
                               extra_curves=None, extra_label=None):
    """
    4-panel (2×2) version of make_figure_geotail: identical geotail/dayside
    logic and data, but each Mach bin gets its own panel.

    Every panel shows:
      • the gray envelope (span of lowest→highest MA bin curves),
      • that bin's own Merka 2005 GPE model curve (solid, bin color),
      • the lowest-MA-bin curve shifted so its nose sits at X=0 (dotted black)
        — the geotail boundary reference, same in every panel,
      • that bin's crossings: filled markers = dayside, open markers = geotail
        (X<0 inside the shifted boundary),
      • that bin's dayside binned-mean line (thick solid).
    """
    df = df.copy()
    df["R_RE"] = np.sqrt(df["Y_RE"] ** 2 + df["Z_RE"] ** 2)

    # Geotail criterion — identical to make_figure_geotail
    xs0, Rs0 = _merka_shifted_xr(MACH_BINS[0])
    x_vals = df["X_RE"].to_numpy()
    r_vals = df["R_RE"].to_numpy()
    r_boundary = np.interp(x_vals, xs0, Rs0, left=np.nan, right=np.nan)
    df["_geotail"] = (x_vals < 0) & (r_vals < r_boundary) & np.isfinite(r_boundary)

    x_nose = _merka_nose_x(MACH_BINS[0])

    n = len(MACH_BINS)
    n_cols = 2
    n_rows = int(np.ceil(n / n_cols))
    # Figure sized so each equal-aspect panel exactly fills its grid cell, so
    # panels touch with zero spacing.  Axes span X∈(-30,20)=50 R_E, R∈(0,45)=45.
    panel_h = 7.0
    x_range, y_range = 50.0, 45.0
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * panel_h * x_range / y_range, n_rows * panel_h),
        sharex=True, sharey=True,
        gridspec_kw=dict(wspace=0.0, hspace=0.0),
    )
    axes_flat = np.atleast_1d(axes).flatten()

    for ax, idx, ((lo, hi), label, style, color, marker) in zip(
            axes_flat, range(n),
            zip(MACH_BINS, BIN_LABELS, BIN_STYLES, BIN_COLORS, BIN_MARKERS)):
        bk = (lo, hi)

        # Gray reference envelope
        _envelope_fill(ax, MACH_BINS[0], MACH_BINS[-1], "xr", symmetric=False)

        # This bin's own Merka curve (solid)
        if bk in MERKA_GPE:
            xb, yp, _ = merka_xy(bk)
            ax.plot(xb, yp, **dict(style, alpha=0.9), zorder=2)

        # Optional curve supplied by the caller for this bin, e.g. a surface fitted to its
        # crossings.  Drawn only over the x range the caller provides.
        if extra_curves is not None and bk in extra_curves:
            ex, ey = extra_curves[bk]
            ax.plot(ex, ey, color=EXTRA_CURVE_COLOR, ls=EXTRA_CURVE_STYLE,
                    lw=EXTRA_CURVE_WIDTH, alpha=EXTRA_CURVE_ALPHA, zorder=6)

        # Shifted lowest-bin nose curve = geotail boundary reference
        ax.plot(xs0, Rs0, color="black", lw=1.5, ls=":", alpha=0.90, zorder=3)
        ax.axvline(0, color="gray", lw=1.0, ls=":", alpha=0.8, zorder=1)
        ax.plot(0, 0, marker="*", color=BIN_COLORS[0], markersize=11, zorder=5)

        # This bin's crossings (identical masking to make_figure_geotail)
        if idx == n - 1:
            mask = (df["alfven_mach"] >= lo) & (df["alfven_mach"] <= hi)
        else:
            mask = (df["alfven_mach"] >= lo) & (df["alfven_mach"] < hi)
        sub  = df[mask]
        day  = sub[~sub["_geotail"]]
        tail = sub[sub["_geotail"]]

        if len(day) > 0:
            ax.scatter(
                day["X_RE"].to_numpy(), day["R_RE"].to_numpy(),
                marker=marker, s=30, color=color, edgecolors=color,
                linewidths=0.8, alpha=0.75, zorder=3
            )
        if len(tail) > 0:
            ax.scatter(
                tail["X_RE"].to_numpy(), tail["R_RE"].to_numpy(),
                marker=marker, s=30, facecolors="none", edgecolors=color,
                linewidths=0.9, alpha=0.75, zorder=3
            )

        # Binned mean from dayside only
        if len(day) >= 4:
            xc, mn, _md = _binned_stats(
                day["X_RE"].to_numpy(), day["R_RE"].to_numpy(), n_bins=n_bins
            )
            if len(xc) > 1:
                ax.plot(xc, mn, color=color, lw=3.0, ls="-", zorder=4)

        print(f"    {label}: day n={len(day)}, tail n={len(tail)}")

        # Larger axis labels and ticks for publication (matches make_figure_geotail).
        # Labels only on the bottom row (X) and left column (R); shared axes
        # already hide the inner tick labels.
        _ax_format(ax, r"$R = \sqrt{Y^2+Z^2}\ [\mathrm{R}_\oplus]$", symmetric=False)
        row, col = idx // n_cols, idx % n_cols
        ax.set_xlabel(r"$X\ [\mathrm{R}_\oplus]$" if row == n_rows - 1 else "",
                      fontsize=18)
        ax.set_ylabel(r"$R = \sqrt{Y^2+Z^2}\ [\mathrm{R}_\oplus]$" if col == 0 else "",
                      fontsize=18)
        ax.tick_params(axis="both", which="major", labelsize=15)

        # Per-panel legend (this bin only) — no title, mirrors single-panel style
        handles = [mlines.Line2D([], [], marker=marker, color=color,
                                 markerfacecolor=color, markersize=8,
                                 ls="-", lw=style["lw"], label=label)]
        if extra_curves is not None and bk in extra_curves and extra_label:
            handles.append(mlines.Line2D([], [], color=EXTRA_CURVE_COLOR,
                                         ls=EXTRA_CURVE_STYLE, lw=EXTRA_CURVE_WIDTH,
                                         label=extra_label))
        ax.legend(handles=handles, fontsize=15, loc="upper right", framealpha=0.9)

    # Prune the tick labels that sit on the shared seams so adjacent panels'
    # edge numbers don't overprint each other (e.g. 0/45 at the row seam).
    from matplotlib.ticker import MaxNLocator
    for idx, ax in enumerate(axes_flat[:n]):
        row, col = idx // n_cols, idx % n_cols
        ax.yaxis.set_major_locator(
            MaxNLocator(nbins=9, steps=[1, 2, 2.5, 5, 10],
                        prune=("lower" if row == 0 else "upper")))
        ax.xaxis.set_major_locator(
            MaxNLocator(nbins=9, steps=[1, 2, 2.5, 5, 10],
                        prune=("upper" if col == 0 else "lower")))

    # Hide any unused panels
    for ax in axes_flat[n:]:
        ax.set_visible(False)

    fig.text(0.99, 0.005,
             "Bow shock model: Merka et al. (2005) GPE, JGR 110, A04202",
             ha="right", va="bottom", fontsize=7.5, color="gray", style="italic")
    # Keep panels touching (gridspec wspace/hspace=0); no tight_layout so the
    # zero spacing is preserved.
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"✓ Saved: {out_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 65)
    print("BOW SHOCK POSITIONS BY ALFVÉN MACH NUMBER")
    print("Bow shock model: Merka et al. (2005) GPE, JGR 110, A04202")
    print("=" * 65)

    df_mach = load_mach_crossings()

    print(f"\nMach bin counts:")
    for lo, hi in MACH_BINS:
        n = ((df_mach["alfven_mach"] >= lo) & (df_mach["alfven_mach"] < hi)).sum()
        print(f"  MA {lo:2d}–{hi:2d}: {n} crossings")

    # --- Clean up old exploratory figures ---
    for _p in [OUT_GSE_LINES, OUT_GSE_STATS, OUT_GSE_ALL, OUT_GSE_BINTEST]:
        try:
            _p.unlink()
            print(f"  Removed {_p.name}")
        except FileNotFoundError:
            pass

    df_gse_matched = None
    df_gpe_raw     = None
    df_gpe_plot    = None
    df_gsm         = None

    # --- GSE ---
    print(f"\nLoading GSE positions...")
    try:
        df_gse_pos     = load_positions_gse()
        df_gse_matched = match_positions(df_mach, df_gse_pos, "GSE")
        if df_gse_matched is not None:
            make_figure(df_gse_matched, "GSE", OUT_GSE)
    except Exception as e:
        import traceback
        print(f"⚠  GSE figure failed: {e}")
        traceback.print_exc()

    # --- GSM ---
    print(f"\nLoading GSM positions...")
    try:
        df_gsm_pos = load_positions_gsm()
        if df_gsm_pos is not None:
            df_gsm = match_positions(df_mach, df_gsm_pos, "GSM")
            if df_gsm is not None:
                make_figure(df_gsm, "GSM", OUT_GSM)
    except Exception as e:
        print(f"⚠  GSM figure failed: {e}")

    # --- GPE ---
    print(f"\nComputing GPE positions (aberration-corrected GSE)...")
    try:
        if df_gse_matched is None:
            df_gse_pos     = load_positions_gse()
            df_gse_matched = match_positions(df_mach, df_gse_pos, "GSE")
        if df_gse_matched is not None:
            df_gpe_raw  = gse_to_gpe(df_gse_matched)
            df_gpe_plot = df_gpe_raw.rename(columns={
                "X_RE": "X_RE_GSE", "Y_RE": "Y_RE_GSE", "Z_RE": "Z_RE_GSE",
                "X_RE_GPE": "X_RE",  "Y_RE_GPE": "Y_RE",  "Z_RE_GPE": "Z_RE",
            })
            make_figure(df_gpe_plot, "GPE", OUT_GPE)
    except Exception as e:
        print(f"⚠  GPE figure failed: {e}")

    # --- Aberration comparison ---
    if df_gse_matched is not None and df_gpe_raw is not None:
        print(f"\nGenerating aberration comparison figure...")
        try:
            plot_aberration_comparison(df_gse_matched, df_gpe_raw, OUT_ABER)
        except Exception as e:
            print(f"⚠  Aberration comparison failed: {e}")

    # --- Geotail-aware figures (saved) ---
    print(f"\nGeotail-aware figures ...")
    for label, df_coord, out_path in [
            ("GSE", df_gse_matched, OUT_GSE_GEOTAIL),
            ("GSM", df_gsm,         OUT_GSM_GEOTAIL),
            ("GPE", df_gpe_plot,    OUT_GPE_GEOTAIL),
    ]:
        if df_coord is None:
            print(f"  {label}: skipped (no data)")
            continue
        try:
            make_figure_geotail(df_coord, label, out_path)
        except Exception as e:
            import traceback
            print(f"⚠  Geotail figure {label} failed: {e}")
            traceback.print_exc()

    # --- Equal-count Mach bins (6 figures: scatter + geotail, × 3 coords) ---
    print(f"\nBuilding equal-count Mach bins ...")
    mach_vals = df_mach["alfven_mach"].dropna().to_numpy()
    eq_bins, eq_labels, eq_merka = build_equal_count_bins(mach_vals, n_bins=4)
    print(f"  Equal-count bin edges (4 bins):")
    for (lo, hi), lbl in zip(eq_bins, eq_labels):
        n = ((df_mach["alfven_mach"] >= lo) & (df_mach["alfven_mach"] <= hi)).sum()
        print(f"    {lbl}:  {n} crossings  (median MA = "
              f"{df_mach.loc[(df_mach['alfven_mach'] >= lo) & (df_mach['alfven_mach'] <= hi), 'alfven_mach'].median():.2f})")

    with _override_bin_config(eq_bins, eq_labels, eq_merka):
        for label, df_coord, out_sc, out_geo in [
                ("GSE", df_gse_matched, OUT_GSE_EQBIN,         OUT_GSE_GEOTAIL_EQBIN),
                ("GSM", df_gsm,         OUT_GSM_EQBIN,         OUT_GSM_GEOTAIL_EQBIN),
                ("GPE", df_gpe_plot,    OUT_GPE_EQBIN,         OUT_GPE_GEOTAIL_EQBIN),
        ]:
            if df_coord is None:
                print(f"  {label}: skipped (no data)")
                continue
            try:
                make_figure(df_coord, f"{label} — equal-count bins", out_sc)
                make_figure_geotail(df_coord, f"{label} — equal-count bins", out_geo)
                if label == "GPE":
                    make_figure_geotail_4panel(
                        df_coord, f"{label} — equal-count bins",
                        OUT_GPE_GEOTAIL_EQBIN_4PANEL)
            except Exception as e:
                import traceback
                print(f"⚠  Equal-bin figures {label} failed: {e}")
                traceback.print_exc()

    # --- Bin-size sensitivity (inline only, X-R panel, all three coords) ---
    print(f"\nBin-size sensitivity plots (inline) ...")
    for label, df_coord in [("GSE", df_gse_matched),
                             ("GSM", df_gsm),
                             ("GPE", df_gpe_plot)]:
        if df_coord is None:
            print(f"  {label}: skipped (no data)")
            continue
        try:
            make_figure_binsize_test(df_coord, label)
        except Exception as e:
            print(f"⚠  Bin-size test {label} failed: {e}")

    print("\n" + "=" * 65)
    print("Done.")
    print("=" * 65)


if __name__ == "__main__":
    main()
