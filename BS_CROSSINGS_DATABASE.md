# Wind Bow Shock Crossings Database

## Overview

This document describes the Wind spacecraft bow shock crossing database and the methodology for generating interactive magnetosphere visualizations.

**Database Summary:**
- **635 individual crossings** across **70 passes** of the Wind spacecraft
- **Temporal Range:** August 1995 to April 2004 (~9 years)
- **Pass Legs:** 140 total (70 passes × 2 directions), with 137 rendered and 3 skipped

A **pass** is the portion of the Wind trajectory traversing the bow shock region, classified as **inbound** (entering magnetosphere) or **outbound** (exiting). A **crossing** is a discrete event where Wind crosses the boundary.

---

## Database Structure

### Crossing Database File
- **Location:** `database_files/`
- **Formats:** Excel (.xlsx) and HDF5 (.h5)
- **Rows:** 635 crossings, one per row
- **Columns:** Pass ID, direction, crossing time (UTC), crossing properties

### Provenance Record
**File:** `pass_parameters_v2.csv`

One row per pass leg (140 total, including 3 skipped):

```
pass_id, sc_direction, t0, n_crossings, status, omni_epoch, offset_h, 
Pdyn, Bz, By, Dst, Ma, r0_mp, alpha, delta, r_bs0, clamped, n_field_lines, seconds, html
```

**Key columns:**
- `status`: "rendered" (137 legs) or "skipped: no complete OMNI hour within 4.0 h" (3 legs)
- `t0`: UTC timestamp of first crossing in the leg
- `n_crossings`: Number of crossings in this leg
- `omni_epoch`: OMNI2 hour used (may be offset from t0)
- `offset_h`: Time offset from t0 to nearest OMNI hour (±4 h max)
- Model parameters: `Pdyn`, `Bz`, `By`, `Dst`, `Ma` (Alfvén Mach number)
- Boundary distances: `r0_mp` (magnetopause nose), `delta` (standoff), `r_bs0` (shock nose)
- Computation: `n_field_lines` (96), `seconds` (runtime), `html` (output filename)

---

## How Surfaces Are Computed

### Magnetopause (Shue et al. 1998)

The magnetopause subsolar distance depends on solar wind dynamic pressure and IMF Bz:

$$r_0 = \left(10.22 + 1.29 \tanh(0.184(B_z + 8.14))\right) P_{dyn}^{-1/6.6}$$

$$\alpha = (0.58 - 0.007 B_z)(1 + 0.024 \ln P_{dyn})$$

$$r(\theta) = r_0 \left(\frac{2}{1 + \cos\theta}\right)^\alpha$$

**Implementation:** Shue et al. Equations 10 and 11, verified line-by-line against the original paper.

**Surface rendering:** 100 azimuths at opacity 0.35 (semi-transparent mesh).

### Bow Shock (Gas-Dynamic Standoff)

The bow shock is a conic section of eccentricity *e* = 0.9, with subsolar distance equal to the magnetopause nose plus the gas-dynamic standoff:

$$\Delta = r_0 \frac{(\gamma - 1)M_A^2 + 2}{(\gamma + 1)(M_A^2 - 1)}, \quad \gamma = \frac{5}{3}$$

$$r_{bs0} = r_0 + \Delta$$

$$r(\theta) = r_{bs0} \frac{1 + e}{1 + e\cos\theta}$$

**Critical note on Mach number:** This relation uses the **Alfvén Mach number** (M_A), not the magnetosonic Mach number. Farris and Russell (1994) published this relation with M_ms; at β = 1, M_ms ≈ 0.74 M_A. Above M_A ≈ 4, the two choices differ by under 1 R⊕, so figures are largely insensitive—but documentation must specify which Mach number is used.

**Surface rendering:** 80 azimuths at opacity 0.20 (semi-transparent surface).

### Magnetic Field Lines

- **Model:** T96 (Tsyganenko 1996) with IGRF magnetic field coefficients
- **Seeds:** 96 points uniformly distributed in longitude
- **Direction:** Traced in both directions from each seed
- **Limit:** rlim = 60 R⊕
- **Purpose:** Visualization only; boundaries are not derived from field lines

---

## Upstream Parameter Selection

**Data source:** OMNI2 hourly data from `D:\Data\OMNI\omni2_cache.pkl`

**Selection rule:** For each pass leg, the generator finds the OMNI2 hour with the **nearest timestamp to the first crossing** at which all five required parameters are present:

1. Dynamic pressure (Pdyn)
2. IMF Bz component
3. IMF By component
4. Disturbance Storm Time (Dst)
5. Alfvén Mach number (Ma)

**Constraint:** The selected hour must lie within ±4 hours of the first crossing.

**Distribution (137 rendered legs):**
- 132 legs: OMNI hour within ±1 h of crossing
- 2 legs: 1–2 h offset
- 3 legs: 2–4 h offset (largest: 3.19 h)

### Parameter Ranges

| Parameter | Min | Median | Max |
|---|---|---|---|
| M_A | 1.40 | 8.10 | 39.70 |
| P_dyn (nPa) | 0.32 | 1.96 | 9.23 |
| r₀ magnetopause (R⊕) | 7.42 | 10.23 | 13.62 |
| α (exponent) | 0.54 | 0.59 | 0.71 |
| Δ standoff (R⊕) | 1.91 | 2.70 | 16.85* |
| r_bs0 shock nose (R⊕) | 9.45 | 13.00 | 29.89* |

*Pass 47 inbound is an outlier; see below.

---

## Skipped Legs

Three pass legs have no OMNI2 hour within ±4 h of the first crossing containing all five required parameters:

| Leg | First Crossing | Crossings Lost |
|---|---|---|
| Pass 4 inbound | 1995-11-28 08:40 | 3 |
| Pass 55 inbound | 2001-02-20 02:45 | 1 |
| Pass 56 inbound | 2001-09-08 12:05 | 7 |

**Total:** 11 of 635 crossings are not rendered (624 rendered, 11 lost to data gaps).

These legs appear in `pass_parameters_v2.csv` with empty model columns, so the absence is recorded rather than silent.

---

## Known Issues

### Pass 47 Inbound — Low Alfvén Mach Number

This leg has M_A = 1.40 at an hour 15 minutes from the crossing, with P_dyn = 0.32 nPa. The gas-dynamic relation diverges as M_A → 1:

$$\Delta = 16.85 \text{ R}_\oplus, \quad r_{bs0} = 29.9 \text{ R}_\oplus$$

compared to the typical 10–16 R⊕ for shock nose across all 137 legs. It is the only leg exceeding 10 R⊕; 136 of 137 have Δ ≤ 5 R⊕.

Pass 47 is already noted in `00_pass_notes.txt` as "problematic", making this divergence explicit. The leg is rendered and included, with the understanding that the model surfaces are not physically meaningful at M_A = 1.4. The rendering is valuable for showing the spacecraft trajectory and crossing times; readers should disregard the boundary surfaces for this leg.

---

## Changes from v1.0.0 to v1.1.0

### Parameter Source (Breaking Change)

**v1.0.0** fell back to Wind's own MFI and SWE measurements when the OMNI2 cache was missing. Wind spends part of every pass in the magnetosheath (inside the bow shock), so these parameters were not upstream solar wind values.

**v1.1.0** uses OMNI2 hourly data exclusively, with no fallback. This:
- Removed 18 files labeled `_WindData` (passes 3o, 4i, 5i, 5o, 6o, 8i, 9i, 13o, 15o, 16o, 19o, 20o, 22i, 24o, 31i, 42o, 55i, 56i)
- Fixed severe errors in seven of those passes
  - **Pass 9 inbound:** previously r₀ = 78 R⊕ (P_dyn fallback = 0), now r₀ = 9.34 R⊕ (OMNI P_dyn = 3.55 nPa)
  - **Pass 19 outbound:** Δ now 3.32 R⊕ (was clipped at 10)
  - **Pass 24 outbound:** Δ now 4.99 R⊕ (was clipped at 10)

### Standoff Distance Clamping (Removed)

**v1.0.0** enforced 1 ≤ Δ ≤ 10 R⊕ silently. The `clamped` column in v1.1.0's CSV is `False` for all 137 legs—the clamp is gone.

This exposes pass 47 inbound's diverging behavior (see above), which was previously hidden by the 10 R⊕ cap.

### Filename Corrections

Three old files carried the wrong timestamp:

| Leg | v1.0.0 stamp | Correct t0 |
|---|---|---|
| Pass 3 outbound | 1995-09-16T1017 (pass 3 in time) | 1995-09-18T0420 |
| Pass 34 inbound | 1999-03-27T2235 | 1999-03-28T0915 |
| Pass 53 outbound | 2000-08-05T2154 | 2000-08-06T0637 |

v1.1.0 filenames use correct crossing times from the crossing list.

---

## Rendering Pipeline

**Generator:** `render_passes_v2.ipynb` (resumable batch notebook)

**Inputs:**
- `pass_parameters_v2.csv` (or regenerated from crossing list)
- OMNI2 cache: `D:\Data\OMNI\omni2_cache.pkl` (built by `build_omni2_cache.ipynb`)

**Process:**
1. For each pass leg:
   - Check if HTML already exists (if `OVERWRITE = False`, skip)
   - Query OMNI2 cache for parameters within ±4 h of first crossing
   - If incomplete, log as skipped
   - If complete, compute Shue and standoff, trace field lines
   - Render Plotly 3D scene with surfaces and trajectory
   - Save as interactive HTML
2. Log all runs to `pass_parameters_v2.csv`

**Runtime:** ~3 h 55 min total (median 104 s per leg).

**Resume:** Interrupt the notebook at any time. Re-run to skip completed legs and continue from the last incomplete one.

---

## References

- Farris, M. H., and C. T. Russell (1994), Bow shock and magnetosheath dynamic pressure variations relative to the solar wind dynamic pressure, *J. Geophys. Res.*, 99, 17681–17689.
- Shue, J. H., et al. (1998), Magnetopause location under extreme solar wind conditions, *J. Geophys. Res.*, 103, 17691–17700.
- Tsyganenko, N. A. (1996), Modeling the Earth's magnetospheric magnetic field confined within a realistic magnetopause, *J. Geophys. Res.*, 101, 27187–27198.

---

**Last updated:** October 2026 (v1.1.0)

For usage examples and installation, see **README.md**.
