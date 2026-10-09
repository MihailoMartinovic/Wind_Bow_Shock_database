# The Wind bow shock crossing list

This document describes `database_files/bs_crossings_V02.xlsx` and its HDF5 copy: what each column holds,
how the rows are organised, and what changed in release 1.2.1. The renderings built from the list are
described in [README.md](README.md).

## Summary

| | |
|---|---|
| Crossings | 644: 322 shock entries, 322 shock exits |
| Passes | 70, numbered 1–70 in time order |
| Pass legs | 140: an inbound and an outbound leg per pass, holding 354 and 290 crossings |
| Crossings per leg | 1 to 65, median 3; 59 legs hold a single crossing |
| Time span | 1995-08-01 21:16:41 to 2004-04-27 04:48:23 UT |
| Files | `database_files/bs_crossings_V02.xlsx`, sheet `bs_crossings_V02`; `database_files/Wind_bow_shock_database.h5`, key `all` |

## The file

The sheet `bs_crossings_V02` holds a header row and 644 data rows in columns A–H. The rows of each leg are
contiguous, but the file is not in pass order (it opens with pass 6), so group by `pass` and `SC direction`,
or sort by `pass` and `time`, rather than relying on row order.

The workbook's second sheet, `Sheet1`, is a leftover scratch copy of 22 rows that also appear in the main
sheet. It is not part of the database.

## Columns

| Column | Meaning | Values and format |
|---|---|---|
| `pass` | Pass number | Integer, 1–70 |
| `SC direction` | Leg of the pass | `in` (inbound) or `out` (outbound) |
| `time` | Crossing time | UT, text `YYYY-MM-DD HH:MM:SS.f`, given to 0.1 s; 69 values have a non-zero tenth |
| `direction` | Crossing direction | `in`: shock entry, from the solar wind into the magnetosheath. `out`: shock exit |
| `start time`, `end time` | Broad window of the leg | UT, text `YYYY-MM-DD HH:MM:SS` |
| `start time zoom`, `end time zoom` | Zoom window of the leg | UT, text `YYYY-MM-DD HH:MM:SS` |

The accompanying paper calls `SC direction` the *pass direction* and `direction` the *crossing direction*.
The list holds times only; positions come from the *Wind* ephemeris.

## The windows

The four window columns describe the leg, not the crossing, and since release 1.2.1 they carry the same
values on every row of the leg.

- The **broad window** is the time span of the data-product figure `figs_data/pass_NN_<leg>.png`. Lengths run
  from 0.68 h to 69 h (median 6.5 h). The median padding is 3.06 h before the leg's first crossing and 3.13 h
  after its last.
- The **zoom window** is the time span of `figs_data/pass_NN_<leg>_zoom.png`. Lengths run from 2 minutes to
  57.5 h (median 40 minutes); the longest belongs to pass 45 inbound, whose 17 crossings span more than two
  days.
- The **midpoint of the zoom window** is the reference time t0 of the leg's 3D rendering. It picks the OMNI2
  hour used for the model boundaries and gives the rendering its file name.

Every crossing lies inside its leg's broad and zoom windows.

Within a leg, entries and exits alternate in time. An inbound leg starts and ends with an entry and an
outbound leg starts and ends with an exit, so with 70 legs of each kind the list holds as many entries as
exits.

Before release 1.2.1 the window values were given only on the first row of each leg, and several notebooks
found the legs by looking for those filled cells. Those notebooks (`render_passes_v2.ipynb`,
`build_omni2_cache.ipynb`, `remake_figure4.ipynb`) now take a leg to be a run of rows with the same `pass` and
`SC direction`. On the old file both rules give the same 140 legs.

## Changes in release 1.2.1 (8 October 2026)

The list changed in two steps. The first completed the window columns and the directions of pass 8
inbound. The second corrected six legs and brought the list from 635 to 644 crossings. Cell and row
numbers in step 1 refer to the file as it was before step 2.

### Step 1: windows and directions (1,993 cells)

`Sheet1` and the `pass`, `SC direction` and `time` columns are unchanged in this step.

**Read from the data-product figures (13 cells).** Each figure's time axis was measured against its own
gridlines, and its dashed crossing lines against the crossing times in the list.

| Cells | Leg | Before → after | Evidence |
|---|---|---|---|
| D48–D52 | 8 in | blank → `in`, `out`, `in`, `out`, `in` | `pass_08_in_zoom.png`: five lines within 2 s of the listed times, black (entry) and brown (exit) alternating. \|B\| rises from about 5 to 20 nT at each entry and falls back at each exit |
| E599 | 68 in | 2003-08-12 16:30:00 → 13:00:00 | `pass_68_in.png`: all 29 crossing lines fit the axis at 41.3 s per pixel (11 s rms). The eight 2-hour gridlines fall within 23 s of 13:00–05:00. The old value was the zoom start entered in the broad column |
| F599 | 68 in | blank → 2003-08-13 05:00:00 | as above |
| E10, F10 | 5 in | 1995-12-22 11:50:00–15:30:00 → 1995-12-20 20:30:00–22:30:00 | `pass_05_in.png`: axis labelled 20:30 and 22:30 at its edges, and 15-minute gridlines within 3 s. The old values were pass 5 outbound's window and contained none of this leg's seven crossings |
| G2, H2 | 6 in | 1996-01-12 15:05:51–15:06:10 → 15:05:00–15:07:00 | `pass_06_in_zoom.png`: all three crossing lines exact, 15-second gridlines exact. The old values were the first two crossing times and left out the third |
| G74, H74 | 12 in | 1996-10-05 00:50:00–01:00:00 → 1996-10-04 00:50:00–01:00:00 | `pass_12_in_zoom.png`: date label 1996 Oct 04, minute gridlines within 0.2 s, and the line at 00:54:55 matches the crossing. The old date was a day after the crossing and outside the broad window |

**Repeated down each leg (1,980 cells in 495 rows).** Each row received its leg's four window values,
copied unchanged from the row that held them.

Two renderings depended on the corrected cells and are redone in 1.2.1. Pass 8 inbound had no crossing
markers, because none of its crossings had a direction. Pass 12 inbound was centred on 1996-10-05 00:55, a
day after its crossing, and used that day's OMNI2 values; it becomes
`Magnetosphere_Pass_12_Inbound_1996-10-04T0055.html`. The zoom midpoint of pass 6 inbound moves by half a
second, which changes neither its file name nor its OMNI2 hour.

### Step 2: six legs corrected (644 crossings)

In these six legs the crossings did not alternate as described in [The windows](#the-windows), or the
data-product figures drew crossings that the list lacked. Each change was checked against the \|B\| and
density panels of the leg's data-product figures, and against the one-minute *Wind* list of Merka et al.
that the accompanying paper compares with. The dashed lines in the data-product figures were drawn from
times typed into the plotting notebook. Their colour alternates from line to line, entry first on inbound
legs and exit first on outbound legs, so on its own it is not evidence of a crossing's direction.

| Leg | Rows now | Change | Evidence |
|---|---|---|---|
| 7 in | 40 | 1996-03-28 14:44:29.5 → 1996-03-26 16:10:29.0, entry | The old time lay two days after the leg, outside both its windows, and repeated the outbound exit at 14:44:53.5. Both inbound figures draw an entry at 16:10:28 ± 2 s, where \|B\| rises from about 6 to 25 nT. Merka et al.: 16:10, entry |
| 18 in | — | Row 109 deleted | Identical to row 107 (15:25:06.0, exit) |
| 18 in | 109 | 15:55:35.0: entry → exit | \|B\| falls from about 12 to 3 nT and the density from about 30 to 10 cm⁻³, and the solar wind lasts until 16:42. Merka et al.: 15:55, exit |
| 18 in | 110 | Added: 1997-07-25 16:42:17.0, entry | \|B\| jumps at 16:42:17 (read from the zoom figure, ±3 s) and stays at magnetosheath level from 16:43:03 until the exit at 16:48:59.5. Merka et al.: 16:43 |
| 20 in | 116–118 | 09:43:18.0: exit → entry. 09:51:40.0: entry → exit. 10:15:10.0 (entry) moved from the first row of the leg to the last | At 09:43:18 \|B\| rises from about 5 to 20 nT and the flow speed falls from 420 to 310 km/s; both return at 09:51:40. Merka et al.: 09:43 entry, 09:51 exit |
| 33 in | 207–210 | Added: 1999-03-11 05:13:55.0 exit, 05:36:58.0 entry, 05:40:45.0 exit, 07:56:59.0 entry | Drawn in both figures after the last listed crossing (05:05:54). Times read from the zoom figure, ±10 s. Merka et al.: crossing pairs from 05:06 and from 05:37, and an entry at 07:57 |
| 47 in | 441–444 | Added: 2000-05-02 12:55:12.0 exit, 13:03:12.0 entry, 13:09:52.0 exit, 14:45:10.0 entry | Drawn in both figures after the last listed crossing (12:46:58). Times read from the zoom figure, ±20 s. Merka et al.: 12:55 exit, 13:03 entry, 13:10 exit, 14:45 entry |
| 53 out | 520–528 | Directions reversed for the nine crossings from 12:17:51 to 13:33:28 | Rows 511 and 512 of the old file were consecutive exits. With the reversed directions, \|B\| is at solar-wind level (5–7 nT) after each exit and above 10 nT after each entry; the intervals between 12:22 and 12:31 are only one to three pixels wide in the zoom figure. Merka et al.: 12:17 entry, 12:33 exit, 13:33 entry |
| 53 out | 529 | Added: 2000-08-06 15:06:39.0, exit | Drawn in both figures. \|B\| falls from about 15 to 5 nT, and without this exit the leg would end with an entry. Time ±24 s. Merka et al.: 15:06, exit |

Rows below row 206 of the old file moved down by 4, below row 436 by 8 and below row 520 by 9.

The renderings of these six legs mark the crossings and are redone. Their file names do not change,
because the windows do not. The data-product figures of passes 18 and 20 inbound (`pass_18_in.png`,
`pass_18_in_zoom.png`, `pass_20_in.png`, `pass_20_in_zoom.png`) are redrawn from the corrected times; the
figures of the other four legs already showed them. README.md lists every data-product figure redrawn for
release 1.2.1.

## Open entries

**Two crossing times** differ slightly from the single line in their zoom figures.

| Leg | Row | What the zoom figure shows |
|---|---|---|
| 17 out | 105 | The line is 23 s after the listed time (21:41:45) |
| 22 out | 139 | The line is 10 s after the listed time (13:00:00.5) |

**Tail crossings that may not be bow shock crossings.** The 21 crossings below lie 10 to 44 R⊕ inside the
model bow shock evaluated at each crossing (`distance` in `database_files/crossing_distances_v2.csv`), and
most of them inside or at the model magnetopause (Shue et al. 1998).
None is in the one-minute list of Merka et al. Several zoom figures show a steady field along ±X on one side
of the line, as in a tail lobe; pass 33 outbound, for example, has about 25 nT along +X. Pass 55 outbound has
no magnetic field data. They may be magnetopause or lobe boundaries rather than bow shock crossings. They are
all the crossings of their 11 legs and are kept in the list, and comments are welcome. Figure 5 of the
accompanying paper flags them as tail crossings and leaves them out of its fits.

| Leg | Rows | Times (UT) | x, ρ (GSM, R⊕) | Inside the model bow shock by (R⊕) |
|---|---|---|---|---|
| 30 out | 183 | 1999-01-19 10:00:57 | −14, 20 | 13 |
| 32 out | 197 | 1999-02-23 18:52:05 | −24, 18 | 20 |
| 33 out | 211–213 | 1999-03-12 22:55:30, 22:58:20, 23:11:29 | −26, 14 | 22 |
| 34 out | 245 | 1999-04-02 21:18:07 | −59, 23 | 28 |
| 44 out | 399–405 | 1999-12-20 09:22:42 to 23:31:57, seven crossings | −49, 24 | 24–25 |
| 54 out | 531–533 | 2000-08-17 15:42:30, 15:43:10, 15:44:00 | −31, 25 | 16 |
| 55 out | 535 | 2001-02-24 22:56:52 | −91, 19 | 42 |
| 58 out | 561 | 2001-10-24 20:43:06 | −14, 25 | 10 |
| 61 in | 578 | 2002-07-20 15:14:46 | −21, 23 | 15 |
| 68 out | 637 | 2003-08-19 17:43:37 | −80, 37 | 15 |
| 69 out | 639 | 2004-02-20 18:55:11 | −119, 20 | 44 |

## The HDF5 copy

`database_files/Wind_bow_shock_database.h5` holds the same eight columns and 644 rows as a pandas table under
the key `all`. It is generated from the list by `scripts/build_hdf5_from_crossings.ipynb`, which keeps the
previous file's key, column types and storage format, and reads the result back to check it. Read it with
`pd.read_hdf(path, key="all")` (needs PyTables).

## Reading the list

```python
import pandas as pd

x = pd.read_excel("database_files/bs_crossings_V02.xlsx", sheet_name="bs_crossings_V02")
for c in ["time", "start time", "end time", "start time zoom", "end time zoom"]:
    x[c] = pd.to_datetime(x[c])

legs = x.groupby(["pass", "SC direction"], sort=False)            # 140 legs
w = legs[["start time zoom", "end time zoom"]].first()
t0 = w["start time zoom"] + (w["end time zoom"] - w["start time zoom"]) / 2   # reference time of each rendering
```
