# Wind Bow Shock Database

A list of Earth's bow shock crossings identified in *Wind* spacecraft data, with data-product figures for
every pass and interactive 3D renderings of the model magnetopause and bow shock along the spacecraft
trajectory.

- **644 crossings** in **140 pass legs** (70 passes, each with an inbound and an outbound leg)
- **1 August 1995 to 27 April 2004**
- https://github.com/MihailoMartinovic/Wind_Bow_Shock_database

**Terms.** A *pass* is one excursion of *Wind* through the region of the bow shock; passes are numbered 1–70
in time order. Each pass has an *inbound* leg (toward Earth) and an *outbound* leg. A *crossing* is one
traversal of the shock, either an *entry* (`in`, from the solar wind into the magnetosheath) or an *exit*
(`out`). A leg holds between 1 and 65 crossings; 59 legs hold one.

## Contents

| Path | What it is |
|---|---|
| `database_files/bs_crossings_V02.xlsx` | The crossing list. Column dictionary and conventions in [BS_CROSSINGS_DATABASE.md](BS_CROSSINGS_DATABASE.md) |
| `database_files/Wind_bow_shock_database.h5` | The same table as a pandas HDF5 file under the key `all`, written from the list by `scripts/build_hdf5_from_crossings.ipynb` |
| `database_files/pass_parameters_v2.csv` | Provenance of every rendering in `figs/`: one row per leg with the upstream parameters, the model results and the file name |
| `database_files/crossing_distances_v2.csv` | Per-crossing position and distance from the model bow shock under both distance definitions, from which Figure 4 of the paper is drawn |
| `database_files/crossing_positions_gpe_v3.csv` | Per-crossing GSE and GPE positions (as computed and pressure-normalised), upstream parameters and the tail flag `tail_fixed` (22 crossings), for the 633 crossings with an OMNI2 Alfvén Mach number within 4 h |
| `database_files/mach_per_crossing_v2.csv` | Per-crossing upstream parameters (Alfvén Mach number, solar wind speed, proton density, dynamic pressure) and the time offset of the OMNI2 value used |
| `database_files/trajectory_*.html` | Interactive plots of the *Wind* trajectory with all crossings |
| `figs/` | 137 interactive 3D renderings (HTML), one per rendered leg |
| `figs_data/` | 280 data-product figures (PNG), two per leg: full window and zoom |
| `figs_deprecated/` | The 140 renderings of release 1.0.0, kept for comparison |
| `MagCarto_Master.ipynb`, `bs_crossings_loader.py` | Interactive notebook from the original pipeline, and the loader it imports from the same folder |
| `scripts/` | The notebooks that produce the renderings, the HDF5 file and Figures 4 and 5 of the paper (listed below), plus the plotting and movie scripts of the original pipeline |

## The crossing list

`database_files/bs_crossings_V02.xlsx`, sheet `bs_crossings_V02`, one row per crossing:

| Column | Content |
|---|---|
| `pass` | Pass number, 1–70 |
| `SC direction` | The leg: `in` (inbound) or `out` (outbound) |
| `time` | Crossing time, UT |
| `direction` | `in` = shock entry, `out` = shock exit |
| `start time`, `end time` | The leg's broad window: the time span of `figs_data/pass_NN_<leg>.png` |
| `start time zoom`, `end time zoom` | The leg's zoom window: the time span of `figs_data/pass_NN_<leg>_zoom.png`. Its midpoint is the reference time of the leg's 3D rendering |

The four window columns describe the leg, not the crossing, and carry the same values on every row of a leg.
Times are stored as text. Units, formats, statistics and the changes made in release 1.2.1 are in
[BS_CROSSINGS_DATABASE.md](BS_CROSSINGS_DATABASE.md).

## The 3D renderings (`figs/`)

Produced by `scripts/render_passes_v2.ipynb`, one file per leg:

1. **Reference time.** t0 is the midpoint of the leg's zoom window. The file name is
   `Magnetosphere_Pass_NN_<Inbound|Outbound>_<t0 as YYYY-MM-DDTHHMM>.html`.
2. **Upstream parameters.** OMNI2 hourly dynamic pressure, B<sub>z</sub> and B<sub>y</sub> (GSM), Dst and
   Alfvén Mach number M<sub>A</sub>, taken from the hour nearest t0 at which all five are present, provided it
   lies within 4 h of t0. No other source is used.
3. **Magnetopause.** Shue et al. (1998), eqs. 10 and 11:
   ```
   r0    = (10.22 + 1.29 tanh(0.184 (Bz + 8.14))) Pdyn^(-1/6.6)
   alpha = (0.58 - 0.007 Bz) (1 + 0.024 ln Pdyn)
   r(θ)  = r0 (2 / (1 + cos θ))^alpha
   ```
4. **Bow shock.** A conic of eccentricity 0.9, `r(θ) = r_bs0 (1 + e) / (1 + e cos θ)`, which closes far down
   the tail (it is an ellipse, not a paraboloid). Its subsolar distance is the magnetopause nose plus a
   gas-dynamic standoff evaluated with the **Alfvén** Mach number, with no limits applied:
   ```
   Δ     = r0 ((γ - 1) M_A² + 2) / ((γ + 1) (M_A² - 1)),   γ = 5/3
   r_bs0 = r0 + Δ
   ```
   Farris and Russell (1994) write this relation with the magnetosonic Mach number. A leg with
   M<sub>A</sub> ≤ 1.2 would be skipped; none is.
5. **Field lines.** T96 with IGRF through `geopack`, traced in both directions from 96 seeds (a shell of 72 at
   2.5 R⊕ and two tailward fans of 12 at 3.5 and 5 R⊕) out to 60 R⊕. They are drawn only; neither boundary is
   derived from them.
6. **Trajectory.** *Wind* GSM positions at 10-minute cadence for t0 ± 12 h. Entries and exits are marked
   separately at the nearest trajectory point.
7. **Symmetry.** Both boundaries are surfaces of revolution about the GSM x axis: no aberration, no dawn–dusk
   asymmetry and no dipole tilt.

Each figure's title gives t0, and its annotation gives the OMNI2 hour used and its offset from t0, the five
parameters, r0, α, the standoff Δ, the subsolar distance r<sub>bs0</sub>, e and the number of field lines.
The same values, for all 140 legs, are in `database_files/pass_parameters_v2.csv`.

**Coverage.** 137 legs are rendered, holding 633 of the 644 crossings. Passes 4, 55 and 56 inbound
(11 crossings) have no hour within 4 h of t0 at which all five OMNI2 fields are present; they appear in the
parameter table with `status = "skipped: no complete OMNI hour within 4.0 h"`. For 132 of the 137 rendered
legs the OMNI2 hour lies within 1 h of t0; the largest offset is 3.19 h.

| Over the 137 rendered legs | min | median | max |
|---|---|---|---|
| M<sub>A</sub> | 1.40 | 8.10 | 39.70 |
| P<sub>dyn</sub> (nPa) | 0.32 | 1.95 | 9.23 |
| Magnetopause r0 (R⊕) | 7.42 | 10.25 | 13.62 |
| α | 0.54 | 0.59 | 0.71 |
| Standoff Δ (R⊕) | 1.91 | 2.71 | 16.85 |
| Bow shock nose r<sub>bs0</sub> (R⊕) | 9.45 | 13.00 | 29.89 |

**Inputs that are not in the repository.** The OMNI2 cache `omni2_cache.pkl`, built by
`scripts/build_omni2_cache.ipynb` from the OMNI2 database that spacepy keeps in
`~/.spacepy/data/omni2data.h5`, and the *Wind* ephemeris `Wind_Ephemerids.h5` (key `/GSM2`). The notebooks, and
`MagCarto_Master.ipynb`, read these locations from a `local_paths.json` in the folder they run in. That file is not
part of the repository; `scripts/local_paths.example.json` is the template.

## The data-product figures (`figs_data/`)

Two figures per leg: `pass_NN_<leg>.png` over the broad window and `pass_NN_<leg>_zoom.png` over the zoom
window. Seven panels share the time axis: a radio spectrogram
(V², dB) against frequency; |B| and its GSE components; the angle between **B** and the GSE x axis; electron
densities (coarse and fit) and proton density; solar wind speed; temperatures (60 s medians); and the quality
factor QF. Crossings are dashed vertical lines, black for a shock entry and brown for a shock exit.

The figures were made between February and May 2026. For release 1.2.1, 30 of them were redrawn from the
crossing list:
- `pass_05_in`, `pass_08_in`, `pass_17_in_zoom`, `pass_24_out`, `pass_47_out`, `pass_49_out`, `pass_53_in`
  and `pass_64_out_zoom`, whose lines or time axes did not match the list;
- both figures of passes 18 and 20 inbound, whose crossings were corrected;
- the 18 figures of passes 46, 49, 50, 52 and 54 that had been exported as JPG at lower resolution.

The lines of every figure now match the list, except two that differ from it by 10 and 23 s (see
[Known issues](#known-issues)).

## Scripts and notebooks (`scripts/`)

| File | Produces |
|---|---|
| `render_passes_v2.ipynb` | The renderings in `figs/` and `pass_parameters_v2.csv`. Can re-render a chosen list of legs (`ONLY_LEGS`) |
| `build_omni2_cache.ipynb` | The OMNI2 cache the rendering and the Figure 4 and 5 notebooks read |
| `build_hdf5_from_crossings.ipynb` | `Wind_bow_shock_database.h5` from the crossing list, with the previous file's key, columns and types |
| `remake_figure4.ipynb` | Figure 4 (distance between the crossings and the model bow shock) and `crossing_distances_v2.csv` |
| `remake_figure5.ipynb`, `plot_bs_positions_by_mach.py` | Figure 5 (crossing positions binned by Alfvén Mach number, with the Merka et al. (2005) bow shock surface and fits of it to each Mach bin) |

The remaining scripts belong to the original pipeline (`MagCarto_Master.ipynb`) and are not used for the
1.2.x renderings.

After any change to the crossing list, run `build_hdf5_from_crossings.ipynb`, then `render_passes_v2.ipynb`
for the affected legs, then `remake_figure4.ipynb` and `remake_figure5.ipynb`.

## Installation and use

```bash
git clone https://github.com/MihailoMartinovic/Wind_Bow_Shock_database.git
cd Wind_Bow_Shock_database
python -m venv venv-magnetosphere
venv-magnetosphere\Scripts\activate          # Windows; on macOS/Linux: source venv-magnetosphere/bin/activate
pip install -r requirements.txt
```

Reading the list needs only pandas and openpyxl:

```python
import pandas as pd

x = pd.read_excel("database_files/bs_crossings_V02.xlsx", sheet_name="bs_crossings_V02")
for c in ["time", "start time", "end time", "start time zoom", "end time zoom"]:
    x[c] = pd.to_datetime(x[c])

legs = x.groupby(["pass", "SC direction"])      # 140 legs
entries = x[x["direction"] == "in"]              # 322 shock entries
```

The HDF5 copy needs PyTables (`pip install tables`): `pd.read_hdf("database_files/Wind_bow_shock_database.h5", key="all")`.

`MagCarto_Master.ipynb` is interactive (ipywidgets). It fetches positions from SSCWeb and upstream conditions
through `spacepy.omni`, so it needs network access and spacepy's OMNI2 data
(`python -c "import spacepy.toolbox as tb; tb.update(omni2=True)"`).

## Known issues

- **Pass 47 inbound.** M<sub>A</sub> = 1.40 and P<sub>dyn</sub> = 0.32 nPa, from the OMNI2 hour at 07:30,
  15 minutes before t0. This is a real low-Mach interval, not a data gap. The gas-dynamic relation diverges as M<sub>A</sub> → 1 and gives Δ = 16.85 R⊕ and a shock nose at 29.89 R⊕,
  against 1.91–4.99 R⊕ for every other leg. The model surfaces in that rendering are not meaningful; the
  rendering itself carries no warning.
- **Three legs are not rendered**: passes 4, 55 and 56 inbound (see Coverage above).
- **Axisymmetric model surfaces**, as described above.
- **Two crossing times to confirm.** Pass 17 outbound (21:41:45) and pass 22 outbound (13:00:00.5) differ
  by 23 s and 10 s from the line in their zoom figures. Details in
  [BS_CROSSINGS_DATABASE.md](BS_CROSSINGS_DATABASE.md#open-entries).
- **Tail crossings that may not be bow shock crossings.** 21 crossings lie 10 to 44 R⊕ inside the model bow
  shock, most of them inside or at the model magnetopause. They are all the crossings of passes 30, 32, 33, 34,
  44, 54, 55, 58, 68 and 69 outbound and 61 inbound. None is in the list of Merka et al., and several zoom
  figures show a steady field along ±X on one side, as in a tail lobe. They may be magnetopause or lobe
  boundaries. They are kept in the list, and comments are welcome; details in
  [BS_CROSSINGS_DATABASE.md](BS_CROSSINGS_DATABASE.md#open-entries). Figure 5 of the paper flags them as
  tail crossings, with pass 53 outbound at 21:54:51, and leaves all 22 out of its fits.
- **Trajectory plots.** The `database_files/trajectory_*.html` files were made on 30 May 2026 from the list
  as it was then, and do not include the corrections of release 1.2.1.
- **Loader field names.** `bs_crossings_loader.py` stores column B (`SC direction`) in `crossing_direction`
  and column D (`direction`) in `sc_direction`. `MagCarto_Master.ipynb` uses them consistently, but code
  calling the loader directly should note the swap.

## Releases

| Version | Date | Notes |
|---|---|---|
| 1.0.0 | 22 May 2026 | First release. 18 legs used *Wind*'s own MFI/SWE parameters where OMNI was missing; those renderings, with the rest of that release, are in `figs_deprecated/` |
| 1.0.2 | 30 May 2026 | One pass added |
| 1.2.0 | 2 October 2026 | Renderings from OMNI2 only, standoff unclamped, three file timestamps corrected. The `v1.2.0` tag points at the commit before the release commit, so the tagged snapshot does not contain these renderings; use 1.2.1 |
| 1.2.1 | 8 October 2026 | Crossing list completed and corrected: directions for pass 8 inbound; windows of passes 5, 6, 12 and 68 inbound read from the data-product figures; window values repeated on every row; six legs corrected (passes 7, 18, 20, 33 and 47 inbound, 53 outbound), giving 644 crossings. HDF5 file regenerated, renderings of those eight legs redone, 30 data-product figures redrawn (the 18 JPGs among them now PNG), generator notebooks and analysis code added, documentation corrected |

## Data sources

*Wind* data from NASA CDAWeb; OMNI2 hourly data from NASA SPDF (through spacepy); *Wind* ephemeris from SSCWeb.

## References

- Farris, M. H., and C. T. Russell (1994), Determining the standoff distance of the bow shock: Mach number
  dependence and use of models, *J. Geophys. Res.*, 99(A9), 17681–17689.
- Merka, J., A. Szabo, J. A. Slavin, and M. Peredo (2005), Three-dimensional position and shape of the bow
  shock and their variation with upstream Mach numbers and interplanetary magnetic field orientation,
  *J. Geophys. Res.*, 110, A04202.
- Shue, J.-H., et al. (1998), Magnetopause location under extreme solar wind conditions,
  *J. Geophys. Res.*, 103(A8), 17691–17700.
- Tsyganenko, N. A. (1995), Modeling the Earth's magnetospheric magnetic field confined within a realistic
  magnetopause, *J. Geophys. Res.*, 100(A4), 5599–5612.
- Tsyganenko, N. A., and D. P. Stern (1996), Modeling the global magnetic field of the large-scale Birkeland
  current systems, *J. Geophys. Res.*, 101(A12), 27187–27198.

## License

MIT; see [LICENSE](LICENSE). Contributions: [CONTRIBUTING.md](CONTRIBUTING.md). Development notes:
[DEVELOPMENT.md](DEVELOPMENT.md).
