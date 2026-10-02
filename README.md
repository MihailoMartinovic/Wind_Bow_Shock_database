# Wind Bow Shock Database: Magnetosphere Visualization with Wind Spacecraft Data

A Python-based tool for visualizing magnetosphere boundaries (magnetopause and bow shock) using Wind spacecraft passes through Earth's bow shock. The repository contains **635 bow shock crossings** across **70 passes** spanning August 1995 to April 2004. Interactive 3D visualizations combine real satellite data with magnetosphere models (T96, Shue 1998) to show magnetospheric structure during Wind's crossings.

## Database Overview

- **635 individual crossings** identified across 70 passes of the Wind spacecraft
- **Temporal Range**: 1995-08 to 2004-04 (~9 years of observations)
- **Pass Classification**: Each pass is either inbound or outbound through the bow shock region
- **Data Sources**: OMNI2 hourly solar wind/IMF data, NOAA Dst index, T96 magnetosphere model with IGRF field

**Key files**:
- `pass_parameters_v2.csv` — Provenance record for all 140 pass legs (includes 3 skipped legs with incomplete OMNI data)
- `database_files/` — Crossing database in Excel and HDF5 formats
- `figs/` — Data product plots for each pass with crossing markers  
- `figs_deprecated/` — Previous-generation renderings (kept for comparison)

## Features

- **OMNI2 Upstream Parameters**: All 137 rendered passes use OMNI2 hourly solar wind/IMF data (single, consistent data source)
- **Shue 1998 Magnetopause Model**: Dynamic standoff distance computed from solar wind pressure and Bz
- **Gas-Dynamic Bow Shock**: Alfvén Mach number-dependent standoff with 0.9-eccentricity conic profile
- **Interactive 3D Visualization**: Plotly-based HTML renderings with spacecraft trajectory, magnetic field lines, and boundary surfaces
- **T96 Field Lines**: 96 field lines traced via geopack T96 model in both directions (rlim = 60 R⊕)
- **Batch Processing**: Resumable notebook (`render_passes_v2.ipynb`) processes all passes, skipping completed renderings
- **Complete Provenance**: `pass_parameters_v2.csv` documents all upstream parameters, model inputs, and computation metadata

## Installation

### Quick Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/Wind_Bow_Shock_database.git
cd Wind_Bow_Shock_database

# Create virtual environment
python3 -m venv venv-magnetosphere
source venv-magnetosphere/bin/activate  # On Windows: venv-magnetosphere\Scripts\activate

# Install dependencies
pip install -U pip wheel setuptools
pip install -r requirements.txt

# Update OMNI2 data (one-time setup)
python -c "import spacepy.toolbox as tb; tb.update(omni2=True)"
```

### Optional Dependencies (for video generation)

```bash
pip install pyvista imageio imageio-ffmpeg
```

## Quick Start

### Running the Jupyter Notebook

```bash
jupyter notebook MagCarto_Master.ipynb
```

The notebook includes:
1. **Pass Selection Widget**: Interactive dropdown to select magnetosphere crossings
2. **Batch Processing**: Process multiple passes with automatic visualization generation
3. **3D Plotting**: View magnetosphere boundaries with spacecraft trajectory and field lines
4. **Parameter Extraction**: Automatically extract solar wind parameters and magnetic field values

### Example: Single Pass Visualization

```python
import sys
sys.path.append('.')
from bs_crossings_loader import load_crossings

# Load bow shock crossing database
db = load_crossings()

# Get a specific crossing
pass_id = 42
crossing = db.get_crossings_for_pass(pass_id)[0]

# Run MagCarto_Master.ipynb cells for visualization
```

## Project Structure

```
Wind_Bow_Shock_database/
├── MagCarto_Master.ipynb          # Main Jupyter notebook with interactive widget interface
├── render_passes_v2.ipynb         # Generator for 137 interactive HTML renderings
├── build_omni2_cache.ipynb        # Builds OMNI2 cache from D:\Data\OMNI\omni2_cache.pkl
├── bs_crossings_loader.py         # Crossing database access module
├── pass_parameters_v2.csv         # Provenance record: 140 rows (137 rendered + 3 skipped)
├── requirements.txt               # Python package dependencies
├── README.md                      # This file
├── BS_CROSSINGS_DATABASE.md       # Data schema and methodology documentation
├── .gitignore                     # Git ignore patterns
├── scripts/                       # Standalone utility scripts
│   ├── magnetosphere_wind_plot.py
│   ├── magnetosphere_wind_batch.py
│   └── ...
├── database_files/                # Crossing database (Excel, HDF5, visualization data)
├── figs/                          # 137 data product plots with crossing markers
├── figs_deprecated/               # Previous-generation renderings (v1.0.0)
└── data/                          # Local data cache (not tracked)
```

## Usage

### Interactive Visualization (MagCarto_Master.ipynb)

Run the main notebook to select and visualize individual bow shock crossings:

```bash
jupyter notebook MagCarto_Master.ipynb
```

Features:
- **Pass Selection Widget**: Choose pass number and direction (inbound/outbound)
- **Parameter Display**: Shows upstream OMNI2 parameters, model inputs, and boundary locations
- **3D Interactive Plot**: Rotate, zoom, and pan to explore magnetosphere boundaries and spacecraft trajectory
- **Export**: Save rendered HTML for sharing

### Batch Processing (render_passes_v2.ipynb)

Generate all 137 interactive renderings (takes ~4 hours):

```bash
jupyter notebook render_passes_v2.ipynb
```

The notebook:
- Reads OMNI2 data from cached pickle file (build with `build_omni2_cache.ipynb` if missing)
- Processes 140 pass legs, skipping 3 with incomplete OMNI data
- Generates one HTML file per leg in `figs/`
- Resumes from last completed leg if interrupted
- Logs provenance to `pass_parameters_v2.csv`

## Magnetosphere Models

### Magnetopause (Shue et al. 1998)
The magnetopause standoff distance and shape are computed from solar wind dynamic pressure and Bz:

```
r₀ = (10.22 + 1.29 × tanh(0.184 × (Bz + 8.14))) × Pdyn^(-1/6.6)
α = (0.58 - 0.007 × Bz) × (1 + 0.024 × ln(Pdyn))
r(θ) = r₀ × (2 / (1 + cos(θ)))^α
```

### Bow Shock (Gas-Dynamic Standoff)
The bow shock is a conic section (eccentricity e = 0.9) with subsolar distance determined by the magnetopause nose plus gas-dynamic standoff:

```
Δ = r₀ × ((γ - 1) × M_A² + 2) / ((γ + 1) × (M_A² - 1)),  γ = 5/3
r_bs0 = r₀ + Δ
r(θ) = r_bs0 × (1 + e) / (1 + e × cos(θ))
```

**Important**: The Mach number used is the **Alfvén Mach number (M_A)**, not the magnetosonic Mach number.

### Field Lines
96 field lines are traced using the T96 magnetosphere model with IGRF magnetic field coefficients, traced in both directions from each seed with a radial limit of 60 R⊕.

**All surfaces are axisymmetric** about the GSM x-axis (no aberration, dawn-dusk asymmetry, or dipole tilt effects).

## Data Sources

- **OMNI2 Hourly Data**: Solar wind velocity, dynamic pressure (Pdyn), IMF components (Bz, By)
- **NOAA Dst Index**: Disturbance Storm Time index for magnetospheric activity
- **Magnetosphere Models**: T96 (Tsyganenko 1996) field model with IGRF coefficients
- **Crossing Database**: Wind bow shock crossing times from `bs_crossings_loader.py`

## Key Dependencies

| Package | Purpose |
|---------|---------|
| `spacepy` | Magnetosphere modeling, OMNI2 data access |
| `geopack` | T96 magnetosphere field model |
| `cdflib` | CDF file format reading |
| `sscws` | SSCWeb satellite position queries |
| `plotly` | Interactive 3D visualization |
| `pandas`, `numpy` | Data manipulation and arrays |
| `ipywidgets` | Interactive Jupyter widgets |

## Configuration

### OMNI2 Cache Location

The batch renderer uses OMNI2 hourly data cached locally:
```python
omni_cache = r"D:\Data\OMNI\omni2_cache.pkl"  # Must be built once with build_omni2_cache.ipynb
```

### Parameter Search Window

For each pass leg, the generator finds the nearest OMNI2 hour within ±4 hours of the first crossing at which all five parameters (Pdyn, Bz, By, Dst, M_A) are present. This window can be adjusted in `render_passes_v2.ipynb`:

```python
MAX_OFFSET_H = 4.0  # Hours before/after first crossing
```

### Alfvén Mach Number Guard

Legs with M_A < 1.2 are considered unreliable (gas-dynamic relation diverges as M_A → 1):
```python
MIN_MACH = 1.2
```
Currently, 137 of 140 legs meet this criterion.

## Output

Visualizations are saved as interactive HTML files:
- `figs_wind/pass_NN_DIR.html` - Individual crossing visualizations
- Each plot includes spacecraft trajectory, magnetosphere boundaries, and field lines

## Troubleshooting

### Missing OMNI2 Cache

If `render_passes_v2.ipynb` fails to find OMNI2 data:
1. Run `build_omni2_cache.ipynb` to fetch and cache OMNI2 hourly data
2. Ensure `D:\Data\OMNI\omni2_cache.pkl` is created and readable
3. Check spacepy is installed: `python -c "import spacepy.omni"`

### Incomplete Legs

Three legs lack complete OMNI2 data within the ±4 h search window:
- Pass 4 inbound (1995-11-28 08:40 UTC): 3 crossings skipped
- Pass 55 inbound (2001-02-20 02:45 UTC): 1 crossing skipped
- Pass 56 inbound (2001-09-08 12:05 UTC): 7 crossings skipped

These are logged as `status = "skipped: no complete OMNI hour within 4.0 h"` in `pass_parameters_v2.csv`.

### Import Errors

Ensure all dependencies are installed:
```bash
pip install -r requirements.txt
python -c "import spacepy; import geopack; import plotly"
```

### Notebook Checkpoint Issues

Clear old checkpoints and restart kernel:
```bash
rm -rf .ipynb_checkpoints
jupyter kernel restart
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Wind spacecraft data provided by CDAWeb (NASA)
- Magnetosphere models: T96 (Tsyganenko 1996), Shue 1998 magnetopause
- Bow shock database development and magnetosphere crossing analysis

## References

- Farris, M. H., and C. T. Russell (1994), Bow shock and magnetosheath dynamic pressure variations relative to the solar wind dynamic pressure, J. Geophys. Res., 99, 17681–17689.
  - Gas-dynamic standoff relation (evaluated with Alfvén Mach number)
- Shue, J. H., et al. (1998), Magnetopause location under extreme solar wind conditions, J. Geophys. Res., 103, 17691–17700.
  - Shue et al. magnetopause model (Eqs. 10–11)
- Tsyganenko, N. A. (1996), Modeling the Earth's magnetospheric magnetic field confined within a realistic magnetopause, J. Geophys. Res., 101, 27187–27198.
  - T96 magnetosphere field model

## Known Issues and Limitations

- **Pass 47 inbound** (M_A = 1.40) renders with artificially large standoff distance (Δ = 16.85 R⊕) due to near-divergence of the gas-dynamic relation. Model surfaces are annotated as unreliable for this leg.
- **Axisymmetric model**: No aberration from Earth's orbital motion, no dawn-dusk asymmetry, no dipole-tilt dependence.
- **Skipped legs**: Three pass legs lack complete OMNI2 data and are not rendered (see Troubleshooting).

## Documentation

For detailed information on the crossing database schema, data processing pipeline, and parameter sources, see **BS_CROSSINGS_DATABASE.md**.

## Contact

For questions or issues, please open an issue on GitHub.
