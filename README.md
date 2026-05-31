# Wind Bow Shock Database: Data Repository

This repository contains data from 69 passes of the Wind spacecraft through Earth's bow shock.

A **pass** is defined as the portion of the spacecraft trajectory that traverses the theoretical bow shock region. Each pass may contain one or more bow shock crossings. Passes are classified as either **inbound** or **outbound**.

A **crossing** is the moment in time when the spacecraft crosses the bow shock.

The database is provided in the **database_files/** directory and includes:

- An Excel file containing all identified crossings
- An HDF5 (`.h5`) file containing the same information
- Four visualization files showing all crossings in GSE coordinates

Additional files are organized into two directories:

- **figs_data/** contains plots of all data products for each pass, with crossing times marked.
- **figs_html/** contains interactive visualizations of each pass.


# Wind Bow Shock Database: Magnetosphere Visualization with Wind Spacecraft Data

A Python-based tool for visualizing magnetosphere boundaries (magnetopause and bow shock) using Wind spacecraft magnetic field and solar wind measurements. Combines real satellite data with magnetosphere models (T96, Shue 1998) to create interactive 3D visualizations of magnetospheric structure during magnetosphere-solar wind crossings.

## Features

- **Wind Spacecraft Data Integration**: Fetches MFI (magnetic field) and SWE (solar wind) measurements from CDAWEB
- **Dynamic Magnetosphere Modeling**: Uses T96 field model with Shue 1998 magnetopause and Mach number-dependent bow shock
- **Interactive 3D Visualization**: Plotly-based interactive plots with spacecraft trajectory, field lines, and boundaries
- **Batch Processing**: Process multiple magnetosphere crossings automatically
- **Database-Driven**: Built-in bow shock crossing database with timing windows

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
├── MagCarto_Master.ipynb        # Main Jupyter notebook with all functionality
├── bs_crossings_loader.py       # Crossing database access module
├── requirements.txt             # Python package dependencies
├── README.md                    # This file
├── .gitignore                   # Git ignore patterns
├── scripts/                     # Standalone utility scripts
│   ├── magnetosphere_wind_plot.py
│   ├── magnetosphere_wind_batch.py
│   ├── magnetosphere_hs_plot.py
│   ├── magnetosphere_hs_batch.py
│   └── ...
├── figs/                        # Output figures (generated, not tracked)
├── figs_wind/                   # Wind-derived figures (generated, not tracked)
└── data/                        # Local data cache (not tracked)
```

## Usage

### Processing a Single Crossing

The notebook cell `PASS_DIRECTION_3D_PLOTLY` processes a single magnetosphere crossing:

```python
# Select pass and direction interactively
# Crossing time is extracted automatically
# Wind MFI and SWE data are fetched from cached CDF files
# 3D visualization is generated showing:
# - Wind spacecraft trajectory (blue line)
# - Magnetopause surface (yellow mesh)
# - Bow shock surface (cyan surface)
# - Magnetic field lines (colored by magnitude)
```

### Batch Processing Multiple Crossings

The `BATCH_WIND_CORRECTED` cell processes all 18 Wind-derived passes:

```python
# Automatically:
# - Loads crossing database
# - Maps CDF files to crossing times
# - Extracts Pdyn, Dst, Bz, By, Ma parameters
# - Generates HTML visualizations
# - Saves to figs_wind/ directory
# - Tracks successes and failures
```

## Data Sources

- **Wind Spacecraft Data**: CDAWEB (wi_h1_mfi_*.cdf, wi_h1_swe_*.cdf)
- **Disturbance Storm Time**: NOAA Space Weather Prediction Center (dst*.txt)
- **Magnetosphere Model Parameters**: T96 model with IGRF magnetic field
- **Crossing Database**: bs_crossings_loader.load_crossings()

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

### Wind Data Location

Edit the data path in the notebook cells:
```python
data_dir = r"D:\Data\Wind\bow_shock\data_used"  # Update to your Wind CDF files location
```

### Time Tolerances

For Wind data extraction with time gaps:
```python
mfi_tolerance = 21600   # ±6 hours for MFI data
swe_tolerance = 21600   # ±6 hours for SWE data
```

## Output

Visualizations are saved as interactive HTML files:
- `figs_wind/pass_NN_DIR.html` - Individual crossing visualizations
- Each plot includes spacecraft trajectory, magnetosphere boundaries, and field lines

## Troubleshooting

### Missing CDF Files

If Wind data is unavailable for a crossing:
1. Check CDAWEB for wi_h1_mfi_*.cdf and wi_h1_swe_*.cdf files
2. Download files and place in the configured data directory
3. Ensure filename format: `mfi_YYYY_YYYYMMDD_vXX.cdf`

### Import Errors

Ensure all dependencies are installed:
```bash
pip install -r requirements.txt
python -c "import spacepy; import geopack; import sscws"
```

### Jupyter Notebook Issues

Clear notebook checkpoints and restart:
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

- Tsyganenko, N. A. (1996), Modeling the Earth's magnetospheric magnetic field confined within a realistic magnetopause, J. Geophys. Res., 101, 27187–27198.
- Shue, J. H., et al. (1998), Magnetopause location under extreme solar wind conditions, J. Geophys. Res., 103, 17691–17700.

## Contact

For questions or issues, please open an issue on GitHub.
