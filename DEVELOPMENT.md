# Development Guide

This guide covers setup and development practices for Wind Bow Shock Database.

## Development Environment Setup

### 1. Clone and Initial Setup

```bash
git clone https://github.com/MihailoMartinovic/Wind_Bow_Shock_database.git
cd Wind_Bow_Shock_database
```

### 2. Virtual Environment

```bash
# Create virtual environment
python3 -m venv venv-magnetosphere

# Activate it
source venv-magnetosphere/bin/activate  # macOS/Linux
# or
venv-magnetosphere\Scripts\activate  # Windows
```

### 3. Install Dependencies

```bash
# Upgrade pip, wheel, setuptools
pip install -U pip wheel setuptools

# Install core dependencies
pip install -r requirements.txt

# Install development tools
pip install pytest pytest-cov black flake8 jupyter jupyterlab

# Download spacepy's OMNI2 database (one-time). MagCarto_Master.ipynb reads it through spacepy.omni,
# and scripts/build_omni2_cache.ipynb builds the renderings' OMNI2 cache from it
# (~/.spacepy/data/omni2data.h5).
python -c "import spacepy.toolbox as tb; tb.update(omni2=True)"
```

### 4. Verify Installation

```bash
# Test imports
python -c "import pandas, openpyxl, tables, h5py, spacepy, geopack, sscws, cdflib, plotly; print('✓ All imports successful')"

# Start Jupyter
jupyter lab
# Open MagCarto_Master.ipynb, or a notebook in scripts/
```

## Paths and inputs

The notebooks in `scripts/` and `MagCarto_Master.ipynb` read machine-specific locations from `local_paths.json`
in the folder they run in. Copy `scripts/local_paths.example.json` to `local_paths.json` and fill it in;
`.gitignore` keeps it out of the repository. Two inputs of the rendering are not in the repository:

- the OMNI2 cache `omni2_cache.pkl`, written by `scripts/build_omni2_cache.ipynb`;
- the *Wind* ephemeris `Wind_Ephemerids.h5`, key `/GSM2` (GSM positions at 10-minute cadence).

The scripts `scripts/magnetosphere_hs_*.py` read the HS-RT run folder from the environment variable `HS_RT_DIR`.

## Regenerating the data products

After any change to `database_files/bs_crossings_V02.xlsx`, run in this order:

1. `scripts/build_hdf5_from_crossings.ipynb` rewrites `database_files/Wind_bow_shock_database.h5` from the
   list, keeping the previous file's key, columns and types.
2. `scripts/render_passes_v2.ipynb` re-renders the affected legs: list them in `ONLY_LEGS` in the first code
   cell. Only those legs' rows of `pass_parameters_v2.csv` change.
3. `scripts/remake_figure4.ipynb`, then `scripts/remake_figure5.ipynb`.

These notebooks never delete what they replace. They archive it in an `00_old_versions` folder beside the
repository (path set in their first code cell), under rules described in that folder's README.

A leg is a run of rows with the same `pass` and `SC direction`. Do not identify legs by row order or by
which cells are filled.

## Code Style

### Formatting

Use `black` for automatic code formatting:

```bash
# Format all Python files
black .

# Check formatting without changing
black --check .
```

### Linting

Check code quality with `flake8`:

```bash
flake8 .
```

### Docstrings

Use NumPy-style docstrings:

```python
def example_function(param1, param2):
    """
    Brief description of function.

    Longer description explaining what the function does,
    including any important details.

    Parameters
    ----------
    param1 : type
        Description of param1
    param2 : type
        Description of param2

    Returns
    -------
    result : type
        Description of return value

    Examples
    --------
    >>> example_function(1, 2)
    3
    """
    return param1 + param2
```

## Working with Jupyter Notebooks

### Notebook Cleanup

Before committing notebooks, strip outputs:

```bash
pip install nbstripout
nbstripout MagCarto_Master.ipynb
```

Or configure git hooks:

```bash
nbstripout --install --attributes .gitattributes
```

### Notebook Best Practices

1. **Cell Organization**: Group related cells logically
2. **Markdown Documentation**: Use markdown cells to explain sections
3. **Comments**: Add inline comments for complex logic
4. **Variable Names**: Use descriptive names (avoid `a`, `b`, `c`)
5. **Output Clearing**: Clear outputs before committing
6. **Magic Commands**: Document any IPython magic commands used

## Testing

### Running Tests (future)

```bash
pytest tests/ -v
pytest tests/ --cov=. --cov-report=html
```

### Manual Testing

For now, testing is done through Jupyter notebook cells:

1. Run individual cells to test functions
2. Test with sample data from the crossing database
3. Verify output visualizations visually

## Debugging

### Jupyter Debugging

```python
# Insert in notebook cell to drop to debugger
%pdb on
# or
import pdb; pdb.set_trace()
```

### Print Debugging

```python
# Simple approach
print(f"Debug: var_name = {var_name}")

# Structured approach
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.debug(f"var_name = {var_name}")
```

### Logging

Add logging to scripts:

```python
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

logger.info("Processing started")
logger.warning("Potential issue detected")
logger.error("Error occurred", exc_info=True)
```

## Git Workflow

### Feature Development

```bash
# Create feature branch
git checkout -b feature/my-feature

# Make changes and commit
git add .
git commit -m "Add my feature"

# Push to fork
git push origin feature/my-feature

# Create pull request on GitHub
```

### Commit Message Format

```
Short summary (50 chars or less)

Detailed explanation wrapped at 72 characters. Explain what
the change does and why it's necessary.

Fixes #123
Related to #456
```

## Documentation

- **README.md**: User-facing documentation, including the table of releases
- **BS_CROSSINGS_DATABASE.md**: The crossing list: columns, conventions and changes
- **CONTRIBUTING.md**: Contribution guidelines
- **DEVELOPMENT.md**: This file - development guide
- **Notebook markdown**: Each notebook in `scripts/` states its inputs, outputs and expected results

### Building Documentation

For future Sphinx documentation:

```bash
pip install sphinx sphinx-rtd-theme
sphinx-quickstart docs
cd docs
make html
```

## Common Tasks

### Adding a New Feature

1. Create feature branch: `git checkout -b feature/new-feature`
2. Implement in notebook or create new script
3. Test thoroughly
4. Update README.md if user-facing
5. Update CONTRIBUTING.md if affects workflow
6. Commit and create pull request

### Updating Dependencies

1. Test with new version locally
2. Update `requirements.txt` with new versions
3. Document breaking changes in PR

### Creating a Release

Tag the release commit itself, after it exists, and check where the tag points before pushing it. Zenodo
archives the tagged commit.

```bash
git add -A
git commit -m "Release vX.Y.Z: summary"
git tag -a vX.Y.Z -m "vX.Y.Z: summary"      # only after the commit has succeeded
git rev-list -n 1 vX.Y.Z                    # must print the release commit
git log -1 --format=%H                      # ... which is this one
git push origin main
git push origin vX.Y.Z
git ls-remote origin                        # main and the tag as GitHub sees them
```

Then add the release to the table in README.md and create the release on GitHub from the tag.

## Performance Optimization

### Profiling

```python
import cProfile
import pstats

profiler = cProfile.Profile()
profiler.enable()

# Code to profile
...

profiler.disable()
stats = pstats.Stats(profiler)
stats.sort_stats('cumulative')
stats.print_stats(10)  # Top 10 functions
```

### Benchmarking

```python
import time

start = time.time()
# Code to benchmark
end = time.time()
print(f"Time: {end - start:.4f} seconds")
```

## Troubleshooting

### Import Errors

```bash
# Verify all dependencies installed
pip list | grep -E "pandas|openpyxl|tables|spacepy|geopack|sscws|cdflib"

# Reinstall if needed
pip install --force-reinstall spacepy
```

### Jupyter Kernel Issues

```bash
# Restart kernel in Jupyter
# Use "Restart Kernel" button

# Or from command line
jupyter kernelspec list
python -m ipykernel install --user --name venv-magnetosphere
```

### CDF File Issues

```python
import cdflib

# Inspect CDF file (recent cdflib versions have no close() method)
file = cdflib.CDF('filename.cdf')
print(file.cdf_info())          # List variables
print(file.varget('VAR_NAME'))  # Read variable
```

## Resources

- [Wind Bow Shock Database GitHub](https://github.com/MihailoMartinovic/Wind_Bow_Shock_database)
- [spacepy Documentation](https://spacepy.github.io/)
- [Jupyter Notebook Guide](https://jupyter-notebook.readthedocs.io/)
- [git Documentation](https://git-scm.com/doc)
- [GitHub Guides](https://guides.github.com/)
