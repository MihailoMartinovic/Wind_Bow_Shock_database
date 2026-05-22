# Development Guide

This guide covers setup and development practices for Wind Bow Shock Database.

## Development Environment Setup

### 1. Clone and Initial Setup

```bash
git clone https://github.com/yourusername/Wind_Bow_Shock_database.git
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

# Update OMNI2 data (one-time)
python -c "import spacepy.toolbox as tb; tb.update(omni2=True)"
```

### 4. Verify Installation

```bash
# Test imports
python -c "import spacepy, geopack, sscws, cdflib, plotly; print('✓ All imports successful')"

# Start Jupyter
jupyter lab
# Open MagCarto_Master.ipynb
```

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

## Environment Variables

For local development, create `.env.local`:

```
WIND_DATA_DIR=/path/to/wind/data
CACHE_DIR=/path/to/cache
DEBUG=True
```

Load with:

```python
from dotenv import load_dotenv
import os

load_dotenv('.env.local')
wind_data_dir = os.getenv('WIND_DATA_DIR')
```

## Documentation

- **README.md**: User-facing documentation
- **CONTRIBUTING.md**: Contribution guidelines
- **DEVELOPMENT.md**: This file - development guide
- **Notebook docstrings**: In-code documentation in the Jupyter notebook

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
3. Test in CI/CD
4. Document breaking changes in PR

### Creating a Release

```bash
# Update version in setup.py or pyproject.toml
# Update CHANGELOG
# Create git tag
git tag v1.0.0
git push origin v1.0.0
# Create release on GitHub
```

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
pip list | grep -E "spacepy|geopack|sscws|cdflib"

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

# Inspect CDF file
file = cdflib.CDF('filename.cdf')
print(file.variables)  # List variables
print(file['VAR_NAME'][:])  # Read variable
file.close()
```

## Resources

- [Wind Bow Shock Database GitHub](https://github.com/yourusername/Wind_Bow_Shock_database)
- [spacepy Documentation](https://spacepy.github.io/)
- [Jupyter Notebook Guide](https://jupyter-notebook.readthedocs.io/)
- [git Documentation](https://git-scm.com/doc)
- [GitHub Guides](https://guides.github.com/)
