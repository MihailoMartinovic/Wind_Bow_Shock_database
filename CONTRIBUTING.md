# Contributing to Wind Bow Shock Database

Thank you for your interest in contributing! This document provides guidelines for contributing to the Wind Bow Shock Database project.

## Code of Conduct

Be respectful and inclusive. We welcome contributions from people of all backgrounds and experience levels.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/yourusername/Wind_Bow_Shock_database.git
   cd Wind_Bow_Shock_database
   ```
3. **Create a virtual environment**:
   ```bash
   python3 -m venv venv-magnetosphere
   source venv-magnetosphere/bin/activate
   ```
4. **Install development dependencies**:
   ```bash
   pip install -r requirements.txt
   pip install pytest pytest-cov black flake8
   ```

## Making Changes

### Branch Naming

Use descriptive branch names:
- `feature/description` for new features
- `fix/description` for bug fixes
- `docs/description` for documentation
- `refactor/description` for code refactoring

### Code Style

- Follow PEP 8 style guidelines
- Use meaningful variable and function names
- Add docstrings to functions and classes
- Format code with `black`:
  ```bash
  black .
  ```
- Check with `flake8`:
  ```bash
  flake8 .
  ```

### Jupyter Notebook Guidelines

- Keep notebooks well-commented
- Clear variable names and function definitions
- Group related cells logically
- Include markdown cells explaining functionality
- Clean output before committing (or use .gitignore to exclude)

### Commit Messages

Write clear commit messages:
```
Short summary (50 characters or less)

Detailed explanation if needed. Wrap at 72 characters.
Reference issues with "Fixes #123" or "Related to #123"
```

## Testing

Before submitting a pull request:

1. Test your changes thoroughly
2. Ensure existing functionality still works
3. Add tests for new features if applicable

## Submitting Changes

1. **Push to your fork**:
   ```bash
   git push origin feature/description
   ```

2. **Open a Pull Request** on GitHub with:
   - Clear title describing the change
   - Description of what was changed and why
   - Reference to any related issues
   - Screenshots for visualization changes

3. **Address feedback** from code review

## Reporting Issues

When reporting a bug, include:

- **Description**: Clear summary of the issue
- **Steps to Reproduce**: How to trigger the bug
- **Expected Behavior**: What should happen
- **Actual Behavior**: What actually happens
- **Environment**: Python version, OS, relevant packages
- **Error Message**: Full error traceback if applicable

## Documentation

- Update README.md for user-facing changes
- Add docstrings to code functions
- Update this file if contributing guidelines change

## Questions?

Open an issue or discussion on GitHub. We're happy to help!

## License

By contributing, you agree that your contributions will be licensed under the same license as the project (MIT License).
