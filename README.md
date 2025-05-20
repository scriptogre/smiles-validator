# Smiles-Validator

[![Test Coverage](https://img.shields.io/badge/coverage-95%25-brightgreen.svg)](https://github.com/scriptogre/smiles-validator)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyPI version](https://img.shields.io/pypi/v/smiles-validator.svg)](https://pypi.org/project/smiles-validator/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/scriptogre/smiles-validator/actions/workflows/ci.yml/badge.svg)](https://github.com/scriptogre/smiles-validator/actions)
[![Code Style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Downloads](https://static.pepy.tech/personalized-badge/smiles-validator?period=total&units=international_system&left_color=grey&right_color=blue&left_text=downloads)](https://pepy.tech/project/smiles-validator)

Lightweight Pydantic v2 validator for RDKit SMILES strings that validates, sanitizes, and optionally returns canonical or original input. Built with RDKit for efficient chemical structure handling.

## Install

```bash
pip install smiles-validator
```

> Requires Python ≥ 3.10 · Pydantic v2 · RDKit 2024.9+

## Usage

### Basic Usage

```python
from pydantic import BaseModel
from typing import Annotated
from smiles_validator import SmilesText, SmilesValidator

class Model(BaseModel):
    # Default: Canonical SMILES
    canonical: SmilesText
    # Keep original input
    original: Annotated[str, SmilesValidator(keep_original=True)]

# Create a model instance
m = Model(
    canonical="C1=CC=CC=C1",
    original="C1=CC=CC=C1"
)

# Access validated/canonicalized values
print(m.canonical)  # => "c1ccccc1"
print(m.original)   # => "C1=CC=CC=C1"
```

### Error Handling

The validator raises appropriate exceptions for invalid input:

```python
from smiles_validator import SmilesValidator

try:
    # Invalid SMILES
    SmilesValidator()("invalid_smiles")
except ValueError as e:
    print(f"Error: {e}")  # => "Error: Invalid SMILES: 'invalid_smiles'"

try:
    # Incorrect type
    SmilesValidator()(123)
except TypeError as e:
    print(f"Error: {e}")  # => "Error: SMILES must be a string, got int"
```

### Advanced Usage

You can use the validator directly without Pydantic:

```python
from smiles_validator import SmilesValidator

# Create a validator instance
validator = SmilesValidator()

# Validate and canonicalize
result = validator("C1=CC=CC=C1")
print(result)  # => "c1ccccc1"

# Preserve original input
original_validator = SmilesValidator(keep_original=True)
result = original_validator("C1=CC=CC=C1")
print(result)  # => "C1=CC=CC=C1"
```

## Features

- 🚀 Fast SMILES validation using RDKit
- 🔄 Canonicalization with optional original input preservation
- 📦 Pydantic v2 integration
- 📊 Comprehensive test coverage
- 🧪 CI/CD with multiple Python versions
- 📚 Well-documented API
- 🔄 Automatic caching of validation results
- 🔒 Type-safe validation with Pydantic
- 📈 Efficient memory usage with LRU caching

## API

### SmilesValidator Class

```python
SmilesValidator(keep_original: bool = False)
```

Parameters:
- `keep_original`: If `True`, returns the exact input string after validation. If `False` (default), returns the canonical SMILES.

### SmilesText Type

`SmilesText` is an alias for `Annotated[str, SmilesValidator()]` that automatically uses the default validator settings (canonical output).

### Validation Process

1. Input validation:
   - Checks if input is a string
   - Validates SMILES syntax
   - Ensures chemical validity

2. Sanitization:
   - Uses RDKit's default sanitization
   - Handles aromaticity perception
   - Performs valence checks

3. Output:
   - Returns canonical SMILES (default)
   - Returns original input (if `keep_original=True`)

### Performance Considerations

- Uses LRU caching with a capacity of 4096 entries
- Cache is automatically cleared when the validator is garbage collected
- Caching improves performance for repeated validations

## Development

### Setup

```bash
# Install dependencies
uv sync --all-extras --dev

# Run tests
just test

# Run linting
just check

# Run formatting
just format
```

### Testing

The project includes comprehensive tests covering:
- Valid SMILES validation
- Invalid SMILES handling
- Error cases
- Caching behavior
- Pydantic integration

### Contributing

1. Fork the repository
2. Create a feature branch
3. Run tests: `just test`
4. Run linting: `just check`
5. Format code: `just format`
6. Submit a pull request

## License

[MIT License](LICENSE) · [GitHub](https://github.com/scriptogre/smiles-validator)

## Troubleshooting

### Common Issues

1. **RDKit Installation**
   - Ensure RDKit is properly installed
   - Check RDKit version compatibility

2. **Validation Errors**
   - Verify input is a valid SMILES string
   - Check for unsupported chemical structures

3. **Performance**
   - Clear cache if memory usage is high
   - Consider batch processing for large datasets

## Acknowledgements

- Built with [RDKit](https://github.com/rdkit/rdkit) for chemical structure handling
- Uses [Pydantic v2](https://github.com/pydantic/pydantic) for validation integration
- Inspired by the Pydantic ecosystem's focus on type safety and validation
