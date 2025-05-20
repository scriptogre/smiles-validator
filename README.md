# Smiles-Validator

[![Test Coverage](https://img.shields.io/badge/coverage-95%25-brightgreen.svg)](https://github.com/scriptogre/smiles-validator)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyPI version](https://img.shields.io/pypi/v/smiles-validator.svg)](https://pypi.org/project/smiles-validator/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/scriptogre/smiles-validator/actions/workflows/ci.yml/badge.svg)](https://github.com/scriptogre/smiles-validator/actions)
[![Code Style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A Pydantic validator for SMILES chemical notation that validates, sanitizes, and optionally returns either canonical or original input.

## Install

```bash
pip install smiles-validator
```

> Requires Python ≥ 3.10 · Pydantic v2 · RDKit 2024.9+

## Usage

```python
from pydantic import BaseModel
from typing import Annotated
from smiles_validator import SmilesText, SmilesValidator

class Model(BaseModel):
    # Default: Canonical SMILES
    canonical: SmilesText
    # Keep original input
    original: Annotated[str, SmilesValidator(keep_original=True)]

m = Model(
    canonical="C1=CC=CC=C1",
    original="C1=CC=CC=C1"
)

print(m.canonical)  # => "c1ccccc1"
print(m.original)   # => "C1=CC=CC=C1"
```

## Features

- 🚀 Fast SMILES validation using RDKit
- 🔄 Canonicalization with optional original input preservation
- 📦 Pydantic v2 integration
- 📊 Comprehensive test coverage
- 🧪 CI/CD with multiple Python versions

## API

### SmilesValidator Class

```python
SmilesValidator(keep_original: bool = False)
```

Parameters:
- `keep_original`: If `True`, returns the exact input string after validation. If `False` (default), returns the canonical SMILES.

### SmilesText Type

`SmilesText` is an alias for `Annotated[str, SmilesValidator()]` that automatically uses the default validator settings (canonical output).

## Development

```bash
# Install dependencies
uv sync --locked --all-extras --dev

# Run tests
just test

# Run linting
just check

# Run formatting
just format
```

## License

[MIT License](LICENSE) · [GitHub](https://github.com/scriptogre/smiles-validator)

## Acknowledgements

- Built with [RDKit](https://github.com/rdkit/rdkit) for chemical structure handling
- Uses [Pydantic v2](https://github.com/pydantic/pydantic) for validation integration
