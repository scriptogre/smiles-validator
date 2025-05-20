# Run all checks and tests
default: check test


# Run Ruff linting
check:
    uv tool run ruff check .

# Run Ruff formatting
format:
    uv tool run ruff format .

# Run tests with coverage
test:
    uv run pytest --cov=smiles_validator tests/


# Install dependencies
install:
    uv sync --locked --all-extras --dev
