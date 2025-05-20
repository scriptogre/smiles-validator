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

# Release a new version
release version:
    @echo "Bumping version to {{ version }}"
    sed -i '' 's/version = "[^"]*"/version = "{{ version }}"/' pyproject.toml
    git add pyproject.toml
    git commit -m "Bump version to {{ version }}"
    git tag v{{ version }}
    git push origin master
    git push --tags
