# Run all checks and tests
default: check format test


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

# Release a new version (usage: just release 0.1.17)
release version:
    @echo "Bumping version to {{ version }}"

    # Update pyproject.toml
    sed -i '' 's/version = "[^"]*"/version = "{{ version }}"/' pyproject.toml

    # Stage, commit and create new tag
    git add pyproject.toml
    git commit -m "Bump version to {{ version }}"
    git tag v{{ version }}

    # Push master branch & new tag
    git push origin master
    git push --tags
