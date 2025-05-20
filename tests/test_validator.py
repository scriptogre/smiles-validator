from typing import Annotated

import pytest
from pydantic import BaseModel, ValidationError

from smiles_validator.validator import SmilesText, SmilesValidator


@pytest.fixture
def validator():
    """Fixture for a default SmilesValidator instance"""
    return SmilesValidator()


@pytest.fixture
def validator_keep_original():
    """Fixture for a SmilesValidator with keep_original=True"""
    return SmilesValidator(keep_original=True)


def test_canonical_smiles(validator):
    """Test canonical mode returns canonicalized SMILES"""
    assert validator("C1=CC=CC=C1") == "c1ccccc1"


def test_keep_original_smiles(validator_keep_original):
    """Test keep_original mode returns original SMILES"""
    assert validator_keep_original("C1=CC=CC=C1") == "C1=CC=CC=C1"


def test_invalid_smiles_raises_error(validator):
    """Test that invalid SMILES raises ValueError"""
    with pytest.raises(ValueError):
        validator("invalid_smiles")


def test_non_string_input_raises_error(validator):
    """Test that non-string input raises TypeError"""
    with pytest.raises(TypeError):
        validator(123)  # Passing an integer instead of string


def test_pydantic_integration():
    """Test integration with Pydantic model"""

    class MoleculeModel(BaseModel):
        canonical: SmilesText
        original: Annotated[str, SmilesValidator(keep_original=True)]

    # Valid case
    model = MoleculeModel(canonical="C1=CC=CC=C1", original="C1=CC=CC=C1")
    assert model.canonical == "c1ccccc1"
    assert model.original == "C1=CC=CC=C1"

    # Invalid case
    with pytest.raises(ValidationError):
        MoleculeModel(canonical="invalid_smiles", original="invalid_smiles")


@pytest.fixture(autouse=True)
def clear_cache():
    """Clear the LRU cache before each test"""
    from smiles_validator.validator import SmilesValidator

    SmilesValidator._process.cache_clear()


def test_caching_behavior(validator):
    """Test that validation results are cached using lru_cache"""
    test_smiles = "C1=CC=CC=C1"
    expected_canonical = "c1ccccc1"

    # First call - should process and cache
    result1 = validator(test_smiles)
    assert result1 == expected_canonical

    # Check cache stats after first call
    cache_info = validator._process.cache_info()
    assert cache_info.misses == 1
    assert cache_info.hits == 0

    # Second call with same input - should use cache
    result2 = validator(test_smiles)
    assert result2 == expected_canonical

    # Check cache stats after second call
    cache_info = validator._process.cache_info()
    assert cache_info.misses == 1
    assert cache_info.hits == 1

    # Clear cache and verify it's empty
    validator.clear_cache()
    cache_info = validator._process.cache_info()
    assert cache_info.misses == 0
    assert cache_info.hits == 0

    # Call with same SMILES again - should miss cache
    validator(test_smiles)
    cache_info = validator._process.cache_info()
    assert cache_info.misses == 1
    assert cache_info.hits == 0
