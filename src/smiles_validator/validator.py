from functools import lru_cache
from typing import Annotated, Any

from pydantic import GetCoreSchemaHandler
from pydantic_core import core_schema
from rdkit import Chem


class SmilesValidator:
    """
    A Pydantic-compatible callable that validates SMILES strings using RDKit.

    Usage examples:
      smiles: Annotated[str, SmilesValidator()]
      smiles: SmilesText

    Parameters:
        keep_original: whether to return the user's original SMILES (default False)
    """

    def __init__(self, *, keep_original: bool = False) -> None:
        self.keep_original = keep_original

    # Class-level cache to prevent memory leaks
    _process_cache = lru_cache(maxsize=4096)

    @_process_cache
    def _process(self, value: str, canonical: bool = True) -> str:
        """Process a SMILES string and return its canonical form if requested.

        Args:
            value: The SMILES string to process
            canonical: Whether to return the canonical form (default: True)

        Returns:
            The processed SMILES string

        Raises:
            ValueError: If the SMILES is invalid or sanitization fails
        """
        # 1) Parse without sanitization
        mol = Chem.MolFromSmiles(value, sanitize=False)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {value!r}")

        # 2) Sanitize if requested
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            raise ValueError(f"SMILES sanitization failed for {value!r}: {e}") from e

        # 3) Return canonical form
        if canonical:
            return Chem.MolToSmiles(mol, canonical=True)
        return Chem.MolToSmiles(mol)

    def clear_cache(self) -> None:
        """Clear the cache of processed SMILES strings."""
        self._process.cache_clear()

    def __call__(self, value: Any) -> str:
        if not isinstance(value, str):
            raise TypeError(f"SMILES must be a string, got {type(value).__name__}")

        # If user wants to keep their exact input, skip canonicalization
        if self.keep_original:
            # Still validate by processing, but ignore its output
            self._process(value, canonical=False)
            return value
        # Otherwise, return cached canonical SMILES
        return self._process(value, canonical=True)

    def __get_pydantic_core_schema__(
        self, source: type[str], handler: GetCoreSchemaHandler
    ) -> core_schema.CoreSchema:
        # We need to modify how the validation is done for Pydantic
        # Use a closure to capture self.keep_original
        keep_original = self.keep_original

        def validate_smiles(val, info):
            # Create or use our instance with the correct keep_original value
            validator = SmilesValidator(keep_original=keep_original)
            return validator(val)

        return core_schema.with_info_plain_validator_function(validate_smiles)


# Type alias for SMILES strings
type SmilesText = Annotated[str, SmilesValidator(keep_original=False)]
