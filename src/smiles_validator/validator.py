from functools import lru_cache
from typing import Annotated, Any

from pydantic import GetCoreSchemaHandler
from pydantic_core import core_schema
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover


class SmilesValidator:
    """
    A Pydantic-compatible callable that validates (and optionally canonicalises)
    SMILES strings with RDKit.  The molecule is *desalted by default*.

    Usage examples
    --------------
        smiles: Annotated[str, SmilesValidator()]     # defaults: return_canonical=True
        smiles: SmilesText                             # same as above
        smiles: Annotated[str, SmilesValidator(return_canonical=False)]

    Parameters
    ----------
    return_canonical : bool, default True
        Return RDKit’s canonical SMILES (if False, returns RDKit’s “as written”
        SMILES after validation and desalting).
    """

    def __init__(
        self, *, return_canonical: bool = True, return_desalted: bool = True
    ) -> None:
        self.return_canonical = return_canonical
        self.return_desalted = return_desalted

    # class-level cache
    _process_cache = lru_cache(maxsize=4096)

    @_process_cache
    def _process(self, value: str) -> str:
        # 1) Parse (no sanitisation yet)
        mol = Chem.MolFromSmiles(value, sanitize=False)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {value!r}")

        # 2) Sanitise
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            raise ValueError(f"SMILES sanitisation failed for {value!r}: {e}") from e

        # 3) Desalt
        if self.return_desalted:
            mol = SaltRemover().StripMol(mol)

        # 4) Return in the desired form
        return Chem.MolToSmiles(mol, canonical=self.return_canonical)

    def clear_cache(self) -> None:  # pragma: no cover
        self._process.cache_clear()

    def __call__(self, value: Any) -> str:
        if not isinstance(value, str):
            raise TypeError(f"SMILES must be a string, got {type(value).__name__}")
        return self._process(value)

    # ---- Pydantic integration -------------------------------------------------

    def __get_pydantic_core_schema__(
        self, source: type[str], handler: GetCoreSchemaHandler
    ) -> core_schema.CoreSchema:
        def validate_smiles(val, info):
            return SmilesValidator(return_canonical=self.return_canonical)(val)

        return core_schema.with_info_plain_validator_function(validate_smiles)


# Convenient type alias
type SmilesText = Annotated[str, SmilesValidator()]
