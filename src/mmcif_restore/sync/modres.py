"""Modified residue category synchronization for mmCIF files."""

import logging

import gemmi

logger = logging.getLogger(__name__)

# CIF special values
CIF_INAPPLICABLE = "."
CIF_UNKNOWN = "?"

# Key for identifying a residue in _pdbx_struct_mod_residue
# (auth_asym_id, auth_seq_id, pdbx_PDB_ins_code, comp_id)
ResidueKey = tuple[str, str, str, str]
"""Residue identifier tuple: (auth_asym_id, auth_seq_id, ins_code, comp_id)."""


def _build_residue_set(structure: gemmi.Structure) -> set[ResidueKey]:
    """Build a set of residue keys from the structure.

    Args:
        structure: Gemmi Structure to analyze

    Returns:
        Set of (auth_asym_id, auth_seq_id, ins_code, comp_id) tuples
    """
    residues: set[ResidueKey] = set()

    for model in structure:
        for chain in model:
            auth_asym_id = chain.name
            for residue in chain:
                auth_seq_id = str(residue.seqid.num)
                ins_code = residue.seqid.icode.strip() if residue.seqid.icode else ""
                comp_id = residue.name
                residues.add((auth_asym_id, auth_seq_id, ins_code, comp_id))

    return residues


def sync_modres_categories(block: gemmi.cif.Block, structure: gemmi.Structure) -> None:
    """Synchronize modified residue category to match current structure.

    Keeps only rows where the referenced residue exists in the structure.
    Uses auth_* identifiers and insertion codes for accurate matching.

    Categories synced:
    - _pdbx_struct_mod_residue (MODRES records)

    Args:
        block: CIF block to modify
        structure: Current structure for residue existence checking
    """
    # Build set of existing residues
    existing_residues = _build_residue_set(structure)

    if not existing_residues:
        return

    # Find _pdbx_struct_mod_residue loop
    target_loop = None
    for item in block:
        if item.loop is not None:
            tags = item.loop.tags
            if tags and tags[0].startswith("_pdbx_struct_mod_residue."):
                target_loop = item.loop
                break

    if target_loop is None:
        logger.debug("No _pdbx_struct_mod_residue loop found")
        return

    # Get column indices for residue identifiers
    tags_list = list(target_loop.tags)

    # Required columns
    required_cols = [
        "auth_asym_id",
        "auth_seq_id",
    ]

    # Optional columns
    optional_cols = [
        "PDB_ins_code",
        "auth_comp_id",
        "label_comp_id",
    ]

    try:
        col_indices: dict[str, int | None] = {
            col: tags_list.index(f"_pdbx_struct_mod_residue.{col}")
            for col in required_cols
        }
    except ValueError as e:
        logger.debug("Required column not found in _pdbx_struct_mod_residue: %s", e)
        return

    # Get optional column indices (may not exist)
    for col in optional_cols:
        try:
            col_indices[col] = tags_list.index(f"_pdbx_struct_mod_residue.{col}")
        except ValueError:
            col_indices[col] = None

    # Get all current values
    all_values = list(target_loop.values)
    width = target_loop.width()
    num_rows = target_loop.length()

    # Build new columns keeping only rows with existing residues
    new_columns: list[list[str]] = [[] for _ in range(width)]

    for row_idx in range(num_rows):
        row_start = row_idx * width

        # Get auth_asym_id
        auth_asym_id = all_values[row_start + col_indices["auth_asym_id"]]  # type: ignore[operator]

        # Get auth_seq_id
        auth_seq_id = all_values[row_start + col_indices["auth_seq_id"]]  # type: ignore[operator]

        # Get insertion code (optional)
        ins_code = ""
        ins_idx = col_indices.get("PDB_ins_code")
        if ins_idx is not None:
            ins_code = all_values[row_start + ins_idx]
            if ins_code in (CIF_INAPPLICABLE, CIF_UNKNOWN):
                ins_code = ""

        # Get comp_id (try auth_comp_id first, then label_comp_id)
        comp_id = ""
        auth_comp_idx = col_indices.get("auth_comp_id")
        label_comp_idx = col_indices.get("label_comp_id")
        if auth_comp_idx is not None:
            comp_id = all_values[row_start + auth_comp_idx]
        need_fallback = not comp_id or comp_id in (CIF_INAPPLICABLE, CIF_UNKNOWN)
        if need_fallback and label_comp_idx is not None:
            comp_id = all_values[row_start + label_comp_idx]
        if comp_id in (CIF_INAPPLICABLE, CIF_UNKNOWN):
            comp_id = ""

        # Build residue key
        residue_key: ResidueKey = (auth_asym_id, auth_seq_id, ins_code, comp_id)

        # Keep only if residue exists
        if residue_key in existing_residues:
            for col_idx in range(width):
                new_columns[col_idx].append(all_values[row_start + col_idx])

    # Replace loop values
    target_loop.set_all_values(new_columns)
