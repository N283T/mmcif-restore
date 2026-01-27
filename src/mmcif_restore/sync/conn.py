"""Connection category synchronization for mmCIF files."""

import logging

import gemmi

logger = logging.getLogger(__name__)

# Key for identifying an atom in _struct_conn
# (auth_asym_id, auth_seq_id, pdbx_PDB_ins_code, auth_comp_id, label_atom_id)
AtomKey = tuple[str, str, str, str, str]


def _build_atom_set(structure: gemmi.Structure) -> set[AtomKey]:
    """Build a set of atom keys from the structure.

    Args:
        structure: Gemmi Structure to analyze

    Returns:
        Set of (auth_asym_id, auth_seq_id, ins_code, comp_id, atom_id) tuples
    """
    atoms: set[AtomKey] = set()

    for model in structure:
        for chain in model:
            auth_asym_id = chain.name
            for residue in chain:
                auth_seq_id = str(residue.seqid.num)
                ins_code = residue.seqid.icode.strip() if residue.seqid.icode else ""
                comp_id = residue.name
                for atom in residue:
                    atoms.add((auth_asym_id, auth_seq_id, ins_code, comp_id, atom.name))

    return atoms


def sync_conn_categories(block: gemmi.cif.Block, structure: gemmi.Structure) -> None:
    """Synchronize connection-related categories to match current structure.

    Keeps only rows where both partners' atoms exist in the structure.
    Uses auth_* identifiers and insertion codes for accurate matching.

    Categories synced:
    - _struct_conn (disulfide bonds, hydrogen bonds, metal coordination, etc.)

    Args:
        block: CIF block to modify
        structure: Current structure for atom existence checking
    """
    # Build set of existing atoms
    existing_atoms = _build_atom_set(structure)

    if not existing_atoms:
        return

    # Find _struct_conn loop
    target_loop = None
    for item in block:
        if item.loop is not None:
            tags = item.loop.tags
            if tags and tags[0].startswith("_struct_conn."):
                target_loop = item.loop
                break

    if target_loop is None:
        logger.debug("No _struct_conn loop found")
        return

    # Get column indices for partner identifiers
    tags_list = list(target_loop.tags)

    # Required columns for both partners
    required_cols = [
        "ptnr1_auth_asym_id",
        "ptnr1_auth_seq_id",
        "ptnr1_auth_comp_id",
        "ptnr1_label_atom_id",
        "ptnr2_auth_asym_id",
        "ptnr2_auth_seq_id",
        "ptnr2_auth_comp_id",
        "ptnr2_label_atom_id",
    ]

    # Optional insertion code columns
    ins_code_cols = [
        "pdbx_ptnr1_PDB_ins_code",
        "pdbx_ptnr2_PDB_ins_code",
    ]

    try:
        col_indices = {
            col: tags_list.index(f"_struct_conn.{col}") for col in required_cols
        }
    except ValueError as e:
        logger.debug("Required column not found in _struct_conn: %s", e)
        return

    # Get optional insertion code indices (may not exist)
    for col in ins_code_cols:
        try:
            col_indices[col] = tags_list.index(f"_struct_conn.{col}")
        except ValueError:
            col_indices[col] = None

    # Get all current values
    all_values = list(target_loop.values)
    width = target_loop.width()
    num_rows = target_loop.length()

    # Build new columns keeping only connections where both atoms exist
    new_columns: list[list[str]] = [[] for _ in range(width)]

    for row_idx in range(num_rows):
        row_start = row_idx * width

        # Get partner 1 atom key
        ptnr1_ins_idx = col_indices.get("pdbx_ptnr1_PDB_ins_code")
        ptnr1_ins = ""
        if ptnr1_ins_idx is not None:
            ptnr1_ins = all_values[row_start + ptnr1_ins_idx]
            if ptnr1_ins in (".", "?"):
                ptnr1_ins = ""

        ptnr1_key: AtomKey = (
            all_values[row_start + col_indices["ptnr1_auth_asym_id"]],
            all_values[row_start + col_indices["ptnr1_auth_seq_id"]],
            ptnr1_ins,
            all_values[row_start + col_indices["ptnr1_auth_comp_id"]],
            all_values[row_start + col_indices["ptnr1_label_atom_id"]],
        )

        # Get partner 2 atom key
        ptnr2_ins_idx = col_indices.get("pdbx_ptnr2_PDB_ins_code")
        ptnr2_ins = ""
        if ptnr2_ins_idx is not None:
            ptnr2_ins = all_values[row_start + ptnr2_ins_idx]
            if ptnr2_ins in (".", "?"):
                ptnr2_ins = ""

        ptnr2_key: AtomKey = (
            all_values[row_start + col_indices["ptnr2_auth_asym_id"]],
            all_values[row_start + col_indices["ptnr2_auth_seq_id"]],
            ptnr2_ins,
            all_values[row_start + col_indices["ptnr2_auth_comp_id"]],
            all_values[row_start + col_indices["ptnr2_label_atom_id"]],
        )

        # Keep only if both partner atoms exist
        if ptnr1_key in existing_atoms and ptnr2_key in existing_atoms:
            for col_idx in range(width):
                new_columns[col_idx].append(all_values[row_start + col_idx])

    # Replace loop values
    target_loop.set_all_values(new_columns)
