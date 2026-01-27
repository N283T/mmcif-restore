"""Scheme category synchronization for mmCIF files."""

import gemmi

from mmcif_restore.structure_info import StructureInfo


def _keep_rows_by_asym_id(
    block: gemmi.cif.Block,
    category_prefix: str,
    asym_id_column: str,
    keep_ids: frozenset[str],
) -> None:
    """Keep only rows where asym_id is in keep_ids.

    Args:
        block: CIF block to modify
        category_prefix: Category prefix
        asym_id_column: Column name for asym ID
        keep_ids: Set of asym IDs (chain/subchain) to keep
    """
    if not keep_ids:
        return

    # Find the loop object containing this category
    target_loop = None
    for item in block:
        if item.loop is not None:
            tags = item.loop.tags
            if tags and tags[0].startswith(category_prefix):
                target_loop = item.loop
                break

    if target_loop is None:
        return

    # Get column index for asym_id
    asym_col_tag = f"{category_prefix}{asym_id_column}"
    try:
        asym_col_idx = list(target_loop.tags).index(asym_col_tag)
    except ValueError:
        return

    # Get all current values
    all_values = list(target_loop.values)
    width = target_loop.width()
    num_rows = target_loop.length()

    # Build new columns keeping only rows with valid asym_ids
    new_columns: list[list[str]] = [[] for _ in range(width)]
    for row_idx in range(num_rows):
        row_start = row_idx * width
        asym_id = all_values[row_start + asym_col_idx]
        if asym_id in keep_ids:
            for col_idx in range(width):
                new_columns[col_idx].append(all_values[row_start + col_idx])

    # Replace loop values
    target_loop.set_all_values(new_columns)


def sync_scheme_categories(block: gemmi.cif.Block, info: StructureInfo) -> None:
    """Synchronize scheme-related categories to match current structure.

    Keeps only rows for chains that exist in the structure.

    Categories synced:
    - _pdbx_poly_seq_scheme
    - _pdbx_nonpoly_scheme
    - _pdbx_branch_scheme

    Args:
        block: CIF block to modify
        info: Current structure information
    """
    keep_chain_ids = info.chain_ids

    # Sync _pdbx_poly_seq_scheme
    _keep_rows_by_asym_id(block, "_pdbx_poly_seq_scheme.", "asym_id", keep_chain_ids)

    # Sync _pdbx_nonpoly_scheme
    _keep_rows_by_asym_id(block, "_pdbx_nonpoly_scheme.", "asym_id", keep_chain_ids)

    # Sync _pdbx_branch_scheme (for branched entities like carbohydrates)
    _keep_rows_by_asym_id(block, "_pdbx_branch_scheme.", "asym_id", keep_chain_ids)
