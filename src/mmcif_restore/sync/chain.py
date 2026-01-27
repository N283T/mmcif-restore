"""Chain category synchronization for mmCIF files."""

import gemmi

from mmcif_restore.structure_info import StructureInfo


def _keep_rows_by_id(
    block: gemmi.cif.Block,
    category_prefix: str,
    id_column: str,
    keep_ids: frozenset[str],
) -> None:
    """Keep only rows where the ID column value is in keep_ids.

    Args:
        block: CIF block to modify
        category_prefix: Category prefix (e.g., "_struct_asym.")
        id_column: Column name for ID (e.g., "id")
        keep_ids: Set of IDs to keep
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

    # Get column index for ID
    id_col_tag = f"{category_prefix}{id_column}"
    try:
        id_col_idx = list(target_loop.tags).index(id_col_tag)
    except ValueError:
        return

    # Get all current values
    all_values = list(target_loop.values)
    width = target_loop.width()
    num_rows = target_loop.length()

    # Build new columns keeping only rows with valid IDs
    new_columns: list[list[str]] = [[] for _ in range(width)]
    for row_idx in range(num_rows):
        row_start = row_idx * width
        row_id = all_values[row_start + id_col_idx]
        if row_id in keep_ids:
            for col_idx in range(width):
                new_columns[col_idx].append(all_values[row_start + col_idx])

    # Replace loop values
    target_loop.set_all_values(new_columns)


def sync_chain_categories(block: gemmi.cif.Block, info: StructureInfo) -> None:
    """Synchronize chain-related categories to match current structure.

    Keeps only rows for chains that exist in the structure.

    Categories synced:
    - _struct_asym

    Args:
        block: CIF block to modify
        info: Current structure information
    """
    # Sync _struct_asym (chain ID is in "id" column)
    _keep_rows_by_id(block, "_struct_asym.", "id", info.chain_ids)
