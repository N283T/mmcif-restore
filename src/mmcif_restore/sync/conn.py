"""Connection category synchronization for mmCIF files."""

import gemmi

from mmcif_restore.structure_info import StructureInfo


def sync_conn_categories(block: gemmi.cif.Block, info: StructureInfo) -> None:
    """Synchronize connection-related categories to match current structure.

    Keeps only rows where both partners' chains exist in the structure.

    Categories synced:
    - _struct_conn (disulfide bonds, hydrogen bonds, metal coordination, etc.)

    Args:
        block: CIF block to modify
        info: Current structure information
    """
    keep_chain_ids = info.chain_ids
    if not keep_chain_ids:
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
        return

    # Get column indices for partner asym_ids
    tags_list = list(target_loop.tags)
    try:
        ptnr1_asym_idx = tags_list.index("_struct_conn.ptnr1_label_asym_id")
        ptnr2_asym_idx = tags_list.index("_struct_conn.ptnr2_label_asym_id")
    except ValueError:
        # Required columns not found
        return

    # Get all current values
    all_values = list(target_loop.values)
    width = target_loop.width()
    num_rows = target_loop.length()

    # Build new columns keeping only connections where both partners exist
    new_columns: list[list[str]] = [[] for _ in range(width)]
    for row_idx in range(num_rows):
        row_start = row_idx * width
        ptnr1_asym = all_values[row_start + ptnr1_asym_idx]
        ptnr2_asym = all_values[row_start + ptnr2_asym_idx]

        # Keep only if both partners' chains still exist
        if ptnr1_asym in keep_chain_ids and ptnr2_asym in keep_chain_ids:
            for col_idx in range(width):
                new_columns[col_idx].append(all_values[row_start + col_idx])

    # Replace loop values
    target_loop.set_all_values(new_columns)
