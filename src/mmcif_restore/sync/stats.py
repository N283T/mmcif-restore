"""Statistics category synchronization for mmCIF files."""

import gemmi


def sync_stats_categories(
    block: gemmi.cif.Block,
    structure: gemmi.Structure,
) -> None:
    """Synchronize statistics categories based on current structure.

    Updates counts in:
    - _refine_hist (atom counts)

    Args:
        block: CIF block to modify
        structure: Edited structure with current atom counts
    """
    # Count atoms in current structure
    atom_count = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                atom_count += len(residue)

    # Update _refine_hist if it exists
    # This is optional - not all CIF files have this category
    refine_hist_loop = None
    for item in block:
        if item.loop is not None:
            tags = item.loop.tags
            if tags and tags[0].startswith("_refine_hist."):
                refine_hist_loop = item.loop
                break

    if refine_hist_loop is None:
        return

    # Find the number_atoms_total column
    try:
        atom_col_idx = list(refine_hist_loop.tags).index(
            "_refine_hist.number_atoms_total"
        )
    except ValueError:
        # Column doesn't exist
        return

    # Get current values and update
    all_values = list(refine_hist_loop.values)
    width = refine_hist_loop.width()
    num_rows = refine_hist_loop.length()

    # Build new columns with updated atom count
    new_columns: list[list[str]] = [[] for _ in range(width)]
    for row_idx in range(num_rows):
        row_start = row_idx * width
        for col_idx in range(width):
            if col_idx == atom_col_idx:
                new_columns[col_idx].append(str(atom_count))
            else:
                new_columns[col_idx].append(all_values[row_start + col_idx])

    # Replace loop values
    refine_hist_loop.set_all_values(new_columns)
