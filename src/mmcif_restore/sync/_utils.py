"""Shared utility functions for sync modules."""

import logging

import gemmi

logger = logging.getLogger(__name__)


def keep_rows_by_column(
    block: gemmi.cif.Block,
    category_prefix: str,
    filter_column: str,
    keep_values: frozenset[str],
) -> None:
    """Keep only rows where the filter column value is in keep_values.

    Args:
        block: CIF block to modify
        category_prefix: Category prefix (e.g., "_entity.")
        filter_column: Column name for filtering (e.g., "id" or "entity_id")
        keep_values: Set of values to keep
    """
    if not keep_values:
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
        logger.debug("No loop found for category %s", category_prefix)
        return

    # Get column index for filter
    filter_col_tag = f"{category_prefix}{filter_column}"
    try:
        filter_col_idx = list(target_loop.tags).index(filter_col_tag)
    except ValueError:
        logger.debug(
            "Column %s not found in category %s", filter_column, category_prefix
        )
        return

    # Get all current values
    all_values = list(target_loop.values)
    width = target_loop.width()
    num_rows = target_loop.length()

    # Build new columns keeping only rows with valid values
    new_columns: list[list[str]] = [[] for _ in range(width)]
    for row_idx in range(num_rows):
        row_start = row_idx * width
        row_value = all_values[row_start + filter_col_idx]
        if row_value in keep_values:
            for col_idx in range(width):
                new_columns[col_idx].append(all_values[row_start + col_idx])

    # Replace loop values
    target_loop.set_all_values(new_columns)
