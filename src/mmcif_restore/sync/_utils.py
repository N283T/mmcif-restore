"""Shared utility functions for sync modules."""

import logging
from typing import Any

import gemmi

logger = logging.getLogger(__name__)


def keep_rows_by_column(
    block: gemmi.cif.Block,
    category_prefix: str,
    filter_column: str,
    keep_values: frozenset[str],
) -> None:
    """Keep only rows where the filter column value is in keep_values.

    Uses get_mmcif_category/set_mmcif_category to handle both loop and pair formats.

    Args:
        block: CIF block to modify
        category_prefix: Category prefix (e.g., "_entity.")
        filter_column: Column name for filtering (e.g., "id" or "entity_id")
        keep_values: Set of values to keep
    """
    if not keep_values:
        return

    # Get category data as dictionary (works for both loop and pair formats)
    try:
        data: dict[str, list[Any]] = block.get_mmcif_category(category_prefix, raw=True)
    except RuntimeError:
        logger.debug("Category %s not found", category_prefix)
        return

    if not data:
        logger.debug("No data found for category %s", category_prefix)
        return

    # Check if filter column exists
    if filter_column not in data:
        logger.debug(
            "Column %s not found in category %s", filter_column, category_prefix
        )
        return

    filter_values = data[filter_column]
    num_rows = len(filter_values)

    # Find row indices to keep
    keep_indices = [i for i, val in enumerate(filter_values) if val in keep_values]

    # If all rows are kept, no changes needed
    if len(keep_indices) == num_rows:
        return

    # Build filtered data
    filtered_data: dict[str, list[Any]] = {}
    for col_name, col_values in data.items():
        filtered_data[col_name] = [col_values[i] for i in keep_indices]

    # Write back filtered data (handles both loop and pair formats)
    block.set_mmcif_category(category_prefix, filtered_data, raw=True)
