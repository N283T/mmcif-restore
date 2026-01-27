"""Main restore function for mmcif-restore."""

from pathlib import Path

import gemmi

from mmcif_restore.structure_info import StructureInfo
from mmcif_restore.sync.chain import sync_chain_categories
from mmcif_restore.sync.conn import sync_conn_categories
from mmcif_restore.sync.entity import sync_entity_categories
from mmcif_restore.sync.scheme import sync_scheme_categories


def restore_categories(
    edited_path: str | Path,
    reference_path: str | Path,
    categories: list[str],
) -> gemmi.cif.Document:
    """Restore mmCIF categories from reference CIF to an edited structure.

    This function takes a CIF file that has been edited (e.g., by ChimeraX,
    PyMOL, or Gemmi) and restores specified categories from the original
    reference CIF, synchronized to match the current structure.

    Args:
        edited_path: Path to the edited CIF file (can have minimal categories)
        reference_path: Path to the original reference CIF with all categories
        categories: List of category prefixes to restore
                   (e.g., ["_entity.", "_struct_asym.", "_struct_conn."])

    Returns:
        A new CIF Document with the restored and synchronized categories

    Example:
        # After editing a structure in ChimeraX and saving as minimal CIF:
        doc = restore_categories(
            "edited.cif",
            "original.cif.gz",
            categories=["_entity.", "_struct_asym.", "_struct_conn."]
        )
        doc.write_file("restored.cif")
    """
    edited_path = Path(edited_path)
    reference_path = Path(reference_path)

    # Load reference CIF (needed for entity mapping)
    ref_doc = gemmi.cif.read(str(reference_path))
    ref_block = ref_doc[0]

    # Read structure from edited CIF and extract info using reference
    structure = gemmi.read_structure(str(edited_path))
    info = StructureInfo.from_structure_with_reference(structure, ref_block)

    # Load edited CIF as base document
    doc = gemmi.cif.read(str(edited_path))
    block = doc[0]

    # Normalize category prefixes (ensure they end with ".")
    categories = [c if c.endswith(".") else c + "." for c in categories]

    # Restore each category from reference with sync
    for category_prefix in categories:
        _restore_and_sync_category(ref_block, block, category_prefix, info)

    return doc


def _restore_and_sync_category(
    source_block: gemmi.cif.Block,
    target_block: gemmi.cif.Block,
    category_prefix: str,
    info: StructureInfo,
) -> None:
    """Restore a category from source to target, synced to structure info.

    Args:
        source_block: Reference CIF block to copy from
        target_block: Edited CIF block to copy to
        category_prefix: Category prefix (e.g., "_struct_conn.")
        info: Current structure information for filtering
    """
    # Find the loop in source
    source_loop = None
    for item in source_block:
        if item.loop is not None:
            tags = list(item.loop.tags)
            if tags and tags[0].startswith(category_prefix):
                source_loop = item.loop
                break

    if source_loop is None:
        return

    # Copy to target
    tags = list(source_loop.tags)
    col_names = [tag.replace(category_prefix, "") for tag in tags]

    new_loop = target_block.init_mmcif_loop(category_prefix, col_names)

    values = list(source_loop.values)
    width = len(tags)
    num_rows = len(values) // width
    columns = [
        [values[row * width + col] for row in range(num_rows)] for col in range(width)
    ]
    new_loop.set_all_values(columns)

    # Apply sync based on category type
    if category_prefix == "_struct_conn.":
        sync_conn_categories(target_block, info)
    elif category_prefix == "_entity.":
        sync_entity_categories(target_block, info)
    elif category_prefix == "_struct_asym.":
        sync_chain_categories(target_block, info)
    elif category_prefix in (
        "_pdbx_nonpoly_scheme.",
        "_pdbx_poly_seq_scheme.",
        "_pdbx_branch_scheme.",
    ):
        sync_scheme_categories(target_block, info)
