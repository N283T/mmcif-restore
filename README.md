# mmcif-restore

Restore mmCIF categories from a reference CIF to an edited structure.

## Problem

When you edit a PDB structure using tools like ChimeraX, PyMOL, or Gemmi (e.g., removing waters, deleting ligands), the exported mmCIF often loses important metadata categories like `_entity`, `_struct_asym`, `_struct_conn`, etc. Even when these categories are preserved, they may become out of sync with the actual structure.

## Solution

`mmcif-restore` takes an edited CIF file and a reference CIF (the original), then restores specified categories from the reference while **synchronizing** them to match the edited structure. For example, if you removed water molecules, the restored `_entity` category will not include the water entity.

## Installation

```bash
pip install mmcif-restore
```

Or with uv:

```bash
uv pip install mmcif-restore
```

## Usage

### Command Line

```bash
# Restore _entity and _struct_conn categories
mmcif-restore edited.cif original.cif.gz -o output.cif -c "_entity.,_struct_conn."

# Restore multiple categories
mmcif-restore edited.cif original.cif.gz -o output.cif \
    -c "_entity.,_struct_asym.,_struct_conn.,_pdbx_nonpoly_scheme."
```

### Python API

```python
from mmcif_restore import restore_categories

# Restore categories from reference
doc = restore_categories(
    "edited.cif",
    "original.cif.gz",
    categories=["_entity.", "_struct_conn."]
)

# Write output
doc.write_file("output.cif")
```

## Supported Categories

| Category | Sync Behavior |
|----------|---------------|
| `_entity` | Removes entities not present in edited structure |
| `_entity_poly` | Synced with `_entity` |
| `_entity_poly_seq` | Synced with `_entity` |
| `_pdbx_entity_nonpoly` | Synced with `_entity` |
| `_struct_asym` | Removes chains not present in edited structure |
| `_struct_conn` | Removes connections involving removed atoms/residues |
| `_pdbx_struct_mod_residue` | Removes MODRES records for removed residues |
| `_pdbx_poly_seq_scheme` | Synced with chain IDs |
| `_pdbx_nonpoly_scheme` | Synced with chain IDs |
| `_pdbx_branch_scheme` | Synced with chain IDs |

## Requirements

- Python >= 3.12
- [Gemmi](https://gemmi.readthedocs.io/) >= 0.7.0

## Development

```bash
# Clone the repository
git clone https://github.com/N283T/mmcif-restore.git
cd mmcif-restore

# Install with dev dependencies
uv sync --group dev

# Run tests
pytest

# Run tests with coverage
pytest --cov=mmcif_restore --cov-report=term-missing
```

## License

MIT
