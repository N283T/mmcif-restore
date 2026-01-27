"""mmcif-restore: Restore mmCIF categories from reference CIF to edited structures."""

__version__ = "0.1.0"

from mmcif_restore.restore import restore_categories

__all__ = ["restore_categories"]
