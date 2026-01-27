"""Synchronization modules for mmCIF categories."""

__all__ = [
    "sync_entity_categories",
    "sync_chain_categories",
    "sync_scheme_categories",
    "sync_stats_categories",
    "sync_conn_categories",
]


def __getattr__(name: str):
    """Lazy import synchronization functions."""
    if name == "sync_entity_categories":
        from mmcif_restore.sync.entity import sync_entity_categories

        return sync_entity_categories
    if name == "sync_chain_categories":
        from mmcif_restore.sync.chain import sync_chain_categories

        return sync_chain_categories
    if name == "sync_scheme_categories":
        from mmcif_restore.sync.scheme import sync_scheme_categories

        return sync_scheme_categories
    if name == "sync_stats_categories":
        from mmcif_restore.sync.stats import sync_stats_categories

        return sync_stats_categories
    if name == "sync_conn_categories":
        from mmcif_restore.sync.conn import sync_conn_categories

        return sync_conn_categories
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
