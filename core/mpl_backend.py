from __future__ import annotations

import os
import tempfile


def configure_headless_matplotlib() -> None:
    """Force a non-GUI backend so rendering is safe in worker threads."""

    backend = os.environ.get("MPLBACKEND", "").strip().lower()
    if not backend:
        os.environ["MPLBACKEND"] = "Agg"

    if not os.environ.get("MPLCONFIGDIR"):
        os.environ["MPLCONFIGDIR"] = os.path.join(tempfile.gettempdir(), "pairwise_mplconfig")

    import matplotlib

    current = str(matplotlib.get_backend()).lower()
    if "agg" not in current:
        matplotlib.use("Agg", force=True)
