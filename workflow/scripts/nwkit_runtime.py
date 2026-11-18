"""Shared NWKIT version contract for rooting, selection, and calibration."""

NWKIT_VERSION = "0.46.11"


def nwkit_backend():
    import nwkit

    if nwkit.__version__ != NWKIT_VERSION:
        raise ValueError(f"use the pinned nwkit {NWKIT_VERSION} in workflow/envs/dating.yaml")
    return nwkit
