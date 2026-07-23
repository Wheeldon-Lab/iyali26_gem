"""
gem_annotate — annotation pipeline for the iYali26 GEM of Yarrowia lipolytica.
"""


def main(*args, **kwargs):
    """Load the pipeline lazily so CLI path configuration takes effect first."""
    from .main import main as run_main

    return run_main(*args, **kwargs)

__all__ = ["main"]
