"""Compatibility entry point for the canonical ``gem_annotate`` build."""

import importlib


def main(*args, **kwargs):
    """Forward dynamically so existing callers and tests remain compatible."""

    pipeline = importlib.import_module("scripts.gem_annotate.main")
    return pipeline.main(*args, **kwargs)


if __name__ == "__main__":
    main()
