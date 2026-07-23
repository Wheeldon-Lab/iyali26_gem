"""Compatibility wrapper for the safe EC-overload patch runner."""

from scripts.gem_annotate.patch_runner import main_for_legacy


if __name__ == "__main__":
    raise SystemExit(main_for_legacy("ec-overload-cleanup"))
