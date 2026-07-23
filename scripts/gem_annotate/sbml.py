"""Deterministic, failure-safe SBML serialization helpers.

COBRApy 0.30 stores ``Group.members`` in a :class:`set` and iterates that set
directly in its SBML writer.  Python deliberately randomises set iteration
between processes, so two builds of an unchanged model can otherwise produce
different byte-level SHA-256 values solely because pathway members move around
in the XML.  The evidence workflow uses the raw model SHA for provenance, so
that serialization detail must be stable.
"""

from __future__ import annotations

import os
import stat
import tempfile
from pathlib import Path
from typing import Any

from cobra.io import write_sbml_model


def _group_member_sort_key(member: Any) -> tuple[str, str]:
    """Return a stable key for every COBRA object allowed in an SBML group."""

    return (type(member).__name__, str(member.id))


def write_deterministic_sbml_model(
    model: Any,
    filename: str | Path,
    **writer_kwargs: Any,
) -> None:
    """Write *model* atomically with deterministic SBML Group member order.

    The wrapper temporarily exposes each group's members to COBRApy as a
    sorted tuple, invokes the standard writer, and restores the exact original
    member container in a ``finally`` block.  The completed temporary file is
    then atomically moved into place, so a failed serialization cannot truncate
    the current model.
    """

    target = Path(filename)
    target.parent.mkdir(parents=True, exist_ok=True)
    target_mode = (
        stat.S_IMODE(target.stat().st_mode) if target.exists() else 0o644
    )
    # libSBML selects gzip/bzip2/zip output from the filename suffix, so the
    # temporary path must retain all of the target's suffixes.
    temporary_suffix = "".join(target.suffixes) or ".tmp"
    descriptor, temporary_name = tempfile.mkstemp(
        dir=target.parent,
        prefix=f".{target.name}.",
        suffix=temporary_suffix,
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)

    original_members: list[tuple[Any, Any]] = []
    try:
        try:
            for group in model.groups:
                original_members.append((group, group._members))
                group._members = tuple(
                    sorted(group.members, key=_group_member_sort_key)
                )
            write_sbml_model(model, str(temporary_path), **writer_kwargs)
        finally:
            for group, members in original_members:
                group._members = members

        os.chmod(temporary_path, target_mode)
        os.replace(temporary_path, target)
    finally:
        temporary_path.unlink(missing_ok=True)
