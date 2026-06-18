"""
Python bindings for minibwa — fast BWA-based alignment from Python.

Classes:
    Index   — build or load a minibwa index.
    Opts    — alignment options (preset or manual).
    Hit     — one alignment result.
    Meth    — methylation-strand enum passed to ``map``.

Functions:
    map       — align one read; returns ``list[Hit]``.
    map_pair  — align a read pair; returns ``(list[Hit], list[Hit])``.
    map_many  — align many reads in one batch; returns ``list[list[Hit]]``.
"""

import enum

from ._minibwa import Hit  # noqa: A004
from ._minibwa import Index
from ._minibwa import Opts
from ._minibwa import map  # noqa: A004
from ._minibwa import map_many
from ._minibwa import map_pair


class Meth(str, enum.Enum):
    """
    Methylation conversion mode for a read.

    Because ``Meth`` subclasses ``str``, values pass through directly to the
    Rust ``&str`` matcher — ``minibwa.map(..., meth=Meth.C2T)`` is identical
    to ``minibwa.map(..., meth="c2t")``.
    """

    NONE = "none"
    """No methylation conversion (default)."""
    C2T = "c2t"
    """C→T conversion (read 1 in bisulfite/EM-seq)."""
    G2A = "g2a"
    """G→A conversion (read 2 in bisulfite/EM-seq)."""


__all__ = ["Index", "Opts", "Hit", "Meth", "map", "map_pair", "map_many"]
