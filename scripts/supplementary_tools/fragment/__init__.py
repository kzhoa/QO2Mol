"""Molecular fragmentation with selectable algorithms."""

from .api import process, BACKENDS, DEFAULT_BACKEND
from .modern import FragmentConfig


FAMILIES = (
    "elem_frag",
    "seco_frag",
    "csf_frag",
    "rsf_frag",
    "chain_frag",
    "scaffold_frag",
    "tf_frag",
    "sssr_frag",
    "sketch_frag",
)

__all__ = [
    "process",
    "BACKENDS",
    "DEFAULT_BACKEND",
    "FAMILIES",
    "FragmentConfig",
]
