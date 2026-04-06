"""Utilities for horizontal and vertical extrapolation on the ISMIP grid."""

from importlib import resources as _resources

from i7aof.extrap.horizontal import (
    extrapolate_horizontal,
    horizontal_extrapolation,
)
from i7aof.extrap.vertical import extrapolate_vertical, vertical_extrapolation

__all__ = [
    'extrapolate_horizontal',
    'extrapolate_vertical',
    'get_template_path',
    'horizontal_extrapolation',
    'load_template_text',
    'vertical_extrapolation',
]


def get_template_path() -> str:
    """Return the filesystem path to the bundled Jinja2 namelist template."""
    template = _resources.files(__package__) / 'namelist_template.nml.j2'
    with _resources.as_file(template) as p:
        return str(p)


def load_template_text() -> str:
    """Return the contents of the bundled Jinja2 namelist template."""
    template = _resources.files(__package__) / 'namelist_template.nml.j2'
    return template.read_text(encoding='utf-8')
