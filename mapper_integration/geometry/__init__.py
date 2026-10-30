"""GDML integration helpers for TMS magnetic-field maps."""

from .gdml_tools import (
    AlignmentReport,
    MapperHeader,
    TMSGeometryInfo,
    compute_field_alignment,
    globalize_mapper,
    inject_arb_bfield,
    inspect_tms_geometry,
    validate_injected_gdml,
)

__all__ = [
    "AlignmentReport",
    "MapperHeader",
    "TMSGeometryInfo",
    "compute_field_alignment",
    "globalize_mapper",
    "inject_arb_bfield",
    "inspect_tms_geometry",
    "validate_injected_gdml",
]
