"""Python core routines for GenL scattering simulations."""

from .convolution import gauss_conv, pad_edge
from .crystal import ProjectedStructure, project_structure, write_genl_projection
from .form_factors import (
    ELEMENT_SYMBOLS,
    FormFactorCoefficients,
    form_factors,
    read_form_factor_coefficients,
)
from .kinematic import Control, Instrument, KinematicResult, Layer, calc_kinematic
from .dynamic import (
    DynamicResult,
    DynamicWorkspace,
    calc_dynamic_density,
    validate_density_sampling,
)
from .poscar import PoscarStructure, read_poscar
from .stack import STACK_FORMAT, STACK_VERSION, StackDefinition, StackModel

__all__ = [
    "Control",
    "DynamicResult",
    "DynamicWorkspace",
    "ELEMENT_SYMBOLS",
    "FormFactorCoefficients",
    "Instrument",
    "KinematicResult",
    "Layer",
    "PoscarStructure",
    "ProjectedStructure",
    "STACK_FORMAT",
    "STACK_VERSION",
    "StackDefinition",
    "StackModel",
    "calc_kinematic",
    "calc_dynamic_density",
    "form_factors",
    "gauss_conv",
    "pad_edge",
    "project_structure",
    "read_form_factor_coefficients",
    "read_poscar",
    "validate_density_sampling",
    "write_genl_projection",
]
