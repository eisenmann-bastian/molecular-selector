"""Functional-group selection for molecular structures."""

from .batch import BatchSelector
from .definitions import NamedSelector, SelectorRegistry, carboxylate_sites
from .matching import MolecularSelector
from .output import CsvResultWriter
from .patterns import *
from .patterns import __all__ as _pattern_exports
from .interactions import (
    BOND_TYPE_ROLES, INTERACTION_ROLES, InteractionRole, InteractionType, get_groups_for_role,
)
from .sites import (
    INTERACTION_SITES, MAPPED_FUNCTIONAL_GROUP_SMARTS, InteractionSite, get_interaction_site,
)

__all__ = [
    *_pattern_exports,
    "BatchSelector", "CsvResultWriter", "MolecularSelector", "NamedSelector",
    "SelectorRegistry", "carboxylate_sites",
    "BOND_TYPE_ROLES", "INTERACTION_ROLES", "InteractionRole", "InteractionType",
    "get_groups_for_role", "INTERACTION_SITES", "MAPPED_FUNCTIONAL_GROUP_SMARTS",
    "InteractionSite", "get_interaction_site",
]
