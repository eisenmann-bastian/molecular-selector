"""Named functional-group and interaction-site selector definitions."""

import json
import re
from dataclasses import dataclass

from rdkit import Chem

from .patterns import (
    CARBONYL_SITE_SMARTS,
    CARBOXYLATE_OXYGEN_PAIR_SMARTS,
    METAL_CENTER_SITE_SMARTS,
    NEGATIVE_CHARGE_SITE_SMARTS,
    POSITIVE_CHARGE_SITE_SMARTS,
    FunctionalGroup,
)
from .interactions import InteractionRole, InteractionType
from .sites import INTERACTION_SITES, MAPPED_FUNCTIONAL_GROUP_SMARTS, get_interaction_site


@dataclass(frozen=True)
class SelectorQuery:
    """A SMARTS query and optional atom-map labels for its highlighted sites."""

    smarts: str
    atom_map_numbers: tuple[int, ...] | None = None


@dataclass(frozen=True)
class NamedSelector:
    """A named query, its aliases, and optionally its component queries."""

    smarts: str | None
    aliases: tuple[str, ...] = ()
    queries: tuple[SelectorQuery, ...] = ()


def _named_queries(queries, aliases=()):
    queries = tuple(queries)
    if not queries:
        raise ValueError("a selector requires at least one query")
    description = queries[0].smarts if len(queries) == 1 else json.dumps(
        [query.smarts for query in queries], separators=(",", ":")
    )
    return NamedSelector(description, aliases, queries)


def _group_query(group, labels=None):
    return SelectorQuery(MAPPED_FUNCTIONAL_GROUP_SMARTS[group], labels)


def _role_queries(interaction, role):
    sites = INTERACTION_SITES[interaction][role]
    return tuple(
        SelectorQuery(sites[group].smarts, sites[group].atom_map_numbers)
        for group in sorted(sites, key=lambda group: group.value)
    )


def _group_site(group, interaction, role):
    site = get_interaction_site(interaction, role, group)
    return SelectorQuery(site.smarts, site.atom_map_numbers)


_HBOND = InteractionType.HYDROGEN_BOND
_DONOR = InteractionRole.DONOR
_ACCEPTOR = InteractionRole.ACCEPTOR


def _default_group_query(group):
    labels = set()
    for role in (_DONOR, _ACCEPTOR):
        site = INTERACTION_SITES[_HBOND][role].get(group)
        if site is not None:
            labels.update(site.atom_map_numbers)
    return _group_query(group, tuple(sorted(labels)) if labels else None)


_GROUP_ALIASES = {
    FunctionalGroup.THIOETHER: ("sulfide",),
    FunctionalGroup.CARBOXYLATE: ("carboxylate anion",),
}
SELECTORS = {
    group.value.replace("_", " "): _named_queries(
        [_default_group_query(group)], _GROUP_ALIASES.get(group, ())
    )
    for group in FunctionalGroup
}

SELECTORS["alcohol"] = _named_queries([
    _group_site(group, _HBOND, _DONOR) for group in (
        FunctionalGroup.PRIMARY_ALCOHOL,
        FunctionalGroup.SECONDARY_ALCOHOL,
        FunctionalGroup.TERTIARY_ALCOHOL,
    )
])
SELECTORS["amine"] = _named_queries([
    _group_site(group, _HBOND, _ACCEPTOR) for group in (
        FunctionalGroup.PRIMARY_AMINE,
        FunctionalGroup.SECONDARY_AMINE,
        FunctionalGroup.TERTIARY_AMINE,
    )
])
SELECTORS["amide"] = _named_queries([
    _group_query(group, (2, 3)) for group in (
        FunctionalGroup.PRIMARY_AMIDE,
        FunctionalGroup.SECONDARY_AMIDE,
        FunctionalGroup.TERTIARY_AMIDE,
    )
])

_INTERACTION_NAMES = {
    (InteractionType.HYDROGEN_BOND, _DONOR): ("hydrogen bond donor", "h bond donor", "hbond donor"),
    (InteractionType.HYDROGEN_BOND, _ACCEPTOR): ("hydrogen bond acceptor", "h bond acceptor", "hbond acceptor"),
    (InteractionType.HALOGEN_BOND, InteractionRole.SIGMA_HOLE_DONOR): ("halogen bond donor", "halogen donor"),
    (InteractionType.HALOGEN_BOND, _ACCEPTOR): ("halogen bond acceptor", "halogen acceptor"),
    (InteractionType.CHALCOGEN_BOND, InteractionRole.SIGMA_HOLE_DONOR): ("chalcogen bond donor", "chalcogen donor"),
    (InteractionType.CHALCOGEN_BOND, _ACCEPTOR): ("chalcogen bond acceptor", "chalcogen acceptor"),
    (InteractionType.METAL_COORDINATION, InteractionRole.COORDINATION_LIGAND): ("metal coordination ligand", "metal ligand"),
}
for interaction, roles in INTERACTION_SITES.items():
    for role in roles:
        shared_name = f"interaction {interaction.value.replace('_', ' ')} {role.value.replace('_', ' ')}"
        names = _INTERACTION_NAMES.get((interaction, role))
        if names is None:
            SELECTORS[shared_name] = _named_queries(_role_queries(interaction, role))
        else:
            SELECTORS[names[0]] = _named_queries(
                _role_queries(interaction, role), (*names[1:], shared_name)
            )

SELECTORS.update({
    "carbonyl": NamedSelector(CARBONYL_SITE_SMARTS),
    "positive charge": NamedSelector(
        POSITIVE_CHARGE_SITE_SMARTS,
        ("cation", "ionic positive", "cation pi cation", "ionic interaction cation"),
    ),
    "negative charge": NamedSelector(
        NEGATIVE_CHARGE_SITE_SMARTS, ("anion", "ionic negative", "ionic interaction anion")
    ),
    "metal center": NamedSelector(METAL_CENTER_SITE_SMARTS, ("metal coordination center",)),
})

CARBOXYLATE_SMARTS = CARBOXYLATE_OXYGEN_PAIR_SMARTS
CARBOXYLATE_QUERY = Chem.MolFromSmarts(CARBOXYLATE_SMARTS)


def carboxylate_sites(mol):
    """Map each oxygen in a carboxylate group to its carbon and oxygen indices."""
    groups = {}
    for carbon, oxygen1, oxygen2 in mol.GetSubstructMatches(CARBOXYLATE_QUERY, maxMatches=0):
        group = (carbon, *sorted((oxygen1, oxygen2)))
        groups[oxygen1] = group
        groups[oxygen2] = group
    return groups


SELECTORS["aromatic ring"] = NamedSelector(
    None, ("aromatic group", "pi pi", "pi pi ring", "pi pi interaction", "cation pi ring")
)


class SelectorRegistry:
    """Resolve canonical selector names and aliases."""

    def __init__(self, definitions=None):
        self.definitions = SELECTORS if definitions is None else definitions

    @staticmethod
    def normalize(name):
        return " ".join(re.sub(r"[-_]+", " ", name).casefold().split())

    def resolve(self, name):
        normalized = self.normalize(name)
        for canonical, definition in self.definitions.items():
            if normalized in (self.normalize(canonical), *map(self.normalize, definition.aliases)):
                return canonical, definition
        raise ValueError(f"unknown selector {name!r}; use --list-selectors to see available names")
