"""Named functional-group and interaction-site selector definitions."""

import json
import re
from dataclasses import dataclass

from molecular_functional_groups import (
    FunctionalGroup,
    INTERACTION_SITES,
    InteractionRole,
    InteractionType,
    MAPPED_FUNCTIONAL_GROUP_SMARTS,
    get_interaction_site,
)
from rdkit import Chem


@dataclass(frozen=True)
class SelectorQuery:
    smarts: str
    atom_map_numbers: tuple[int, ...] | None = None


@dataclass(frozen=True)
class NamedSelector:
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


CARBOXYLATE_SMARTS = "[C;+0;X3;$([C]-[#6]),$([C;H1])](=[O;+0;X1;H0:1])-[O-;X1;H0:2]"
CARBOXYLATE_QUERY = Chem.MolFromSmarts(CARBOXYLATE_SMARTS)


def carboxylate_sites(mol):
    """Map each oxygen in a carboxylate group to its carbon and oxygen indices."""
    groups = {}
    for carbon, oxygen1, oxygen2 in mol.GetSubstructMatches(CARBOXYLATE_QUERY, maxMatches=0):
        group = (carbon, *sorted((oxygen1, oxygen2)))
        groups[oxygen1] = group
        groups[oxygen2] = group
    return groups


SELECTORS = {
    group.value.replace("_", " "): _named_queries([_group_query(group)])
    for group in FunctionalGroup
}

_HBOND = InteractionType.HYDROGEN_BOND
_DONOR = InteractionRole.DONOR
_ACCEPTOR = InteractionRole.ACCEPTOR
_SITE_GROUPS = {
    "thioether": (FunctionalGroup.THIOETHER, _HBOND, _ACCEPTOR, ("sulfide",)),
    "thiol": (FunctionalGroup.THIOL, _HBOND, _DONOR, ()),
    "phenol": (FunctionalGroup.PHENOL, _HBOND, _DONOR, ()),
    "ether": (FunctionalGroup.ETHER, _HBOND, _ACCEPTOR, ()),
    "aldehyde": (FunctionalGroup.ALDEHYDE, _HBOND, _ACCEPTOR, ()),
    "ketone": (FunctionalGroup.KETONE, _HBOND, _ACCEPTOR, ()),
    "carboxylate": (FunctionalGroup.CARBOXYLATE, _HBOND, _ACCEPTOR, ("carboxylate anion",)),
    "ester": (FunctionalGroup.ESTER, _HBOND, _ACCEPTOR, ()),
}
for name, (group, interaction, role, aliases) in _SITE_GROUPS.items():
    SELECTORS[name] = _named_queries([_group_site(group, interaction, role)], aliases)

SELECTORS["carboxylic acid"] = _named_queries([
    _group_query(FunctionalGroup.CARBOXYLIC_ACID, (2, 3)),
])
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

for interaction, roles in INTERACTION_SITES.items():
    for role in roles:
        name = f"interaction {interaction.value.replace('_', ' ')} {role.value.replace('_', ' ')}"
        SELECTORS[name] = _named_queries(_role_queries(interaction, role))

SELECTORS.update({
    "carbonyl": NamedSelector("[C]=[O:1]"),
    "hydrogen bond donor": _named_queries(
        _role_queries(_HBOND, _DONOR), ("h bond donor", "hbond donor")
    ),
    "hydrogen bond acceptor": _named_queries(
        _role_queries(_HBOND, _ACCEPTOR), ("h bond acceptor", "hbond acceptor")
    ),
    "halogen bond donor": _named_queries(
        _role_queries(InteractionType.HALOGEN_BOND, InteractionRole.SIGMA_HOLE_DONOR),
        ("halogen donor",),
    ),
    "halogen bond acceptor": _named_queries(
        _role_queries(InteractionType.HALOGEN_BOND, _ACCEPTOR), ("halogen acceptor",)
    ),
    "positive charge": NamedSelector(
        "[+1,+2,+3,+4:1]", ("cation", "ionic positive", "cation pi cation", "ionic interaction cation")
    ),
    "negative charge": NamedSelector(
        "[-1,-2,-3,-4:1]", ("anion", "ionic negative", "ionic interaction anion")
    ),
    "aromatic ring": NamedSelector(
        None, ("aromatic group", "pi pi", "pi pi ring", "pi pi interaction", "cation pi ring")
    ),
    "chalcogen bond donor": _named_queries(
        _role_queries(InteractionType.CHALCOGEN_BOND, InteractionRole.SIGMA_HOLE_DONOR),
        ("chalcogen donor",),
    ),
    "chalcogen bond acceptor": _named_queries(
        _role_queries(InteractionType.CHALCOGEN_BOND, _ACCEPTOR), ("chalcogen acceptor",)
    ),
    "metal center": NamedSelector(
        "[Li,Na,K,Mg,Ca,Mn,Fe,Co,Ni,Cu,Zn,Cd,Hg,Pt,Pd;+1,+2,+3,+4:1]",
        ("metal coordination center",),
    ),
    "metal coordination ligand": _named_queries(
        _role_queries(InteractionType.METAL_COORDINATION, InteractionRole.COORDINATION_LIGAND),
        ("metal ligand",),
    ),
})


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
