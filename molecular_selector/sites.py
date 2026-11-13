"""Atom-mapped group patterns and potential interaction sites."""

import re
from dataclasses import dataclass
from types import MappingProxyType
from typing import Final, Mapping

from .interactions import INTERACTION_ROLES, InteractionRole, InteractionType
from .patterns import FunctionalGroup


# Labels are local to a pattern; they are not molecular atom indices.
MAPPED_FUNCTIONAL_GROUP_SMARTS: Final[Mapping[FunctionalGroup, str]] = MappingProxyType(
    {
        FunctionalGroup.ALKANE_CARBON: "[CX4:1]",
        FunctionalGroup.ALKENE: "[CX3:1]=[CX3:2]",
        FunctionalGroup.ALKYNE: "[CX2:1]#[CX2:2]",
        FunctionalGroup.ARENE: "[c:1]1[c:2][c:3][c:4][c:5][c:6]1",
        FunctionalGroup.ALKYL_HALIDE: "[CX4:1]-[F,Cl,Br,I:2]",
        FunctionalGroup.ARYL_HALIDE: "[c:1]-[F,Cl,Br,I:2]",
        FunctionalGroup.PRIMARY_ALCOHOL: "[CX4H2:1]-[OX2H:2]",
        FunctionalGroup.SECONDARY_ALCOHOL: "[CX4H1:1]([#6:2])-[OX2H:3]",
        FunctionalGroup.TERTIARY_ALCOHOL: "[CX4H0:1]([#6:2])([#6:3])-[OX2H:4]",
        FunctionalGroup.PHENOL: "[c:1]-[OX2H:2]",
        FunctionalGroup.ETHER: "[#6;!$(C=O):1]-[OX2;H0:2]-[#6;!$(C=O):3]",
        FunctionalGroup.ALDEHYDE: "[CX3H1:1](=[O:2])[#6:3]",
        FunctionalGroup.KETONE: "[CX3:1](=[O:2])([#6:3])[#6:4]",
        FunctionalGroup.ACYL_HALIDE: "[CX3:1](=[O:2])[F,Cl,Br,I:3]",
        FunctionalGroup.CARBOXYLIC_ACID: "[CX3:1](=[O:2])[OX2H1:3]",
        FunctionalGroup.CARBOXYLATE: "[CX3:1](=[O:2])[O-:3]",
        FunctionalGroup.ESTER: "[CX3:1](=[O:2])[OX2H0:3][#6:4]",
        FunctionalGroup.ACID_ANHYDRIDE: "[CX3:1](=[O:2])[OX2:3][CX3:4](=[O:5])",
        FunctionalGroup.PRIMARY_AMINE: "[NX3;H2;!$(N-C=O):1]",
        FunctionalGroup.SECONDARY_AMINE: "[NX3;H1;!$(N-C=O):1]",
        FunctionalGroup.TERTIARY_AMINE: "[NX3;H0;!$(N-C=O):1]",
        FunctionalGroup.AMMONIUM: "[NX4+:1]",
        FunctionalGroup.PRIMARY_AMIDE: "[CX3:1](=[O:2])[NX3;H2:3]",
        FunctionalGroup.SECONDARY_AMIDE: "[CX3:1](=[O:2])[NX3;H1:3]",
        FunctionalGroup.TERTIARY_AMIDE: "[CX3:1](=[O:2])[NX3;H0:3]",
        FunctionalGroup.NITRILE: "[C:1]#[N:2]",
        FunctionalGroup.IMINE: "[CX3:1]=[NX2:2]",
        FunctionalGroup.ENOL: "[CX3:1]=[CX3:2]-[OX2H:3]",
        FunctionalGroup.NITRO: "[N+:1](=[O:2])[O-:3]",
        FunctionalGroup.ALPHA_AMINO_ACID: "[NX3;H2,H1:1]-[CX4:2]-[CX3:3](=[O:4])[OX2H1:5]",
        FunctionalGroup.THIOL: "[SX2H:1]",
        FunctionalGroup.THIOETHER: "[SX2:1]([#6:2])[#6:3]",
        FunctionalGroup.DISULFIDE: "[SX2:1]-[SX2:2]",
        FunctionalGroup.SULFOXIDE: "[SX3:1](=[OX1:2])([#6:3])[#6:4]",
        FunctionalGroup.SULFONE: "[SX4:1](=[OX1:2])(=[OX1:3])([#6:4])[#6:5]",
        FunctionalGroup.SULFONAMIDE: "[SX4:1](=[OX1:2])(=[OX1:3])[NX3:4]",
        FunctionalGroup.SULFONAMIDE_NH: "[SX4:1](=[OX1:2])(=[OX1:3])[NX3;H1,H2:4]",
        FunctionalGroup.SULFONIC_ACID: "[SX4:1](=[OX1:2])(=[OX1:3])[OX2H1:4]",
        FunctionalGroup.PHOSPHATE_ESTER: "[PX4:1](=[OX1:2])([OX2:3])([OX2:4])[OX2:5][#6:6]",
        FunctionalGroup.PHOSPHONATE: "[PX4:1](=[OX1:2])([OX2:3])([OX2:4])[#6:5]",
    }
)


@dataclass(frozen=True)
class InteractionSite:
    """SMARTS and local atom-map labels for a potential interaction site."""

    smarts: str
    atom_map_numbers: tuple[int, ...]


_DONORS = {
    FunctionalGroup.PRIMARY_ALCOHOL: (2,),
    FunctionalGroup.SECONDARY_ALCOHOL: (3,),
    FunctionalGroup.TERTIARY_ALCOHOL: (4,),
    FunctionalGroup.PHENOL: (2,),
    FunctionalGroup.CARBOXYLIC_ACID: (3,),
    FunctionalGroup.PRIMARY_AMINE: (1,),
    FunctionalGroup.SECONDARY_AMINE: (1,),
    FunctionalGroup.AMMONIUM: (1,),
    FunctionalGroup.PRIMARY_AMIDE: (3,),
    FunctionalGroup.SECONDARY_AMIDE: (3,),
    FunctionalGroup.SULFONAMIDE_NH: (4,),
    FunctionalGroup.SULFONIC_ACID: (4,),
    FunctionalGroup.ALPHA_AMINO_ACID: (1, 5),
    FunctionalGroup.THIOL: (1,),
}

_ACCEPTORS = {
    FunctionalGroup.PRIMARY_ALCOHOL: (2,),
    FunctionalGroup.SECONDARY_ALCOHOL: (3,),
    FunctionalGroup.TERTIARY_ALCOHOL: (4,),
    FunctionalGroup.PHENOL: (2,),
    FunctionalGroup.ETHER: (2,),
    FunctionalGroup.ALDEHYDE: (2,),
    FunctionalGroup.KETONE: (2,),
    FunctionalGroup.CARBOXYLIC_ACID: (2,),
    FunctionalGroup.CARBOXYLATE: (2, 3),
    FunctionalGroup.ESTER: (2, 3),
    FunctionalGroup.ACID_ANHYDRIDE: (2, 5),
    FunctionalGroup.PRIMARY_AMINE: (1,),
    FunctionalGroup.SECONDARY_AMINE: (1,),
    FunctionalGroup.TERTIARY_AMINE: (1,),
    FunctionalGroup.PRIMARY_AMIDE: (2,),
    FunctionalGroup.SECONDARY_AMIDE: (2,),
    FunctionalGroup.TERTIARY_AMIDE: (2,),
    FunctionalGroup.NITRILE: (2,),
    FunctionalGroup.IMINE: (2,),
    FunctionalGroup.ENOL: (3,),
    FunctionalGroup.NITRO: (2, 3),
    FunctionalGroup.THIOETHER: (1,),
    FunctionalGroup.DISULFIDE: (1, 2),
    FunctionalGroup.SULFOXIDE: (2,),
    FunctionalGroup.SULFONE: (2, 3),
    FunctionalGroup.SULFONAMIDE: (2, 3),
    FunctionalGroup.SULFONIC_ACID: (2, 3),
    FunctionalGroup.PHOSPHATE_ESTER: (2, 3, 4, 5),
    FunctionalGroup.PHOSPHONATE: (2, 3, 4),
    FunctionalGroup.ALPHA_AMINO_ACID: (1, 4),
}

_LIGANDS = {
    **_ACCEPTORS,
    FunctionalGroup.THIOL: (1,),
}

_SIGMA_HOLES = {
    FunctionalGroup.ALKYL_HALIDE: (2,),
    FunctionalGroup.ARYL_HALIDE: (2,),
    FunctionalGroup.THIOETHER: (1,),
    FunctionalGroup.DISULFIDE: (1, 2),
    FunctionalGroup.SULFOXIDE: (1,),
    FunctionalGroup.SULFONE: (1,),
    FunctionalGroup.ALDEHYDE: (1,),
    FunctionalGroup.KETONE: (1,),
    FunctionalGroup.ACYL_HALIDE: (1,),
    FunctionalGroup.CARBOXYLIC_ACID: (1,),
    FunctionalGroup.CARBOXYLATE: (1,),
    FunctionalGroup.ESTER: (1,),
    FunctionalGroup.ACID_ANHYDRIDE: (1, 4),
    FunctionalGroup.PHOSPHATE_ESTER: (1,),
    FunctionalGroup.PHOSPHONATE: (1,),
}

_NONPOLAR = {
    FunctionalGroup.ALKANE_CARBON: (1,),
    FunctionalGroup.ALKENE: (1, 2),
    FunctionalGroup.ALKYNE: (1, 2),
    FunctionalGroup.ARENE: (1, 2, 3, 4, 5, 6),
    FunctionalGroup.ALKYL_HALIDE: (1, 2),
    FunctionalGroup.ARYL_HALIDE: (1, 2),
}


def _all_labels(group: FunctionalGroup) -> tuple[int, ...]:
    return tuple(map(int, re.findall(r":(\d+)\]", MAPPED_FUNCTIONAL_GROUP_SMARTS[group])))


def _site(role: InteractionRole, group: FunctionalGroup) -> InteractionSite:
    smarts = MAPPED_FUNCTIONAL_GROUP_SMARTS[group]
    if role == InteractionRole.DONOR:
        labels = _DONORS[group]
        if group == FunctionalGroup.AMMONIUM:
            smarts = "[NX4+;!H0:1]"
    elif role == InteractionRole.ACCEPTOR:
        labels = _ACCEPTORS[group]
    elif role == InteractionRole.COORDINATION_LIGAND:
        labels = _LIGANDS[group]
    elif role == InteractionRole.SIGMA_HOLE_DONOR:
        labels = _SIGMA_HOLES[group]
    elif role == InteractionRole.NONPOLAR_SURFACE:
        labels = _NONPOLAR[group]
    elif role == InteractionRole.CATION:
        labels = (1,)
    elif role == InteractionRole.ANION:
        labels = (2, 3)
    elif role in (
        InteractionRole.DIPOLE_PARTNER,
        InteractionRole.PI_SURFACE,
        InteractionRole.CONTACT_PARTNER,
    ):
        labels = _all_labels(group)
    else:
        raise ValueError(f"No site-label definition for interaction role {role!r}")
    return InteractionSite(smarts, labels)


INTERACTION_SITES: Final[
    Mapping[InteractionType, Mapping[InteractionRole, Mapping[FunctionalGroup, InteractionSite]]]
] = MappingProxyType(
    {
        interaction: MappingProxyType(
            {
                role: MappingProxyType(
                    {group: _site(role, group) for group in groups}
                )
                for role, groups in roles.items()
            }
        )
        for interaction, roles in INTERACTION_ROLES.items()
    }
)


def get_interaction_site(
    interaction: InteractionType | str,
    role: InteractionRole | str,
    group: FunctionalGroup | str,
) -> InteractionSite:
    """Return a site's query and map labels; raise KeyError for unmapped roles/groups."""
    interaction_type = InteractionType(interaction)
    interaction_role = InteractionRole(role)
    functional_group = FunctionalGroup(group)
    return INTERACTION_SITES[interaction_type][interaction_role][functional_group]
