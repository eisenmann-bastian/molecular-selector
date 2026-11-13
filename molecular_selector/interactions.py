"""Functional-group roles in non-covalent interactions."""

from enum import Enum
from types import MappingProxyType
from typing import Final, Mapping

from .patterns import FunctionalGroup


class InteractionType(str, Enum):
    HYDROGEN_BOND = "hydrogen_bond"
    IONIC = "ionic"
    DIPOLE_DIPOLE = "dipole_dipole"
    PI_PI = "pi_pi"
    CATION_PI = "cation_pi"
    METAL_COORDINATION = "metal_coordination"
    HALOGEN_BOND = "halogen_bond"
    CHALCOGEN_BOND = "chalcogen_bond"
    TETREL_BOND = "tetrel_bond"
    PNICTOGEN_BOND = "pnictogen_bond"
    VAN_DER_WAALS = "van_der_waals"
    HYDROPHOBIC_EFFECT = "hydrophobic_effect"


class InteractionRole(str, Enum):
    DONOR = "donor"
    ACCEPTOR = "acceptor"
    CATION = "cation"
    ANION = "anion"
    PI_SURFACE = "pi_surface"
    COORDINATION_LIGAND = "coordination_ligand"
    SIGMA_HOLE_DONOR = "sigma_hole_donor"
    DIPOLE_PARTNER = "dipole_partner"
    CONTACT_PARTNER = "contact_partner"
    NONPOLAR_SURFACE = "nonpolar_surface"


_HBOND_DONORS = frozenset(
    {
        FunctionalGroup.PRIMARY_ALCOHOL,
        FunctionalGroup.SECONDARY_ALCOHOL,
        FunctionalGroup.TERTIARY_ALCOHOL,
        FunctionalGroup.PHENOL,
        FunctionalGroup.CARBOXYLIC_ACID,
        FunctionalGroup.PRIMARY_AMINE,
        FunctionalGroup.SECONDARY_AMINE,
        FunctionalGroup.AMMONIUM,
        FunctionalGroup.PRIMARY_AMIDE,
        FunctionalGroup.SECONDARY_AMIDE,
        FunctionalGroup.SULFONAMIDE_NH,
        FunctionalGroup.SULFONIC_ACID,
        FunctionalGroup.ALPHA_AMINO_ACID,
        FunctionalGroup.THIOL,
    }
)

_HBOND_ACCEPTORS = frozenset(
    {
        FunctionalGroup.PRIMARY_ALCOHOL,
        FunctionalGroup.SECONDARY_ALCOHOL,
        FunctionalGroup.TERTIARY_ALCOHOL,
        FunctionalGroup.PHENOL,
        FunctionalGroup.ETHER,
        FunctionalGroup.ALDEHYDE,
        FunctionalGroup.KETONE,
        FunctionalGroup.CARBOXYLIC_ACID,
        FunctionalGroup.CARBOXYLATE,
        FunctionalGroup.ESTER,
        FunctionalGroup.ACID_ANHYDRIDE,
        FunctionalGroup.PRIMARY_AMINE,
        FunctionalGroup.SECONDARY_AMINE,
        FunctionalGroup.TERTIARY_AMINE,
        FunctionalGroup.PRIMARY_AMIDE,
        FunctionalGroup.SECONDARY_AMIDE,
        FunctionalGroup.TERTIARY_AMIDE,
        FunctionalGroup.NITRILE,
        FunctionalGroup.IMINE,
        FunctionalGroup.ENOL,
        FunctionalGroup.NITRO,
        FunctionalGroup.THIOETHER,
        FunctionalGroup.DISULFIDE,
        FunctionalGroup.SULFOXIDE,
        FunctionalGroup.SULFONE,
        FunctionalGroup.SULFONAMIDE,
        FunctionalGroup.SULFONIC_ACID,
        FunctionalGroup.PHOSPHATE_ESTER,
        FunctionalGroup.PHOSPHONATE,
        FunctionalGroup.ALPHA_AMINO_ACID,
    }
)

_LEWIS_BASES = _HBOND_ACCEPTORS - {
    FunctionalGroup.CARBOXYLIC_ACID,
    FunctionalGroup.PHENOL,
    FunctionalGroup.ENOL,
}

_POLAR_GROUPS = frozenset(
    {
        FunctionalGroup.ALKYL_HALIDE,
        FunctionalGroup.ARYL_HALIDE,
        FunctionalGroup.PRIMARY_ALCOHOL,
        FunctionalGroup.SECONDARY_ALCOHOL,
        FunctionalGroup.TERTIARY_ALCOHOL,
        FunctionalGroup.PHENOL,
        FunctionalGroup.ETHER,
        FunctionalGroup.ALDEHYDE,
        FunctionalGroup.KETONE,
        FunctionalGroup.ACYL_HALIDE,
        FunctionalGroup.CARBOXYLIC_ACID,
        FunctionalGroup.CARBOXYLATE,
        FunctionalGroup.ESTER,
        FunctionalGroup.ACID_ANHYDRIDE,
        FunctionalGroup.PRIMARY_AMINE,
        FunctionalGroup.SECONDARY_AMINE,
        FunctionalGroup.TERTIARY_AMINE,
        FunctionalGroup.AMMONIUM,
        FunctionalGroup.PRIMARY_AMIDE,
        FunctionalGroup.SECONDARY_AMIDE,
        FunctionalGroup.TERTIARY_AMIDE,
        FunctionalGroup.NITRILE,
        FunctionalGroup.IMINE,
        FunctionalGroup.NITRO,
        FunctionalGroup.SULFOXIDE,
        FunctionalGroup.SULFONE,
        FunctionalGroup.SULFONAMIDE,
        FunctionalGroup.SULFONIC_ACID,
        FunctionalGroup.PHOSPHATE_ESTER,
        FunctionalGroup.PHOSPHONATE,
    }
)

_COORDINATION_LIGANDS = frozenset(
    {
        FunctionalGroup.CARBOXYLIC_ACID,
        FunctionalGroup.CARBOXYLATE,
        FunctionalGroup.ETHER,
        FunctionalGroup.ALDEHYDE,
        FunctionalGroup.KETONE,
        FunctionalGroup.ESTER,
        FunctionalGroup.ACID_ANHYDRIDE,
        FunctionalGroup.PRIMARY_AMINE,
        FunctionalGroup.SECONDARY_AMINE,
        FunctionalGroup.TERTIARY_AMINE,
        FunctionalGroup.IMINE,
        FunctionalGroup.NITRILE,
        FunctionalGroup.THIOL,
        FunctionalGroup.THIOETHER,
        FunctionalGroup.DISULFIDE,
        FunctionalGroup.PHENOL,
        FunctionalGroup.SULFOXIDE,
        FunctionalGroup.SULFONE,
        FunctionalGroup.PHOSPHATE_ESTER,
        FunctionalGroup.PHOSPHONATE,
    }
)

_ALL_GROUPS = frozenset(FunctionalGroup)
_HYDROPHOBIC_GROUPS = frozenset(
    {
        FunctionalGroup.ALKANE_CARBON,
        FunctionalGroup.ALKENE,
        FunctionalGroup.ALKYNE,
        FunctionalGroup.ARENE,
        FunctionalGroup.ALKYL_HALIDE,
        FunctionalGroup.ARYL_HALIDE,
    }
)


def _role_map(
    **roles: frozenset[FunctionalGroup],
) -> Mapping[InteractionRole, frozenset[FunctionalGroup]]:
    return MappingProxyType(
        {InteractionRole(role): groups for role, groups in roles.items()}
    )


INTERACTION_ROLES: Final[
    Mapping[InteractionType, Mapping[InteractionRole, frozenset[FunctionalGroup]]]
] = MappingProxyType(
    {
        InteractionType.HYDROGEN_BOND: _role_map(
            donor=_HBOND_DONORS,
            acceptor=_HBOND_ACCEPTORS,
        ),
        InteractionType.IONIC: _role_map(
            cation=frozenset({FunctionalGroup.AMMONIUM}),
            anion=frozenset({FunctionalGroup.CARBOXYLATE}),
        ),
        InteractionType.DIPOLE_DIPOLE: _role_map(
            dipole_partner=_POLAR_GROUPS,
        ),
        InteractionType.PI_PI: _role_map(
            pi_surface=frozenset({FunctionalGroup.ARENE}),
        ),
        InteractionType.CATION_PI: _role_map(
            cation=frozenset({FunctionalGroup.AMMONIUM}),
            pi_surface=frozenset({FunctionalGroup.ARENE}),
        ),
        InteractionType.METAL_COORDINATION: _role_map(
            coordination_ligand=_COORDINATION_LIGANDS,
        ),
        InteractionType.HALOGEN_BOND: _role_map(
            sigma_hole_donor=frozenset(
                {FunctionalGroup.ALKYL_HALIDE, FunctionalGroup.ARYL_HALIDE}
            ),
            acceptor=_LEWIS_BASES,
        ),
        InteractionType.CHALCOGEN_BOND: _role_map(
            sigma_hole_donor=frozenset(
                {
                    FunctionalGroup.THIOETHER,
                    FunctionalGroup.DISULFIDE,
                    FunctionalGroup.SULFOXIDE,
                    FunctionalGroup.SULFONE,
                }
            ),
            acceptor=_LEWIS_BASES,
        ),
        InteractionType.TETREL_BOND: _role_map(
            sigma_hole_donor=frozenset(
                {
                    FunctionalGroup.ALDEHYDE,
                    FunctionalGroup.KETONE,
                    FunctionalGroup.ACYL_HALIDE,
                    FunctionalGroup.CARBOXYLIC_ACID,
                    FunctionalGroup.CARBOXYLATE,
                    FunctionalGroup.ESTER,
                    FunctionalGroup.ACID_ANHYDRIDE,
                }
            ),
            acceptor=_LEWIS_BASES,
        ),
        InteractionType.PNICTOGEN_BOND: _role_map(
            sigma_hole_donor=frozenset(
                {FunctionalGroup.PHOSPHATE_ESTER, FunctionalGroup.PHOSPHONATE}
            ),
            acceptor=_LEWIS_BASES,
        ),
        InteractionType.VAN_DER_WAALS: _role_map(contact_partner=_ALL_GROUPS),
        InteractionType.HYDROPHOBIC_EFFECT: _role_map(
            nonpolar_surface=_HYDROPHOBIC_GROUPS,
        ),
    }
)

BOND_TYPE_ROLES = INTERACTION_ROLES


def get_groups_for_role(
    interaction: InteractionType | str,
    role: InteractionRole | str,
) -> frozenset[FunctionalGroup]:
    """Return groups mapped to an interaction role, or an empty set if absent."""
    interaction_type = InteractionType(interaction)
    interaction_role = InteractionRole(role)
    return INTERACTION_ROLES.get(interaction_type, {}).get(interaction_role, frozenset())


__all__ = [
    "InteractionType",
    "InteractionRole",
    "INTERACTION_ROLES",
    "BOND_TYPE_ROLES",
    "get_groups_for_role",
]
