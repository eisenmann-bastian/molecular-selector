"""Functional-group and atom-site SMARTS patterns used by the selector."""

from enum import Enum
from types import MappingProxyType
from typing import Final, Mapping


class FunctionalGroup(str, Enum):
    ALKANE_CARBON = "alkane_carbon"
    ALKENE = "alkene"
    ALKYNE = "alkyne"
    ARENE = "arene"
    ALKYL_HALIDE = "alkyl_halide"
    ARYL_HALIDE = "aryl_halide"
    PRIMARY_ALCOHOL = "primary_alcohol"
    SECONDARY_ALCOHOL = "secondary_alcohol"
    TERTIARY_ALCOHOL = "tertiary_alcohol"
    PHENOL = "phenol"
    ETHER = "ether"
    ALDEHYDE = "aldehyde"
    KETONE = "ketone"
    ACYL_HALIDE = "acyl_halide"
    CARBOXYLIC_ACID = "carboxylic_acid"
    CARBOXYLATE = "carboxylate"
    ESTER = "ester"
    ACID_ANHYDRIDE = "acid_anhydride"
    PRIMARY_AMINE = "primary_amine"
    SECONDARY_AMINE = "secondary_amine"
    TERTIARY_AMINE = "tertiary_amine"
    AMMONIUM = "ammonium"
    PRIMARY_AMIDE = "primary_amide"
    SECONDARY_AMIDE = "secondary_amide"
    TERTIARY_AMIDE = "tertiary_amide"
    NITRILE = "nitrile"
    IMINE = "imine"
    ENOL = "enol"
    NITRO = "nitro"
    ALPHA_AMINO_ACID = "alpha_amino_acid"
    THIOL = "thiol"
    THIOETHER = "thioether"
    DISULFIDE = "disulfide"
    SULFOXIDE = "sulfoxide"
    SULFONE = "sulfone"
    SULFONAMIDE = "sulfonamide"
    SULFONAMIDE_NH = "sulfonamide_nh"
    SULFONIC_ACID = "sulfonic_acid"
    PHOSPHATE_ESTER = "phosphate_ester"
    PHOSPHONATE = "phosphonate"


ALKANE_CARBON_SMARTS: Final = "[CX4]"
ALKENE_SMARTS: Final = "[CX3]=[CX3]"
ALKYNE_SMARTS: Final = "[CX2]#[CX2]"
ARENE_SMARTS: Final = "c1ccccc1"
ALKYL_HALIDE_SMARTS: Final = "[CX4]-[F,Cl,Br,I]"
ARYL_HALIDE_SMARTS: Final = "[c]-[F,Cl,Br,I]"
PRIMARY_ALCOHOL_SMARTS: Final = "[CX4H2]-[OX2H]"
SECONDARY_ALCOHOL_SMARTS: Final = "[CX4H1]([#6])-[OX2H]"
TERTIARY_ALCOHOL_SMARTS: Final = "[CX4H0]([#6])([#6])-[OX2H]"
PHENOL_SMARTS: Final = "[c]-[OX2H]"
ETHER_SMARTS: Final = "[#6;!$(C=O)]-[OX2;H0]-[#6;!$(C=O)]"
ALDEHYDE_SMARTS: Final = "[CX3H1](=O)[#6]"
KETONE_SMARTS: Final = "[CX3](=O)([#6])[#6]"
ACYL_HALIDE_SMARTS: Final = "[CX3](=O)[F,Cl,Br,I]"
CARBOXYLIC_ACID_SMARTS: Final = "[CX3](=O)[OX2H1]"
CARBOXYLATE_SMARTS: Final = "[CX3](=O)[O-]"
ESTER_SMARTS: Final = "[CX3](=O)[OX2H0][#6]"
ACID_ANHYDRIDE_SMARTS: Final = "[CX3](=O)[OX2][CX3](=O)"
PRIMARY_AMINE_SMARTS: Final = "[NX3;H2;!$(N-C=O)]"
SECONDARY_AMINE_SMARTS: Final = "[NX3;H1;!$(N-C=O)]"
TERTIARY_AMINE_SMARTS: Final = "[NX3;H0;!$(N-C=O)]"
AMMONIUM_SMARTS: Final = "[NX4+]"
PRIMARY_AMIDE_SMARTS: Final = "[CX3](=O)[NX3;H2]"
SECONDARY_AMIDE_SMARTS: Final = "[CX3](=O)[NX3;H1]"
TERTIARY_AMIDE_SMARTS: Final = "[CX3](=O)[NX3;H0]"
NITRILE_SMARTS: Final = "[C]#[N]"
IMINE_SMARTS: Final = "[CX3]=[NX2]"
ENOL_SMARTS: Final = "[CX3]=[CX3]-[OX2H]"
NITRO_SMARTS: Final = "[N+](=O)[O-]"
ALPHA_AMINO_ACID_SMARTS: Final = "[NX3;H2,H1]-[CX4]-[CX3](=O)[OX2H1]"
THIOL_SMARTS: Final = "[SX2H]"
THIOETHER_SMARTS: Final = "[SX2]([#6])[#6]"
DISULFIDE_SMARTS: Final = "[SX2]-[SX2]"
SULFOXIDE_SMARTS: Final = "[SX3](=[OX1])([#6])[#6]"
SULFONE_SMARTS: Final = "[SX4](=[OX1])(=[OX1])([#6])[#6]"
SULFONAMIDE_SMARTS: Final = "[SX4](=[OX1])(=[OX1])[NX3]"
SULFONAMIDE_NH_SMARTS: Final = "[SX4](=[OX1])(=[OX1])[NX3;H1,H2]"
SULFONIC_ACID_SMARTS: Final = "[SX4](=[OX1])(=[OX1])[OX2H1]"
PHOSPHATE_ESTER_SMARTS: Final = "[PX4](=[OX1])([OX2])([OX2])[OX2][#6]"
PHOSPHONATE_SMARTS: Final = "[PX4](=[OX1])([OX2])([OX2])[#6]"

# General atom-site queries; these are broader than individual group patterns.
CARBONYL_SITE_SMARTS: Final = "[C]=[O:1]"
POSITIVE_CHARGE_SITE_SMARTS: Final = "[+1,+2,+3,+4:1]"
NEGATIVE_CHARGE_SITE_SMARTS: Final = "[-1,-2,-3,-4:1]"
METAL_CENTER_SITE_SMARTS: Final = (
    "[Li,Na,K,Mg,Ca,Mn,Fe,Co,Ni,Cu,Zn,Cd,Hg,Pt,Pd;+1,+2,+3,+4:1]"
)
CARBOXYLATE_OXYGEN_PAIR_SMARTS: Final = (
    "[C;+0;X3;$([C]-[#6]),$([C;H1])](=[O;+0;X1;H0:1])-[O-;X1;H0:2]"
)

FUNCTIONAL_GROUP_SMARTS: Final[Mapping[FunctionalGroup, str]] = MappingProxyType(
    {
        FunctionalGroup.ALKANE_CARBON: ALKANE_CARBON_SMARTS,
        FunctionalGroup.ALKENE: ALKENE_SMARTS,
        FunctionalGroup.ALKYNE: ALKYNE_SMARTS,
        FunctionalGroup.ARENE: ARENE_SMARTS,
        FunctionalGroup.ALKYL_HALIDE: ALKYL_HALIDE_SMARTS,
        FunctionalGroup.ARYL_HALIDE: ARYL_HALIDE_SMARTS,
        FunctionalGroup.PRIMARY_ALCOHOL: PRIMARY_ALCOHOL_SMARTS,
        FunctionalGroup.SECONDARY_ALCOHOL: SECONDARY_ALCOHOL_SMARTS,
        FunctionalGroup.TERTIARY_ALCOHOL: TERTIARY_ALCOHOL_SMARTS,
        FunctionalGroup.PHENOL: PHENOL_SMARTS,
        FunctionalGroup.ETHER: ETHER_SMARTS,
        FunctionalGroup.ALDEHYDE: ALDEHYDE_SMARTS,
        FunctionalGroup.KETONE: KETONE_SMARTS,
        FunctionalGroup.ACYL_HALIDE: ACYL_HALIDE_SMARTS,
        FunctionalGroup.CARBOXYLIC_ACID: CARBOXYLIC_ACID_SMARTS,
        FunctionalGroup.CARBOXYLATE: CARBOXYLATE_SMARTS,
        FunctionalGroup.ESTER: ESTER_SMARTS,
        FunctionalGroup.ACID_ANHYDRIDE: ACID_ANHYDRIDE_SMARTS,
        FunctionalGroup.PRIMARY_AMINE: PRIMARY_AMINE_SMARTS,
        FunctionalGroup.SECONDARY_AMINE: SECONDARY_AMINE_SMARTS,
        FunctionalGroup.TERTIARY_AMINE: TERTIARY_AMINE_SMARTS,
        FunctionalGroup.AMMONIUM: AMMONIUM_SMARTS,
        FunctionalGroup.PRIMARY_AMIDE: PRIMARY_AMIDE_SMARTS,
        FunctionalGroup.SECONDARY_AMIDE: SECONDARY_AMIDE_SMARTS,
        FunctionalGroup.TERTIARY_AMIDE: TERTIARY_AMIDE_SMARTS,
        FunctionalGroup.NITRILE: NITRILE_SMARTS,
        FunctionalGroup.IMINE: IMINE_SMARTS,
        FunctionalGroup.ENOL: ENOL_SMARTS,
        FunctionalGroup.NITRO: NITRO_SMARTS,
        FunctionalGroup.ALPHA_AMINO_ACID: ALPHA_AMINO_ACID_SMARTS,
        FunctionalGroup.THIOL: THIOL_SMARTS,
        FunctionalGroup.THIOETHER: THIOETHER_SMARTS,
        FunctionalGroup.DISULFIDE: DISULFIDE_SMARTS,
        FunctionalGroup.SULFOXIDE: SULFOXIDE_SMARTS,
        FunctionalGroup.SULFONE: SULFONE_SMARTS,
        FunctionalGroup.SULFONAMIDE: SULFONAMIDE_SMARTS,
        FunctionalGroup.SULFONAMIDE_NH: SULFONAMIDE_NH_SMARTS,
        FunctionalGroup.SULFONIC_ACID: SULFONIC_ACID_SMARTS,
        FunctionalGroup.PHOSPHATE_ESTER: PHOSPHATE_ESTER_SMARTS,
        FunctionalGroup.PHOSPHONATE: PHOSPHONATE_SMARTS,
    }
)


def get_smarts_pattern(group: FunctionalGroup | str) -> str:
    """Return the SMARTS pattern for a functional-group enum or identifier."""
    return FUNCTIONAL_GROUP_SMARTS[FunctionalGroup(group)]


__all__ = [
    "FunctionalGroup",
    "FUNCTIONAL_GROUP_SMARTS",
    "get_smarts_pattern",
    *(name for name in globals() if name.endswith("_SMARTS")),
]
