import json
import unittest

from molecular_selector import (
    AMMONIUM_SMARTS,
    BOND_TYPE_ROLES,
    CARBONYL_SITE_SMARTS,
    CARBOXYLATE_OXYGEN_PAIR_SMARTS,
    CARBOXYLATE_SMARTS,
    FunctionalGroup,
    FUNCTIONAL_GROUP_SMARTS,
    INTERACTION_ROLES,
    INTERACTION_SITES,
    InteractionRole,
    InteractionType,
    MAPPED_FUNCTIONAL_GROUP_SMARTS,
    METAL_CENTER_SITE_SMARTS,
    NEGATIVE_CHARGE_SITE_SMARTS,
    POSITIVE_CHARGE_SITE_SMARTS,
    get_groups_for_role,
    get_interaction_site,
)
from molecular_selector.definitions import SELECTORS
from molecular_selector import interactions, patterns, sites


class FunctionalGroupPatternsTests(unittest.TestCase):
    def test_public_definitions_are_used_directly_by_the_selector(self) -> None:
        self.assertIs(FunctionalGroup, patterns.FunctionalGroup)
        self.assertIs(INTERACTION_ROLES, interactions.INTERACTION_ROLES)
        self.assertIs(INTERACTION_SITES, sites.INTERACTION_SITES)
        self.assertIs(MAPPED_FUNCTIONAL_GROUP_SMARTS, sites.MAPPED_FUNCTIONAL_GROUP_SMARTS)
        self.assertEqual(
            SELECTORS["phenol"].queries[0].smarts,
            MAPPED_FUNCTIONAL_GROUP_SMARTS[FunctionalGroup.PHENOL],
        )

    def test_patterns_cover_every_functional_group(self) -> None:
        self.assertEqual(set(FUNCTIONAL_GROUP_SMARTS), set(FunctionalGroup))
        self.assertTrue(all(FUNCTIONAL_GROUP_SMARTS.values()))

    def test_named_constants_are_exported(self) -> None:
        self.assertEqual(CARBOXYLATE_SMARTS, FUNCTIONAL_GROUP_SMARTS[FunctionalGroup.CARBOXYLATE])
        self.assertEqual(AMMONIUM_SMARTS, "[NX4+]")

    def test_general_site_queries(self) -> None:
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("RDKit is optional")
        cases = (
            (CARBONYL_SITE_SMARTS, "CC=O", "CCO", (1,)),
            (POSITIVE_CHARGE_SITE_SMARTS, "C[NH3+]", "CN", (1,)),
            (NEGATIVE_CHARGE_SITE_SMARTS, "CC(=O)[O-]", "CC(=O)O", (1,)),
            (METAL_CENTER_SITE_SMARTS, "[Zn+2]", "C[NH3+]", (1,)),
            (CARBOXYLATE_OXYGEN_PAIR_SMARTS, "CC(=O)[O-]", "C(=O)([O-])OC", (1, 2)),
        )
        for smarts, positive, negative, labels in cases:
            with self.subTest(smarts=smarts):
                query = Chem.MolFromSmarts(smarts)
                self.assertIsNotNone(query)
                self.assertEqual(tuple(a.GetAtomMapNum() for a in query.GetAtoms()
                                       if a.GetAtomMapNum()), labels)
                self.assertTrue(Chem.MolFromSmiles(positive).HasSubstructMatch(query))
                self.assertFalse(Chem.MolFromSmiles(negative).HasSubstructMatch(query))

    def test_rdkit_parses_all_patterns_when_available(self) -> None:
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("RDKit is optional")

        invalid = [
            (group.value, smarts)
            for group, smarts in FUNCTIONAL_GROUP_SMARTS.items()
            if Chem.MolFromSmarts(smarts) is None
        ]
        self.assertEqual(invalid, [])

        examples = (
            (FunctionalGroup.CARBOXYLATE, "CC(=O)[O-]"),
            (FunctionalGroup.PRIMARY_AMIDE, "CC(N)=O"),
            (FunctionalGroup.ARENE, "c1ccccc1"),
            (FunctionalGroup.PHENOL, "Oc1ccccc1"),
            (FunctionalGroup.PHOSPHONATE, "CP(=O)(O)O"),
        )
        for group, smiles in examples:
            with self.subTest(group=group):
                molecule = Chem.MolFromSmiles(smiles)
                query = Chem.MolFromSmarts(FUNCTIONAL_GROUP_SMARTS[group])
                self.assertIsNotNone(molecule)
                self.assertIsNotNone(query)
                self.assertTrue(molecule.HasSubstructMatch(query))


class NamedSelectorTests(unittest.TestCase):
    def test_group_and_interaction_selectors_are_derived_from_shared_sites(self) -> None:
        for group in FunctionalGroup:
            selector = SELECTORS[group.value.replace("_", " ")]
            labels = set()
            for role in (InteractionRole.DONOR, InteractionRole.ACCEPTOR):
                site = INTERACTION_SITES[InteractionType.HYDROGEN_BOND][role].get(group)
                if site is not None:
                    labels.update(site.atom_map_numbers)
            self.assertEqual(selector.queries[0].smarts, MAPPED_FUNCTIONAL_GROUP_SMARTS[group])
            self.assertEqual(
                selector.queries[0].atom_map_numbers, tuple(sorted(labels)) if labels else None
            )

        for interaction, roles in INTERACTION_SITES.items():
            for role, sites in roles.items():
                alias = f"interaction {interaction.value.replace('_', ' ')} {role.value.replace('_', ' ')}"
                matches = [
                    selector for name, selector in SELECTORS.items()
                    if name == alias or alias in selector.aliases
                ]
                self.assertEqual(len(matches), 1, alias)
                self.assertEqual(
                    {(query.smarts, query.atom_map_numbers) for query in matches[0].queries},
                    {(site.smarts, site.atom_map_numbers) for site in sites.values()},
                )

    def test_composite_queries_and_aliases(self) -> None:
        for name in ("alcohol", "amine", "amide"):
            selector = SELECTORS[name]
            self.assertEqual(len(selector.queries), 3)
            self.assertEqual(json.loads(selector.smarts), [query.smarts for query in selector.queries])
        self.assertEqual(
            [query.atom_map_numbers for query in SELECTORS["amide"].queries],
            [(2, 3), (2, 3), (2, 3)],
        )
        self.assertIn("sulfide", SELECTORS["thioether"].aliases)
        self.assertIn("hbond donor", SELECTORS["hydrogen bond donor"].aliases)


class InteractionRoleTests(unittest.TestCase):
    def test_hydrogen_bond_donor_and_acceptor_roles(self) -> None:
        donors = get_groups_for_role(InteractionType.HYDROGEN_BOND, InteractionRole.DONOR)
        acceptors = get_groups_for_role("hydrogen_bond", "acceptor")

        self.assertIn(FunctionalGroup.PHENOL, donors)
        self.assertIn(FunctionalGroup.CARBOXYLATE, acceptors)
        self.assertNotIn(FunctionalGroup.TERTIARY_AMIDE, donors)
        self.assertIn(FunctionalGroup.TERTIARY_AMIDE, acceptors)

    def test_specific_interaction_roles(self) -> None:
        halogen_donors = get_groups_for_role(
            InteractionType.HALOGEN_BOND,
            InteractionRole.SIGMA_HOLE_DONOR,
        )
        coordination_ligands = get_groups_for_role(
            InteractionType.METAL_COORDINATION,
            InteractionRole.COORDINATION_LIGAND,
        )

        self.assertEqual(
            halogen_donors,
            frozenset({FunctionalGroup.ALKYL_HALIDE, FunctionalGroup.ARYL_HALIDE}),
        )
        self.assertIn(FunctionalGroup.CARBOXYLATE, coordination_ligands)

    def test_non_pair_interactions_have_no_donor_acceptor_mapping(self) -> None:
        self.assertEqual(
            get_groups_for_role(InteractionType.HYDROPHOBIC_EFFECT, InteractionRole.DONOR),
            frozenset(),
        )
        self.assertIs(BOND_TYPE_ROLES, INTERACTION_ROLES)

    def test_unknown_interaction_or_role_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            get_groups_for_role("unknown", InteractionRole.DONOR)
        with self.assertRaises(ValueError):
            get_groups_for_role(InteractionType.IONIC, "unknown")


class InteractionSiteTests(unittest.TestCase):
    def test_sites_cover_existing_mappings_and_are_immutable(self) -> None:
        self.assertEqual(set(MAPPED_FUNCTIONAL_GROUP_SMARTS), set(FunctionalGroup))
        self.assertEqual(set(INTERACTION_SITES), set(INTERACTION_ROLES))
        for interaction, roles in INTERACTION_ROLES.items():
            self.assertEqual(set(INTERACTION_SITES[interaction]), set(roles))
            for role, groups in roles.items():
                self.assertEqual(set(INTERACTION_SITES[interaction][role]), set(groups))
        with self.assertRaises(TypeError):
            INTERACTION_SITES[InteractionType.HYDROGEN_BOND][InteractionRole.DONOR][
                FunctionalGroup.PHENOL
            ] = None

    def test_invalid_or_unmapped_sites_raise(self) -> None:
        with self.assertRaises(ValueError):
            get_interaction_site("unknown", "donor", "phenol")
        with self.assertRaises(ValueError):
            get_interaction_site("hydrogen_bond", "unknown", "phenol")
        with self.assertRaises(ValueError):
            get_interaction_site("hydrogen_bond", "donor", "unknown")
        with self.assertRaises(KeyError):
            get_interaction_site("hydrogen_bond", "donor", "tertiary_amide")
        with self.assertRaises(KeyError):
            get_interaction_site("hydrogen_bond", "cation", "ammonium")

    def test_all_mapped_patterns_preserve_original_queries(self) -> None:
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("RDKit is optional")

        for group, original in FUNCTIONAL_GROUP_SMARTS.items():
            with self.subTest(group=group):
                mapped = Chem.MolFromSmarts(MAPPED_FUNCTIONAL_GROUP_SMARTS[group])
                query = Chem.MolFromSmarts(original)
                self.assertIsNotNone(mapped)
                self.assertIsNotNone(query)
                labels = [atom.GetAtomMapNum() for atom in mapped.GetAtoms()]
                self.assertTrue(all(labels))
                self.assertEqual(len(labels), len(set(labels)))
                for atom in mapped.GetAtoms():
                    atom.SetAtomMapNum(0)
                self.assertEqual(Chem.MolToSmarts(mapped), Chem.MolToSmarts(query))

    def test_all_site_labels_exist_in_their_queries(self) -> None:
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("RDKit is optional")

        for interaction, roles in INTERACTION_SITES.items():
            for role, groups in roles.items():
                for group, site in groups.items():
                    with self.subTest(interaction=interaction, role=role, group=group):
                        query = Chem.MolFromSmarts(site.smarts)
                        self.assertIsNotNone(query)
                        labels = [atom.GetAtomMapNum() for atom in query.GetAtoms()]
                        self.assertEqual(len(labels), len(set(labels)))
                        self.assertTrue(site.atom_map_numbers)
                        self.assertTrue(set(site.atom_map_numbers).issubset(labels))

    def test_role_specific_sites_and_conditional_ammonium_donor(self) -> None:
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("RDKit is optional")

        def find_sites(smiles, interaction, role, group):
            site = get_interaction_site(interaction, role, group)
            query = Chem.MolFromSmarts(site.smarts)
            molecule = Chem.MolFromSmiles(smiles)
            positions = [
                atom.GetIdx() for atom in query.GetAtoms()
                if atom.GetAtomMapNum() in site.atom_map_numbers
            ]
            return [
                tuple(match[position] for position in positions)
                for match in molecule.GetSubstructMatches(query)
            ]

        cases = (
            ("CC(=O)O", "hydrogen_bond", "donor", "carboxylic_acid", [(3,)]),
            ("CC(=O)O", "hydrogen_bond", "acceptor", "carboxylic_acid", [(2,)]),
            ("CC(=O)N", "hydrogen_bond", "donor", "primary_amide", [(3,)]),
            ("CC(=O)N", "hydrogen_bond", "acceptor", "primary_amide", [(2,)]),
            ("CC(=O)[O-]", "hydrogen_bond", "acceptor", "carboxylate", [(2, 3)]),
            ("CC([O-])=O", "hydrogen_bond", "acceptor", "carboxylate", [(3, 2)]),
            ("CSC", "chalcogen_bond", "sigma_hole_donor", "thioether", [(1,)]),
            ("CCl", "halogen_bond", "sigma_hole_donor", "alkyl_halide", [(1,)]),
            ("CN", "hydrogen_bond", "donor", "primary_amine", [(1,)]),
            ("C[NH3+]", "hydrogen_bond", "donor", "ammonium", [(1,)]),
            ("C[N+](C)(C)C", "hydrogen_bond", "donor", "ammonium", []),
            ("C[N+](C)(C)C", "ionic", "cation", "ammonium", [(1,)]),
        )
        for smiles, interaction, role, group, expected in cases:
            with self.subTest(smiles=smiles, interaction=interaction, role=role, group=group):
                self.assertEqual(find_sites(smiles, interaction, role, group), expected)


if __name__ == "__main__":
    unittest.main()
