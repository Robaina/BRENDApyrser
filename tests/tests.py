#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for brendapyrser (BRENDA 2026 JSON format).
"""

import gzip
import json
import os
import shutil
import tarfile
import tempfile
import unittest

from brendapyrser import BRENDA, Reaction

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "brenda_sample.json")


class TestBRENDA(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db = BRENDA(FIXTURE)

    def test_metadata(self):
        self.assertEqual(self.db.release, "2026.1")
        self.assertEqual(self.db.schema_version, "1")

    def test_number_of_reactions(self):
        # "spontaneous" is excluded; the fixture holds two EC entries.
        self.assertEqual(len(self.db.reactions), 2)

    def test_get_by_id(self):
        rxn = self.db.reactions.get_by_id("1.1.1.304")
        self.assertEqual(rxn.ec_number, "1.1.1.304")

    def test_get_by_id_missing_raises(self):
        with self.assertRaises(ValueError):
            self.db.reactions.get_by_id("9.9.9.9")

    def test_filter_by_organism(self):
        hits = self.db.reactions.filter_by_organism("Staphylococcus")
        self.assertEqual([r.ec_number for r in hits], ["1.1.1.304"])

    def test_getKMcompounds(self):
        self.assertIn("NADH", self.db.getKMcompounds())

    def test_getOrganisms(self):
        organisms = self.db.getOrganisms()
        self.assertIn("Staphylococcus aureus", organisms)
        self.assertNotIn("", organisms)


class TestReaction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db = BRENDA(FIXTURE)
        cls.rxn = cls.db.reactions.get_by_id("1.1.1.304")

    def test_ec_number(self):
        self.assertEqual(self.rxn.ec_number, "1.1.1.304")

    def test_name_verbatim(self):
        # Name is returned verbatim from BRENDA (chemical casing preserved),
        # unlike the legacy parser which lower-cased it via str.capitalize().
        self.assertEqual(self.rxn.name, "diacetyl reductase [(S)-acetoin forming]")
        self.assertIn("(S)", self.rxn.name)

    def test_sysname(self):
        self.assertEqual(self.rxn.systematic_name, "(S)-acetoin:NAD+ oxidoreductase")

    def test_reaction_str(self):
        self.assertEqual(
            self.rxn.reaction_str, "(S)-acetoin + NAD+ = diacetyl + NADH + H+"
        )

    def test_KMvalues(self):
        self.assertEqual(
            self.rxn.KMvalues.get_values()[:4], [0.045, 0.095, 0.025, 0.11]
        )

    def test_KKMvalues(self):
        self.assertEqual(
            self.rxn.KKMvalues.get_values()[:4], [210.0, 5072.0, 8.519, 36.4]
        )

    def test_Kcatvalues(self):
        self.assertEqual(
            self.rxn.Kcatvalues.get_values()[:4], [110.0, 1274.0, 163.0, 1222.0]
        )

    def test_temperature_optimum(self):
        self.assertEqual(self.rxn.temperature["optimum"][0]["value"], 50.0)

    def test_temperature_condition_keys(self):
        with self.assertRaises(KeyError):
            self.rxn.temperature.filter_by_condition("not_a_condition")

    def test_organisms(self):
        self.assertEqual(len(self.rxn.organisms), 10)
        self.assertIn("Staphylococcus aureus", self.rxn.organisms)

    def test_proteins_enriched(self):
        protein = self.rxn.proteins["7"]
        self.assertEqual(protein["organism"], "Klebsiella pneumoniae")
        self.assertEqual(protein["name"], "Klebsiella pneumoniae")
        self.assertEqual(protein["source"], "uniprot")
        self.assertEqual(protein["accessions"], ["Q48436"])
        self.assertEqual(protein["proteinID"], "Q48436")

    def test_references_and_bibliography(self):
        self.assertIn("Vidal", self.rxn.references["4"])
        self.assertEqual(self.rxn.bibliography["4"]["pmid"], 3041963)
        self.assertEqual(self.rxn.bibliography["4"]["year"], 1988)

    def test_substrates_products_reversibility(self):
        first = self.rxn.substrates_products[0]
        self.assertEqual(first["reversibility"], "irreversible")
        self.assertNotIn("{ir}", first["value"])

    def test_synonyms(self):
        self.assertIn("ADS1", self.rxn.synonyms)

    def test_cofactors(self):
        self.assertIn("NADH", self.rxn.cofactors.keys())

    def test_field_accessor(self):
        cloned = self.rxn.field("cloned")
        self.assertTrue(all("value" in rec and "organisms" in rec for rec in cloned))

    def test_direct_construction_from_dict(self):
        with open(FIXTURE, encoding="utf-8") as fh:
            entry = json.load(fh)["data"]["1.1.1.304"]
        rxn = Reaction(entry)
        self.assertEqual(rxn.ec_number, "1.1.1.304")
        self.assertEqual(rxn.KMvalues.get_values()[:2], [0.045, 0.095])


class TestEdgeCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db = BRENDA(FIXTURE)

    def test_missing_recommended_name(self):
        # EC 6.6.99.99 is a transferred entry with no recommended_name.
        rxn = self.db.reactions.get_by_id("6.6.99.99")
        self.assertEqual(rxn.name, "")
        self.assertEqual(rxn.systematic_name, "thermophilic enzyme")

    def test_empty_kinetics_and_organisms(self):
        rxn = self.db.reactions.get_by_id("6.6.99.99")
        self.assertEqual(dict(rxn.KMvalues), {})
        self.assertEqual(rxn.organisms, [])


class TestLoaders(unittest.TestCase):
    """The loader must accept plain JSON, gzip, and tar.gz inputs."""

    def _assert_loads(self, path):
        db = BRENDA(path)
        self.assertEqual(db.release, "2026.1")
        self.assertEqual(len(db.reactions), 2)

    def test_plain_json(self):
        self._assert_loads(FIXTURE)

    def test_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            gz_path = os.path.join(tmp, "brenda_sample.json.gz")
            with open(FIXTURE, "rb") as src, gzip.open(gz_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
            self._assert_loads(gz_path)

    def test_tar_gz(self):
        with tempfile.TemporaryDirectory() as tmp:
            tar_path = os.path.join(tmp, "brenda_sample.json.tar.gz")
            with tarfile.open(tar_path, "w:gz") as tar:
                tar.add(FIXTURE, arcname="brenda_sample.json")
            self._assert_loads(tar_path)


if __name__ == "__main__":
    unittest.main()
