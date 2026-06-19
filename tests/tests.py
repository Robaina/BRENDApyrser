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

    def test_substrates_and_products_preserve_ions(self):
        # PR #13 issue #1: splitting reactions on a bare "+" dropped ion charges
        # (e.g. "H+" -> "H", "NAD+" -> "NAD"). Split on " + " keeps them intact.
        sp = self.rxn.substratesAndProducts
        self.assertEqual(len(sp), 1)
        self.assertIn("H+", sp[0]["substrates"])
        self.assertIn("NAD+", sp[0]["products"])
        self.assertNotIn("H", sp[0]["substrates"])
        self.assertNotIn("NAD", sp[0]["products"])

    def test_filter_by_compound_matches_ions(self):
        # Ionic compounds are now filterable because their "+" is preserved.
        self.assertEqual(
            [r.ec_number for r in self.db.reactions.filter_by_compound("NAD+")],
            ["1.1.1.304"],
        )

    def test_substrates_products_strip_annotations(self):
        # PR #13 issue #2: |...| pipe comments and <refs> must not leak into the
        # compound string; the pipe note is folded into `comment`, not discarded.
        for sp in self.rxn.substrates_products:
            self.assertNotIn("|", sp["value"])
            self.assertNotIn("<", sp["value"])
            self.assertNotIn("{", sp["value"])
        rec = next(
            sp
            for sp in self.rxn.substrates_products
            if sp["value"] == "2,3-pentanedione + NADH + H+ = ? + NAD+"
        )
        self.assertIn("69% of the activity with diacetyl", rec["comment"])

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


class TestCompoundNameParsing(unittest.TestCase):
    """Regression test for issue #12: compound names containing parentheses must
    not be truncated at the first '('. The legacy text parser extracted the
    first ``(...)`` group as a comment, cutting the EC 5.6.2.2 inhibitor below
    down to ``6-methoxy-4-``."""

    # Real EC 5.6.2.2 inhibitor; the nested ()/[] are the trap.
    COMPOUND = (
        "6-methoxy-4-(2-[4-[([1,3]oxathiolo[5,4-c]pyridin-6-ylmethyl)amino]"
        "piperidin-1-yl]ethyl)quinoline-3-carbonitrile"
    )

    def setUp(self):
        self.rxn = Reaction(
            {
                "id": "5.6.2.2",
                "protein": {
                    "1": {
                        "id": "1",
                        "organism": "Homo sapiens",
                        "references": ["1"],
                        "comment": "",
                    }
                },
                "reference": {
                    "1": {
                        "id": "1",
                        "title": "t",
                        "authors": ["a"],
                        "journal": "j",
                        "year": 2020,
                    }
                },
                "inhibitor": [
                    {
                        "value": self.COMPOUND,
                        "proteins": ["1"],
                        "references": ["1"],
                        "comment": "",
                    }
                ],
                "ki_value": [
                    {
                        "value": f"0.5 {{{self.COMPOUND}}}",
                        "proteins": ["1"],
                        "references": ["1"],
                        "comment": "",
                    }
                ],
            }
        )

    def test_inhibitor_name_not_truncated(self):
        names = list(self.rxn.inhibitors.keys())
        self.assertIn(self.COMPOUND, names)
        # the legacy bug truncated this to "6-methoxy-4-"
        self.assertNotIn("6-methoxy-4-", names)

    def test_kinetic_substrate_name_not_truncated(self):
        # the parenthesised compound is the {substrate} key of the Ki value
        self.assertIn(self.COMPOUND, self.rxn.KIvalues.keys())


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


class TestBRENDAExtras(unittest.TestCase):
    """Database-level accessors not covered by TestBRENDA."""

    @classmethod
    def setUpClass(cls):
        cls.db = BRENDA(FIXTURE)

    def test_fields_and_units_exposed(self):
        self.assertTrue(self.db.fields)
        self.assertTrue(self.db.units)

    def test_copyright(self):
        self.assertIn("Braunschweig", self.db.copyright)

    def test_repr_html(self):
        html = self.db._repr_html_()
        self.assertIn("2026.1", html)
        self.assertIn("Number of Enzymes", html)


class TestReactionList(unittest.TestCase):
    """ReactionList lookups, filters and container semantics."""

    @classmethod
    def setUpClass(cls):
        cls.db = BRENDA(FIXTURE)

    def test_get_by_name_case_insensitive(self):
        rxn = self.db.reactions.get_by_name("DIACETYL REDUCTASE [(S)-ACETOIN FORMING]")
        self.assertEqual(rxn.ec_number, "1.1.1.304")

    def test_get_by_name_missing_raises(self):
        with self.assertRaises(ValueError):
            self.db.reactions.get_by_name("no such enzyme")

    def test_filter_by_substrate(self):
        hits = self.db.reactions.filter_by_substrate("NADH")
        self.assertEqual([r.ec_number for r in hits], ["1.1.1.304"])

    def test_filter_by_product(self):
        hits = self.db.reactions.filter_by_product("NAD+")
        self.assertEqual([r.ec_number for r in hits], ["1.1.1.304"])

    def test_filter_no_match_returns_empty(self):
        self.assertEqual(self.db.reactions.filter_by_substrate("unobtainium"), [])

    def test_slicing_returns_reactionlist(self):
        from brendapyrser.parser import ReactionList

        self.assertIsInstance(self.db.reactions[:1], ReactionList)
        self.assertEqual(len(self.db.reactions[:1]), 1)

    def test_filter_by_organism_returns_reactionlist(self):
        from brendapyrser.parser import ReactionList

        self.assertIsInstance(
            self.db.reactions.filter_by_organism("Staphylococcus"), ReactionList
        )


class TestReactionProperties(unittest.TestCase):
    """Reaction properties beyond the core kinetics already covered."""

    @classmethod
    def setUpClass(cls):
        cls.db = BRENDA(FIXTURE)
        cls.rxn = cls.db.reactions.get_by_id("1.1.1.304")

    def test_metals_excludes_more(self):
        self.assertIn("Mn2+", self.rxn.metals)
        self.assertNotIn("more", self.rxn.metals)

    def test_activators(self):
        self.assertIn("DMSO", self.rxn.activators)

    def test_inhibitors(self):
        self.assertIn("diacetyl", self.rxn.inhibitors)
        self.assertNotIn("more", self.rxn.inhibitors)

    def test_specific_activities(self):
        first = self.rxn.specificActivities[0]
        self.assertEqual(first["value"], 72.6)
        self.assertIn("species", first)
        self.assertIn("refs", first)

    def test_mechanism(self):
        self.assertEqual(
            self.rxn.mechanism, ["(S)-acetoin + NAD+ = diacetyl + NADH + H+"]
        )

    def test_reaction_type_absent_is_empty(self):
        # The fixture entry carries no `reaction_type`; must degrade to [].
        self.assertEqual(self.rxn.reaction_type, [])

    def test_field_backed_properties(self):
        for records in (
            self.rxn.molecular_weight,
            self.rxn.subunits,
            self.rxn.pi_value,
            self.rxn.source_tissue,
        ):
            self.assertIsInstance(records, list)
            self.assertTrue(records)
            self.assertTrue(all("value" in r and "comment" in r for r in records))

    def test_molecular_weight_enriched_with_citation(self):
        first = self.rxn.molecular_weight[0]
        self.assertEqual(first["value"], "68000")
        self.assertIn("Staphylococcus aureus", first["organisms"])

    def test_references_embed_pubmed(self):
        self.assertIn("Pubmed:3041963", self.rxn.references["4"])


class TestConditions(unittest.TestCase):
    """pH / temperature EnzymeConditionDict behaviour."""

    @classmethod
    def setUpClass(cls):
        cls.rxn = BRENDA(FIXTURE).reactions.get_by_id("1.1.1.304")

    def test_ph_optimum(self):
        self.assertEqual(self.rxn.PH["optimum"][0]["value"], 6.0)

    def test_ph_range_is_pair(self):
        self.assertEqual(self.rxn.PH["range"][0]["value"], [5.0, 8.0])

    def test_filter_by_condition_valid(self):
        only = self.rxn.temperature.filter_by_condition("optimum")
        self.assertEqual(list(only.keys()), ["optimum"])

    def test_get_values(self):
        self.assertTrue(self.rxn.PH.get_values())


class TestPropertyDicts(unittest.TestCase):
    """EnzymePropertyDict filtering helpers."""

    @classmethod
    def setUpClass(cls):
        cls.rxn = BRENDA(FIXTURE).reactions.get_by_id("1.1.1.304")

    def test_filter_by_compound_hit(self):
        self.assertEqual(
            list(self.rxn.KMvalues.filter_by_compound("NADH").keys()), ["NADH"]
        )

    def test_filter_by_compound_miss_returns_empty_list(self):
        self.assertEqual(self.rxn.KMvalues.filter_by_compound("zzz"), {"zzz": []})

    def test_filter_by_organism(self):
        filtered = self.rxn.KMvalues.filter_by_organism("Staphylococcus aureus")
        self.assertTrue(
            all(
                any("Staphylococcus aureus" in s for s in v["species"])
                for recs in filtered.values()
                for v in recs
            )
        )


class TestRangeEvaluation(unittest.TestCase):
    """Numeric/range coercion via constructed entries (deterministic)."""

    def _reaction(self, **fields):
        base = {"id": "9.9.9.9", "protein": {}, "reference": {}}
        base.update(fields)
        return Reaction(base)

    def test_km_range_is_averaged(self):
        rxn = self._reaction(
            km_value=[{"value": "1.0-3.0 {X}", "proteins": [], "references": []}]
        )
        self.assertEqual(rxn.KMvalues.get_values(), [2.0])

    def test_ph_optimum_range_is_averaged(self):
        rxn = self._reaction(
            ph_optimum=[{"value": "6-8", "proteins": [], "references": []}]
        )
        self.assertEqual(rxn.PH["optimum"][0]["value"], 7.0)

    def test_ph_range_kept_as_pair(self):
        rxn = self._reaction(
            ph_range=[{"value": "5-8", "proteins": [], "references": []}]
        )
        self.assertEqual(rxn.PH["range"][0]["value"], [5.0, 8.0])

    def test_unparseable_scalar_is_sentinel(self):
        rxn = self._reaction(
            temperature_optimum=[{"value": "oops", "proteins": [], "references": []}]
        )
        self.assertEqual(rxn.temperature["optimum"][0]["value"], -999)

    def test_unparseable_range_is_sentinel_pair(self):
        rxn = self._reaction(
            ph_range=[{"value": "low-high", "proteins": [], "references": []}]
        )
        self.assertEqual(rxn.PH["range"][0]["value"], [-999, -999])


class TestDisplay(unittest.TestCase):
    """Jupyter/summary rendering helpers."""

    @classmethod
    def setUpClass(cls):
        cls.rxn = BRENDA(FIXTURE).reactions.get_by_id("1.1.1.304")

    def test_summary_dataframe(self):
        summary = self.rxn.summary
        self.assertEqual(summary.loc["EC number"].iloc[0], "1.1.1.304")

    def test_repr_html(self):
        html = self.rxn._repr_html_()
        self.assertIn("1.1.1.304", html)
        self.assertIn("Enzyme identifier", html)


if __name__ == "__main__":
    unittest.main()
