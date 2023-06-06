import unittest

from taxonomy import Taxonomy, TaxonomyError

JSON_DATA = """
{
  "multigraph": false,
  "directed": true,
  "graph": [],
  "nodes": [
    {
      "name": "root",
      "rank": "no rank",
      "id": 1
    },
    {
      "name": "genus 2",
      "rank": "genus",
      "id": 9
    },
    {
      "name": "superkingdom 1",
      "rank": "superkingdom",
      "id": 2
    },
    {
      "name": "species 2.1",
      "rank": "species",
      "id": 11
    },
    {
      "name": "genus 1",
      "rank": "genus",
      "id": 8
    },
    {
      "name": "class 1",
      "rank": "class",
      "id": 5
    },
    {
      "name": "kingdom 1",
      "rank": "kingdom",
      "id": 3
    },
    {
      "name": "phylum 1",
      "rank": "phylum",
      "id": 4
    },
    {
      "name": "order 1",
      "rank": "order",
      "id": 6
    },
    {
      "name": "family 1",
      "rank": "family",
      "id": 7
    },
    {
      "name": "species 1.1",
      "rank": "species",
      "id": 10
    },
    {
      "name": "species 1.1",
      "rank": "species",
      "id": 12
    }
  ],
  "links": [
    {
      "source": 3,
      "target": 1
    },
    {
      "source": 10,
      "target": 4
    },
    {
      "source": 1,
      "target": 9
    },
    {
      "source": 4,
      "target": 9
    },
    {
      "source": 9,
      "target": 8
    },
    {
      "source": 8,
      "target": 5
    },
    {
      "source": 5,
      "target": 7
    },
    {
      "source": 7,
      "target": 6
    },
    {
      "source": 6,
      "target": 2
    },
    {
      "source": 2,
      "target": 0
    },
    {
      "source": 11,
      "target": 4
    }
  ]
}        """


class JsonTestCase(unittest.TestCase):
    def _create_tax(self):
        return Taxonomy.from_json(JSON_DATA)

    def setUp(self) -> None:
        self.tax = self._create_tax()

    def test_internal_index(self):
        self.assertEqual(
            [
                self.tax.internal_index(x)
                for x in ["1", "9", "2", "11", "8", "5", "3", "4", "6", "7", "10"]
            ],
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        )

    def test_find_all_by_name(self):
        self.assertEqual(
            sorted([n.id for n in self.tax.find_all_by_name("species 1.1")]), ["10", "12"]
        )

    def test_edit_node_parent_updates_children(self):
        assert self.tax["5"].parent == "4"

        # original parent = 4
        self.tax.edit_node("5", parent_id="1")

        node = self.tax["5"]
        assert node.parent == "1"

        assert "5" not in {n.id for n in self.tax.children("4")}
        assert "5" in {n.id for n in self.tax.children("1")}


class NewickTestCase(unittest.TestCase):
    def _create_tax(self):
        # https://en.wikipedia.org/wiki/Newick_format#Examples
        return Taxonomy.from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")

    def setUp(self) -> None:
        self.tax = self._create_tax()

    def test_root(self):
        root = self.tax.root
        self.assertEqual(root.id, "F")
        self.assertIsNone(root.parent)

    def test_find_node_by_id(self):
        node = self.tax.node("A")
        self.assertEqual(node, self.tax.node("A"))
        self.assertEqual(node.id, "A")
        self.assertEqual(node.parent, "F")

        node = self.tax.node("D")
        self.assertEqual(node.id, "D")
        self.assertEqual(node.parent, "E")

        node = self.tax.node("unknown")
        self.assertIsNone(node)

    def test_index(self):
        node = self.tax["A"]
        self.assertEqual(node.id, "A")
        self.assertEqual(node.parent, "F")

        with self.assertRaises(TaxonomyError):
            _ = self.tax["unknown"]

    def test_find_all_by_name(self):
        # They are not named so we can't find anything by name
        nodes = self.tax.find_all_by_name("A")
        self.assertEqual(nodes, [])

    def test_parent(self):
        parent = self.tax.parent("D")
        self.assertEqual(parent.id, "E")

    def test_parent_with_distance(self):
        parent, distance = self.tax.parent_with_distance("D")
        self.assertEqual(parent.id, "E")
        # Float precision issue
        # 0.4 becomes 0.4000000059604645
        self.assertAlmostEqual(distance, 0.4)

    def test_children(self):
        children = self.tax.children("E")
        self.assertEqual(len(children), 2)
        self.assertEqual(children[0].id, "C")
        self.assertEqual(children[1].id, "D")

    def test_lineage(self):
        lineage = self.tax.lineage("D")
        self.assertEqual(len(lineage), 3)
        self.assertEqual(lineage[0].id, "D")
        self.assertEqual(lineage[1].id, "E")
        self.assertEqual(lineage[2].id, "F")

    def test_parents(self):
        lineage = self.tax.parents("D")
        self.assertEqual(len(lineage), 2)
        self.assertEqual(lineage[0].id, "E")
        self.assertEqual(lineage[1].id, "F")

    def test_lca(self):
        lca = self.tax.lca("A", "D")
        self.assertEqual(lca.id, "F")

    def test_prune(self):
        new_tax = self.tax.prune(remove=["E"])
        self.assertIsNone(new_tax.node("D"))
        self.assertIsNone(new_tax.node("E"))
        # We removed a leaf
        self.assertEqual(len(new_tax), 3)

        new_tax = self.tax.prune(keep=["E", "D"])
        self.assertEqual(len(new_tax), 3)
        self.assertIsNotNone(new_tax.node("F"))

    def test_remove(self):
        tax = self._create_tax()
        tax.remove_node("E")
        self.assertIsNotNone(tax.node("D"))
        self.assertIsNone(tax.node("E"))
        self.assertEqual(len(tax), 5)

    def test_add(self):
        tax = self._create_tax()
        tax.add_node("D", "G", "something", "species")
        node = tax["G"]
        self.assertEqual(node.parent, "D")

        tax.add_node("G", "H", "something else", "species")
        node = tax["H"]
        self.assertEqual(node.parent, "G")

    def test_edit_node(self):
        tax = self._create_tax()
        tax.edit_node("D", parent_distance=3)
        node, distance = tax.parent_with_distance("D")
        self.assertEqual(distance, 3)

    def test_can_clone(self):
        tax = self._create_tax()
        tax2 = tax.clone()
        tax.remove_node("E")

        self.assertIsNotNone(tax.node("D"))
        self.assertIsNone(tax.node("E"))
        self.assertEqual(len(tax), 5)

        self.assertIsNotNone(tax2.node("D"))
        self.assertIsNotNone(tax2.node("E"))
        self.assertEqual(len(tax2), 6)


class NCBITestCase(unittest.TestCase):
    def _create_tax(self):
        return Taxonomy.from_ncbi("tests/data/")

    def setUp(self) -> None:
        self.tax = self._create_tax()

    def test_root(self):
        root = self.tax.root
        self.assertEqual(root.id, "1")
        self.assertIsNone(root.parent)

    def test_find_node_by_id(self):
        node = self.tax.node("1236")
        self.assertEqual(node.id, "1236")
        self.assertEqual(node.name, "Gammaproteobacteria")
        self.assertEqual(node.parent, "1224")

        node = self.tax.node("unknown")
        self.assertIsNone(node)

    def test_index(self):
        node = self.tax["1236"]
        self.assertEqual(node.id, "1236")
        self.assertEqual(node.name, "Gammaproteobacteria")
        self.assertEqual(node.parent, "1224")

        with self.assertRaises(TaxonomyError):
            _ = self.tax["unknown"]

    def test_find_all_by_name(self):
        nodes = self.tax.find_all_by_name("Escherichia coli")
        self.assertEqual([n.id for n in nodes], ["562"])
        self.assertEqual([n.name for n in nodes], ["Escherichia coli"])
        self.assertEqual([n.parent for n in nodes], ["561"])

    def test_parent(self):
        parent = self.tax.parent("562")
        self.assertEqual(parent.id, "561")

    def test_parent_with_distance(self):
        parent, distance = self.tax.parent_with_distance("562")
        self.assertEqual(parent.id, "561")
        # Float precision issue
        # 0.4 becomes 0.4000000059604645
        self.assertAlmostEqual(distance, 1.0)

    def test_children(self):
        children = self.tax.children("561")
        self.assertEqual(len(children), 1)
        self.assertEqual(children[0].id, "562")

    def test_lineage(self):
        lineage = self.tax.lineage("562")
        self.assertEqual(len(lineage), 9)
        self.assertEqual(lineage[0].id, "562")
        self.assertEqual(lineage[1].id, "561")
        self.assertEqual(lineage[-1].id, "1")

    def test_parents(self):
        lineage = self.tax.parents("562")
        self.assertEqual(len(lineage), 8)
        self.assertEqual(lineage[0].id, "561")
        self.assertEqual(lineage[-1].id, "1")

    def test_lca(self):
        lca = self.tax.lca("562", "91347")
        self.assertEqual(lca.id, "91347")

    def test_prune(self):
        new_tax = self.tax.prune(remove=["561"])
        self.assertIsNone(new_tax.node("561"))
        self.assertIsNone(new_tax.node("562"))
        self.assertEqual(len(new_tax), 7)

        new_tax = self.tax.prune(keep=["561"])
        self.assertEqual(len(new_tax), 8)
        self.assertIsNotNone(new_tax.node("561"))

    @unittest.skip("tax.remove doesn't work on truncated taxonomies?")
    def test_remove(self):
        tax = self._create_tax()
        tax.remove_node("561")
        self.assertIsNotNone(tax.node("562"))
        self.assertIsNone(tax.node("561"))
        self.assertEqual(len(tax), 8)

    def test_add(self):
        tax = self._create_tax()
        tax.add_node("561", "563", "Listeria", "species")
        node = tax["563"]
        self.assertEqual(node.parent, "561")
        self.assertEqual(node.name, "Listeria")
        self.assertEqual(node.rank, "species")
        tax.add_node("563", "100000001", "Pizzeria", "genus")
        node = tax["100000001"]
        self.assertEqual(node.parent, "563")
        self.assertEqual(node.name, "Pizzeria")
        self.assertEqual(node.rank, "genus")

    def test_cannot_add_duplicate_tax_id(self):
        tax = self._create_tax()
        tax.add_node("561", "563", "Listeria", "species")

        with self.assertRaises(TaxonomyError) as context:
            tax.add_node("561", "563", "Listeria", "species")
        self.assertTrue("563" in str(context.exception))

    def test_edit_node(self):
        tax = self._create_tax()
        tax.edit_node("562", parent_distance=3)
        node, distance = tax.parent_with_distance("562")
        self.assertEqual(distance, 3)

    def test_edit_node_parent(self):
        tax = self._create_tax()
        self.assertEqual(tax["562"].parent, "561")
        tax.edit_node("562", parent_id="1")
        self.assertEqual(tax["562"].parent, "1")

    def test_repr(self):
        tax = self._create_tax()
        self.assertEqual(
            tax["562"].__repr__(),
            '<TaxonomyNode (id="562" rank="species" name="Escherichia coli")>',
        )


class GtdbTestCase(unittest.TestCase):
    def setUp(self):
        with open("tests/data/gtdb_sample.tsv") as file:
            self.tax = Taxonomy.from_gtdb(file.read())

    def test_root(self):
        root = self.tax.root
        self.assertEqual(root.id, "d__Bacteria")
        self.assertEqual(root.rank, "domain")
        self.assertIsNone(root.parent)

    def test_lineage(self):
        self.assertEqual([n.id for n in self.tax.lineage("d__Bacteria")], ["d__Bacteria"])

        self.assertEqual(
            [n.id for n in self.tax.lineage("c__Bacilli")],
            ["c__Bacilli", "p__Firmicutes", "d__Bacteria"],
        )

        self.assertEqual(
            [n.id for n in self.tax.lineage("s__Escherichia coli")],
            [
                "s__Escherichia coli",
                "g__Escherichia",
                "f__Enterobacteriaceae",
                "o__Enterobacterales",
                "c__Gammaproteobacteria",
                "p__Proteobacteria",
                "d__Bacteria",
            ],
        )

    def test_invalid_format(self):
        with open("tests/data/gtdb_invalid.tsv") as file:
            with self.assertRaises(TaxonomyError):
                Taxonomy.from_gtdb(file.read())


if __name__ == "__main__":
    unittest.main()
