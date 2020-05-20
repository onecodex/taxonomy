import unittest

from taxonomy import Taxonomy, TaxonomyError


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
        node = self.tax.find_by_id("A")
        self.assertEqual(node.id, "A")
        self.assertEqual(node.parent, "F")

        node = self.tax.find_by_id("D")
        self.assertEqual(node.id, "D")
        self.assertEqual(node.parent, "E")

        node = self.tax.find_by_id("unknown")
        self.assertIsNone(node)

    def test_index(self):
        node = self.tax["A"]
        self.assertEqual(node.id, "A")
        self.assertEqual(node.parent, "F")

        with self.assertRaises(TaxonomyError):
            _ = self.tax["unknown"]

    def test_find_by_name(self):
        # They are not named so we can't find anything by name
        node = self.tax.find_by_name("A")
        self.assertIsNone(node)

    def test_parent(self):
        parent, distance = self.tax.parent("D")
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

    def test_lca(self):
        lca = self.tax.lca("A", "D")
        self.assertEqual(lca.id, "F")

    def test_prune(self):
        new_tax = self.tax.prune(remove=["E"])
        self.assertIsNone(new_tax.find_by_id("D"))
        self.assertIsNone(new_tax.find_by_id("E"))
        # We removed a leaf
        self.assertEqual(len(new_tax), 3)

        new_tax = self.tax.prune(keep=["E", "D"])
        self.assertEqual(len(new_tax), 3)
        self.assertIsNotNone(new_tax.find_by_id("F"))

    def test_remove(self):
        tax = self._create_tax()
        tax.remove("E")
        self.assertIsNotNone(tax.find_by_id("D"))
        self.assertIsNone(tax.find_by_id("E"))
        self.assertEqual(len(tax), 5)

    def test_add(self):
        tax = self._create_tax()
        tax.add("D", "G")
        node = tax["G"]
        self.assertEqual(node.parent, "D")

    def test_edit_node(self):
        tax = self._create_tax()
        tax.edit_node("D", parent_distance=3)
        node, distance = tax.parent("D")
        self.assertEqual(distance, 3)


class NCBITestCase(unittest.TestCase):
    def _create_tax(self):
        return Taxonomy.from_ncbi("tests/data/ncbi_subset_tax.nodes.dmp", "tests/data/ncbi_subset_tax.names.dmp")

    def setUp(self) -> None:
        self.tax = self._create_tax()

    def test_root(self):
        root = self.tax.root
        self.assertEqual(root.id, "1")
        self.assertIsNone(root.parent)

    def test_find_node_by_id(self):
        node = self.tax.find_by_id("1236")
        self.assertEqual(node.id, "1236")
        self.assertEqual(node.name, "Gammaproteobacteria")
        self.assertEqual(node.parent, "1224")

        node = self.tax.find_by_id("unknown")
        self.assertIsNone(node)

    def test_index(self):
        node = self.tax["1236"]
        self.assertEqual(node.id, "1236")
        self.assertEqual(node.name, "Gammaproteobacteria")
        self.assertEqual(node.parent, "1224")

        with self.assertRaises(TaxonomyError):
            _ = self.tax["unknown"]

    def test_find_by_name(self):
        node = self.tax.find_by_name("Escherichia coli")
        self.assertEqual(node.id, "562")
        self.assertEqual(node.name, "Escherichia coli")
        self.assertEqual(node.parent, "561")

    def test_parent(self):
        parent, distance = self.tax.parent("562")
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

    def test_lca(self):
        lca = self.tax.lca("562", "91347")
        self.assertEqual(lca.id, "91347")

    def test_prune(self):
        new_tax = self.tax.prune(remove=["561"])
        self.assertIsNone(new_tax.find_by_id("561"))
        self.assertIsNone(new_tax.find_by_id("562"))
        self.assertEqual(len(new_tax), 7)

        new_tax = self.tax.prune(keep=["561"])
        self.assertEqual(len(new_tax), 8)
        self.assertIsNotNone(new_tax.find_by_id("561"))

    @unittest.skip("tax.remove doesn't work on truncated taxonomies?")
    def test_remove(self):
        tax = self._create_tax()
        tax.remove("561")
        self.assertIsNotNone(tax.find_by_id("562"))
        self.assertIsNone(tax.find_by_id("561"))
        self.assertEqual(len(tax), 8)

    def test_add(self):
        tax = self._create_tax()
        tax.add("561", "563")
        node = tax["563"]
        self.assertEqual(node.parent, "561")

    def test_edit_node(self):
        tax = self._create_tax()
        tax.edit_node("562", parent_distance=3)
        node, distance = tax.parent("562")
        self.assertEqual(distance, 3)


if __name__ == '__main__':
    unittest.main()
