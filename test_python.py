import json
import pytest

from taxonomy import Taxonomy, TaxonomyError
from downloads import download
import os
import subprocess

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


@pytest.fixture
def json_tax():
    return Taxonomy.from_json(JSON_DATA)


@pytest.fixture
def newick_tax():
    return Taxonomy.from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")


@pytest.fixture
def ncbi_tax():
    return Taxonomy.from_ncbi("tests/data/")


@pytest.fixture
def gtdb_tax():
    with open("tests/data/gtdb_sample.tsv") as file:
        return Taxonomy.from_gtdb(file.read())


def test_json_internal_index(json_tax: Taxonomy):
    assert [
        json_tax.internal_index(x)
        for x in ["1", "9", "2", "11", "8", "5", "3", "4", "6", "7", "10"]
    ] == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]


def test_json_find_all_by_name(json_tax: Taxonomy):
    assert sorted([n.id for n in json_tax.find_all_by_name("species 1.1")]) == ["10", "12"]


def test_json_edit_node_parent_updates_children(json_tax: Taxonomy):
    assert json_tax["5"].parent == "4"
    json_tax.edit_node("5", parent_id="1")
    node = json_tax["5"]
    assert node.parent == "1"
    assert "5" not in {n.id for n in json_tax.children("4")}
    assert "5" in {n.id for n in json_tax.children("1")}


def test_json_prune_works_after_editing_tree(json_tax: Taxonomy):
    tax = json_tax.clone()
    tax.edit_node("5", parent_id="1")
    pruned = tax.prune(keep=["5"])
    assert pruned["5"].parent == "1"


def test_json_to_json_tree(json_tax: Taxonomy):
    small_tax = json_tax.prune(remove=[str(i) for i in range(3, 12)])
    actual = json.loads(small_tax.to_json_tree())
    expected = {
        "id": "1",
        "name": "root",
        "rank": "no rank",
        "children": [
            {
                "id": "2",
                "name": "superkingdom 1",
                "rank": "superkingdom",
                "children": [],
            }
        ],
    }
    assert actual == expected


def test_json_to_json_tree_with_empty_tree(json_tax):
    empty_tax = json_tax.prune(keep=[])
    with pytest.raises(TaxonomyError):
        empty_tax.to_json_tree()


def test_json_to_json_node_links_empty_tree(json_tax):
    empty_tax = json_tax.prune(keep=[])
    actual = json.loads(empty_tax.to_json_node_links())
    expected = {
        "directed": True,
        "graph": [],
        "links": [],
        "multigraph": False,
        "nodes": [],
    }
    assert actual == expected


def test_newick_root(newick_tax: Taxonomy):
    root = newick_tax.root
    assert root.id == "F"
    assert root.parent is None


def test_newick_find_node_by_id(newick_tax: Taxonomy):
    node = newick_tax.node("A")
    assert node == newick_tax.node("A")
    assert node is not None
    assert node.id == "A"
    assert node.parent == "F"

    node = newick_tax.node("D")
    assert node is not None
    assert node.id == "D"
    assert node.parent == "E"

    node = newick_tax.node("unknown")
    assert node is None


def test_newick_index(newick_tax: Taxonomy):
    node = newick_tax["A"]
    assert node.id == "A"
    assert node.parent == "F"

    with pytest.raises(TaxonomyError):
        _ = newick_tax["unknown"]


def test_newick_find_all_by_name(newick_tax: Taxonomy):
    nodes = newick_tax.find_all_by_name("A")
    assert nodes == []


def test_newick_parent(newick_tax: Taxonomy):
    parent = newick_tax.parent("D")
    assert parent is not None
    assert parent.id == "E"


def test_newick_parent_with_distance(newick_tax: Taxonomy):
    parent, distance = newick_tax.parent_with_distance("D")
    assert parent is not None
    assert distance is not None
    assert parent.id == "E"
    assert abs(distance - 0.4) < 1e-6


def test_newick_children(newick_tax: Taxonomy):
    children = newick_tax.children("E")
    assert len(children) == 2
    assert children[0].id == "C"
    assert children[1].id == "D"


def test_newick_lineage(newick_tax: Taxonomy):
    lineage = newick_tax.lineage("D")
    assert len(lineage) == 3
    assert lineage[0].id == "D"
    assert lineage[1].id == "E"
    assert lineage[2].id == "F"


def test_newick_parents(newick_tax: Taxonomy):
    lineage = newick_tax.parents("D")
    assert len(lineage) == 2
    assert lineage[0].id == "E"
    assert lineage[1].id == "F"


def test_newick_lca(newick_tax: Taxonomy):
    lca = newick_tax.lca("A", "D")
    assert lca is not None
    assert lca.id == "F"


def test_newick_prune(newick_tax: Taxonomy):
    new_tax = newick_tax.prune(remove=["E"])
    assert new_tax.node("D") is None
    assert new_tax.node("E") is None
    assert len(new_tax) == 3

    new_tax = newick_tax.prune(keep=["E", "D"])
    assert len(new_tax) == 3
    assert new_tax.node("F") is not None


def test_newick_remove(newick_tax: Taxonomy):
    newick_tax.remove_node("E")
    assert newick_tax.node("D") is not None
    assert newick_tax.node("E") is None
    assert len(newick_tax) == 5


def test_newick_add(newick_tax: Taxonomy):
    newick_tax.add_node("D", "G", "something", "species")
    node = newick_tax["G"]
    assert node.parent == "D"

    newick_tax.add_node("G", "H", "something else", "species")
    node = newick_tax["H"]
    assert node.parent == "G"


def test_newick_edit_node(newick_tax: Taxonomy):
    newick_tax.edit_node("D", parent_distance=3)
    _, distance = newick_tax.parent_with_distance("D")
    assert distance == 3


def test_newick_can_clone(newick_tax: Taxonomy):
    tax2 = newick_tax.clone()
    newick_tax.remove_node("E")

    assert newick_tax.node("D") is not None
    assert newick_tax.node("E") is None
    assert len(newick_tax) == 5

    assert tax2.node("D") is not None
    assert tax2.node("E") is not None
    assert len(tax2) == 6


def test_newick_output_uses_tax_ids(newick_tax: Taxonomy):
    res = newick_tax.to_newick().decode("utf-8")
    for tax_id in ["A", "B", "C", "D", "E", "F"]:
        assert tax_id in res


def test_ncbi_root(ncbi_tax: Taxonomy):
    root = ncbi_tax.root
    assert root.id == "1"
    assert root.parent is None


def test_ncbi_find_node_by_id(ncbi_tax: Taxonomy):
    node = ncbi_tax.node("1236")
    assert node is not None
    assert node.id == "1236"
    assert node.name == "Gammaproteobacteria"
    assert node.parent == "1224"

    node = ncbi_tax.node("unknown")
    assert node is None


def test_ncbi_index(ncbi_tax: Taxonomy):
    node = ncbi_tax["1236"]
    assert node.id == "1236"
    assert node.name == "Gammaproteobacteria"
    assert node.parent == "1224"

    with pytest.raises(TaxonomyError):
        _ = ncbi_tax["unknown"]


def test_ncbi_find_all_by_name(ncbi_tax: Taxonomy):
    nodes = ncbi_tax.find_all_by_name("Escherichia coli")
    assert [n.id for n in nodes] == ["562"]
    assert [n.name for n in nodes] == ["Escherichia coli"]
    assert [n.parent for n in nodes] == ["561"]


def test_ncbi_parent(ncbi_tax: Taxonomy):
    parent = ncbi_tax.parent("562")
    assert parent is not None
    assert parent.id == "561"


def test_ncbi_parent_with_distance(ncbi_tax: Taxonomy):
    parent, distance = ncbi_tax.parent_with_distance("562")
    assert parent is not None
    assert distance is not None
    assert parent.id == "561"
    assert abs(distance - 1.0) < 1e-6


def test_ncbi_children(ncbi_tax: Taxonomy):
    children = ncbi_tax.children("561")
    assert len(children) == 1
    assert children[0].id == "562"


def test_ncbi_lineage(ncbi_tax: Taxonomy):
    lineage = ncbi_tax.lineage("562")
    assert len(lineage) == 9
    assert lineage[0].id == "562"
    assert lineage[1].id == "561"
    assert lineage[-1].id == "1"


def test_ncbi_parents(ncbi_tax: Taxonomy):
    lineage = ncbi_tax.parents("562")
    assert len(lineage) == 8
    assert lineage[0].id == "561"
    assert lineage[-1].id == "1"


def test_ncbi_lca(ncbi_tax: Taxonomy):
    lca = ncbi_tax.lca("562", "91347")
    assert lca is not None
    assert lca.id == "91347"


def test_ncbi_prune(ncbi_tax: Taxonomy):
    new_tax = ncbi_tax.prune(remove=["561"])
    assert new_tax.node("561") is None
    assert new_tax.node("562") is None
    assert len(new_tax) == 8

    new_tax = ncbi_tax.prune(keep=["561"])
    assert len(new_tax) == 8
    assert new_tax.node("561") is not None


def test_taxonomy_is_iterable(ncbi_tax: Taxonomy):
    for _ in ncbi_tax:
        break


@pytest.mark.skip(reason="tax.remove doesn't work on truncated taxonomies?")
def test_ncbi_remove():
    tax = Taxonomy.from_ncbi("tests/data/")
    tax.remove_node("561")
    assert tax.node("562") is not None
    assert tax.node("561") is None
    assert len(tax) == 8


def test_ncbi_add():
    tax = Taxonomy.from_ncbi("tests/data/")
    tax.add_node("561", "563", "Listeria", "species")
    node = tax["563"]
    assert node.parent == "561"
    assert node.name == "Listeria"
    assert node.rank == "species"

    tax.add_node("563", "100000001", "Pizzeria", "genus")
    node = tax["100000001"]
    assert node.parent == "563"
    assert node.name == "Pizzeria"
    assert node.rank == "genus"


def test_ncbi_cannot_add_duplicate_tax_id():
    tax = Taxonomy.from_ncbi("tests/data/")
    tax.add_node("561", "563", "Listeria", "species")

    with pytest.raises(TaxonomyError) as context:
        tax.add_node("561", "563", "Listeria", "species")
    assert "563" in str(context.value)


def test_ncbi_edit_node():
    tax = Taxonomy.from_ncbi("tests/data/")
    tax.edit_node("562", parent_distance=3)
    _, distance = tax.parent_with_distance("562")
    assert distance == 3


def test_ncbi_edit_node_parent():
    tax = Taxonomy.from_ncbi("tests/data/")
    assert tax["562"].parent == "561"
    tax.edit_node("562", parent_id="1")
    assert tax["562"].parent == "1"


def test_ncbi_repr():
    tax = Taxonomy.from_ncbi("tests/data/")
    assert (
        tax["562"].__repr__() == '<TaxonomyNode (id="562" rank="species" name="Escherichia coli")>'
    )


def test_gtdb_root(gtdb_tax: Taxonomy):
    root = gtdb_tax.root
    assert root.id == "d__Bacteria"
    assert root.rank == "domain"
    assert root.parent is None


def test_gtdb_lineage(gtdb_tax: Taxonomy):
    assert [n.id for n in gtdb_tax.lineage("d__Bacteria")] == ["d__Bacteria"]
    assert [n.id for n in gtdb_tax.lineage("c__Bacilli")] == [
        "c__Bacilli",
        "p__Firmicutes",
        "d__Bacteria",
    ]
    assert [n.id for n in gtdb_tax.lineage("s__Escherichia coli")] == [
        "s__Escherichia coli",
        "g__Escherichia",
        "f__Enterobacteriaceae",
        "o__Enterobacterales",
        "c__Gammaproteobacteria",
        "p__Proteobacteria",
        "d__Bacteria",
    ]


def test_gtdb_invalid_format():
    with open("tests/data/gtdb_invalid.tsv") as file:
        with pytest.raises(TaxonomyError):
            Taxonomy.from_gtdb(file.read())


@pytest.mark.skipif(
    not os.getenv("TAXONOMY_TEST_NCBI"), reason="Define TAXONOMY_TEST_NCBI to run NCBI test"
)
def test_latestncbi_load_latest_ncbi_taxonomy():
    download("https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
    subprocess.check_output(["tar", "-zxvf", "taxdump.tar.gz"])
    Taxonomy.from_ncbi(".")
