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
    return Taxonomy.from_newick(
        "(Nephrops_norvegicus:0.1,Homarinus_capensis:0.2,"
        "(Homarus_americanus:0.3,Homarus_gammarus:0.4)Homarus:0.5)Homaridae;"
    )


def newick_id(tax: Taxonomy, name: str) -> str:
    (node,) = tax.find_all_by_name(name)
    return node.id


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


def test_taxonomy_iterator_is_iterable(newick_tax):
    it = iter(newick_tax)
    assert iter(it) is it


def test_json_find_all_by_name(json_tax: Taxonomy):
    assert sorted([n.id for n in json_tax.find_all_by_name("species 1.1")]) == [
        "10",
        "12",
    ]


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


@pytest.mark.parametrize("to_json", ["to_json_tree", "to_json_node_links"])
def test_from_dict_matches_from_json(json_tax: Taxonomy, to_json: str):
    serialized = getattr(json_tax, to_json)().decode()
    from_dict = Taxonomy.from_dict(json.loads(serialized))
    from_json = Taxonomy.from_json(serialized)
    assert list(from_dict) == list(from_json)
    for tax_id in from_json:
        assert from_dict[tax_id] == from_json[tax_id]
        assert from_dict[tax_id].parent == from_json[tax_id].parent
        assert from_dict[tax_id].rank == from_json[tax_id].rank


def test_from_dict_integer_ids_and_extra_data():
    tax = Taxonomy.from_dict(
        {
            "id": 1,
            "name": "root",
            "rank": "no rank",
            "readcount": 5,
            "children": [{"id": 2, "name": "child", "rank": "species", "abundance": 0.5}],
        }
    )
    assert tax["2"].parent == "1"
    assert tax["1"]["readcount"] == 5
    assert tax["2"]["abundance"] == 0.5


@pytest.mark.parametrize(
    "value",
    [
        {"id": "1", "name": "root", "children": [{"id": "1", "name": "dup"}]},
        {"name": "no id"},
        {"id": -1, "name": "root"},
        {"id": "1", "name": "root", "rank": "not a rank"},
        {"nodes": None, "links": []},
        {"nodes": [{"id": "1", "name": "root"}]},
        {"nodes": [{"id": "abc", "name": "root"}, {"id": "2"}], "links": []},
        {"nodes": [{"id": "abc", "name": "root"}]},
        {"nodes": [{"id": "1", "name": "root"}], "links": [{"source": 0}]},
        {"nodes": [{"id": "1", "name": "root"}], "links": [{"source": 5, "target": 0}]},
    ],
)
def test_from_dict_errors_match_from_json(value):
    with pytest.raises(TaxonomyError) as from_json_err:
        Taxonomy.from_json(json.dumps(value))
    with pytest.raises(TaxonomyError) as from_dict_err:
        Taxonomy.from_dict(value)
    assert str(from_dict_err.value) == str(from_json_err.value)


@pytest.mark.parametrize(
    "value,message",
    [
        ({"nodes": 5, "links": []}, "'nodes' must be an array"),
        ({"nodes": None, "links": []}, "'nodes' must be an array"),
        ({"nodes": []}, "'links' not in JSON"),
        ({"nodes": [], "links": 5}, "'links' must be an array"),
    ],
)
@pytest.mark.parametrize("load", ["from_json", "from_dict"])
def test_node_link_array_errors(value, message, load):
    arg = json.dumps(value) if load == "from_json" else value
    with pytest.raises(TaxonomyError, match=message):
        getattr(Taxonomy, load)(arg)


def test_newick_root(newick_tax: Taxonomy):
    root = newick_tax.root
    assert root.id == "0"
    assert root.name == "Homaridae"
    assert root.parent is None


def test_newick_labels_are_names_not_ids(newick_tax: Taxonomy):
    assert newick_tax.node("Homarus_gammarus") is None
    assert "Homarus_gammarus" not in newick_tax
    assert sorted(newick_tax) == ["0", "1", "2", "3", "4", "5"]


def test_newick_find_node_by_id(newick_tax: Taxonomy):
    node = newick_tax.node(newick_id(newick_tax, "Nephrops_norvegicus"))
    assert node is not None
    assert node.name == "Nephrops_norvegicus"
    assert node.parent == newick_id(newick_tax, "Homaridae")

    node = newick_tax.node(newick_id(newick_tax, "Homarus_gammarus"))
    assert node is not None
    assert node.name == "Homarus_gammarus"
    assert node.parent == newick_id(newick_tax, "Homarus")

    assert newick_tax.node("unknown") is None


def test_newick_index(newick_tax: Taxonomy):
    node = newick_tax[newick_id(newick_tax, "Nephrops_norvegicus")]
    assert node.name == "Nephrops_norvegicus"
    assert node.parent == newick_id(newick_tax, "Homaridae")

    with pytest.raises(TaxonomyError):
        _ = newick_tax["unknown"]


def test_newick_find_all_by_name(newick_tax: Taxonomy):
    assert [n.name for n in newick_tax.find_all_by_name("Homarus")] == ["Homarus"]
    assert newick_tax.find_all_by_name("Metanephrops_japonicus") == []


def test_newick_parent(newick_tax: Taxonomy):
    parent = newick_tax.parent(newick_id(newick_tax, "Homarus_gammarus"))
    assert parent is not None
    assert parent.name == "Homarus"


def test_newick_parent_with_distance(newick_tax: Taxonomy):
    parent, distance = newick_tax.parent_with_distance(
        newick_id(newick_tax, "Homarus_gammarus")
    )
    assert parent is not None
    assert distance is not None
    assert parent.name == "Homarus"
    assert abs(distance - 0.4) < 1e-6


def test_newick_children(newick_tax: Taxonomy):
    children = newick_tax.children(newick_id(newick_tax, "Homarus"))
    assert [n.name for n in children] == ["Homarus_americanus", "Homarus_gammarus"]


def test_newick_lineage(newick_tax: Taxonomy):
    lineage = newick_tax.lineage(newick_id(newick_tax, "Homarus_gammarus"))
    assert [n.name for n in lineage] == ["Homarus_gammarus", "Homarus", "Homaridae"]


def test_newick_parents(newick_tax: Taxonomy):
    parents = newick_tax.parents(newick_id(newick_tax, "Homarus_gammarus"))
    assert [n.name for n in parents] == ["Homarus", "Homaridae"]


def test_newick_lca(newick_tax: Taxonomy):
    lca = newick_tax.lca(
        newick_id(newick_tax, "Nephrops_norvegicus"),
        newick_id(newick_tax, "Homarus_gammarus"),
    )
    assert lca is not None
    assert lca.name == "Homaridae"


def test_newick_prune(newick_tax: Taxonomy):
    homarus = newick_id(newick_tax, "Homarus")
    gammarus = newick_id(newick_tax, "Homarus_gammarus")

    new_tax = newick_tax.prune(remove=[homarus])
    assert new_tax.node(gammarus) is None
    assert new_tax.node(homarus) is None
    assert len(new_tax) == 3

    new_tax = newick_tax.prune(keep=[homarus, gammarus])
    assert len(new_tax) == 3
    assert new_tax.node(newick_id(newick_tax, "Homaridae")) is not None


def test_newick_remove(newick_tax: Taxonomy):
    homarus = newick_id(newick_tax, "Homarus")
    gammarus = newick_id(newick_tax, "Homarus_gammarus")
    newick_tax.remove_node(homarus)
    assert newick_tax.node(gammarus) is not None
    assert newick_tax.node(homarus) is None
    assert len(newick_tax) == 5


def test_newick_add(newick_tax: Taxonomy):
    gammarus = newick_id(newick_tax, "Homarus_gammarus")
    newick_tax.add_node(gammarus, "100", "Homarus gammarus (Norway)", "no rank")
    assert newick_tax["100"].parent == gammarus

    newick_tax.add_node("100", "101", "Homarus gammarus (Oslofjord)", "no rank")
    assert newick_tax["101"].parent == "100"


def test_newick_edit_node(newick_tax: Taxonomy):
    gammarus = newick_id(newick_tax, "Homarus_gammarus")
    newick_tax.edit_node(gammarus, parent_distance=3)
    _, distance = newick_tax.parent_with_distance(gammarus)
    assert distance == 3


def test_newick_can_clone(newick_tax: Taxonomy):
    homarus = newick_id(newick_tax, "Homarus")
    gammarus = newick_id(newick_tax, "Homarus_gammarus")
    tax2 = newick_tax.clone()
    newick_tax.remove_node(homarus)

    assert newick_tax.node(gammarus) is not None
    assert newick_tax.node(homarus) is None
    assert len(newick_tax) == 5

    assert tax2.node(gammarus) is not None
    assert tax2.node(homarus) is not None
    assert len(tax2) == 6


def test_newick_output_uses_names(newick_tax: Taxonomy):
    res = newick_tax.to_newick().decode("utf-8")
    for name in [
        "Nephrops_norvegicus",
        "Homarinus_capensis",
        "Homarus_americanus",
        "Homarus_gammarus",
        "Homarus",
        "Homaridae",
    ]:
        assert name in res


def test_newick_output_replaces_spaces_in_names(newick_tax: Taxonomy):
    homarus = newick_id(newick_tax, "Homarus")
    newick_tax.edit_node(homarus, name="Homarus Weber 1795")
    assert "Homarus_Weber_1795" in newick_tax.to_newick().decode("utf-8")


def test_newick_round_trip(newick_tax: Taxonomy):
    tax2 = Taxonomy.from_newick(newick_tax.to_newick().decode("utf-8"))
    assert sorted(n.name for n in [tax2[i] for i in tax2]) == sorted(
        n.name for n in [newick_tax[i] for i in newick_tax]
    )
    gammarus = newick_id(tax2, "Homarus_gammarus")
    assert [n.name for n in tax2.lineage(gammarus)] == [
        "Homarus_gammarus",
        "Homarus",
        "Homaridae",
    ]


def test_newick_unlabelled_nodes_get_unique_ids():
    tax = Taxonomy.from_newick(
        "((Homarus_americanus,Homarus_gammarus),"
        "(Nephrops_norvegicus,Metanephrops_japonicus))Homaridae;"
    )
    assert len(tax) == 7
    descendant_ids = tax.descendant_ids(tax.root.id)
    assert len(descendant_ids) == len(set(descendant_ids)) == 6
    assert sorted(n.name for n in tax.descendants(tax.root.id)) == [
        "",
        "",
        "Homarus_americanus",
        "Homarus_gammarus",
        "Metanephrops_japonicus",
        "Nephrops_norvegicus",
    ]
    assert "((" in tax.to_newick().decode("utf-8")


def test_newick_duplicate_labels_are_separate_nodes():
    tax = Taxonomy.from_newick("(Homarus_americanus,Homarus_americanus)Homarus;")
    nodes = tax.find_all_by_name("Homarus_americanus")
    assert len(nodes) == 2
    assert nodes[0].id != nodes[1].id
    assert len(tax.descendant_ids(tax.root.id)) == 2


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


def test_ncbi_contains(ncbi_tax: Taxonomy):
    assert "562" in ncbi_tax
    assert "unknown" not in ncbi_tax


def test_ncbi_descendants_unknown_id(ncbi_tax: Taxonomy):
    with pytest.raises(TaxonomyError):
        ncbi_tax.descendants("unknown")


def test_ncbi_descendant_ids(ncbi_tax: Taxonomy):
    descendant_ids = ncbi_tax.descendant_ids("543")
    assert descendant_ids == [n.id for n in ncbi_tax.descendants("543")]
    assert "562" in descendant_ids
    assert "543" not in descendant_ids


def test_ncbi_descendants_are_depth_first(ncbi_tax: Taxonomy):
    seen = {"543"}
    for node in ncbi_tax.descendants("543"):
        assert node.parent in seen
        seen.add(node.id)


def test_ncbi_descendant_ids_unknown_id(ncbi_tax: Taxonomy):
    with pytest.raises(TaxonomyError):
        ncbi_tax.descendant_ids("unknown")


def test_ncbi_lineage(ncbi_tax: Taxonomy):
    lineage = ncbi_tax.lineage("562")
    assert len(lineage) == 9
    assert lineage[0].id == "562"
    assert lineage[1].id == "561"
    assert lineage[-1].id == "1"


def test_ncbi_lineage_ids(ncbi_tax: Taxonomy):
    lineage_ids = ncbi_tax.lineage_ids("562")
    assert lineage_ids == [n.id for n in ncbi_tax.lineage("562")]
    assert lineage_ids[0] == "562"
    assert lineage_ids[-1] == "1"


def test_ncbi_lineage_ids_unknown_id(ncbi_tax: Taxonomy):
    with pytest.raises(TaxonomyError):
        ncbi_tax.lineage_ids("unknown")


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
        tax["562"].__repr__()
        == '<TaxonomyNode (id="562" rank="species" name="Escherichia coli")>'
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
    not os.getenv("TAXONOMY_TEST_NCBI"),
    reason="Define TAXONOMY_TEST_NCBI to run NCBI test",
)
def test_latestncbi_load_latest_ncbi_taxonomy():
    download("https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
    subprocess.check_output(["tar", "-zxvf", "taxdump.tar.gz"])
    Taxonomy.from_ncbi(".")
