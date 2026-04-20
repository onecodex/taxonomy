"""Tests for taxonomy data storage and tree aggregation functions.

Taxonomy: Cetacea, as catalogued by Ishmael in Moby Dick.
"Sightings" data represents whale encounters recorded in the novel.

Tree structure:
    Cetacea
    ├── Mysticeti (baleen whales)
    │   ├── Eubalaena_glacialis       (right whale,     5 sightings)
    │   └── Balaenopteridae
    │       ├── Balaenoptera_musculus  (blue whale,      2 sightings)
    │       └── Megaptera_novaeangliae (humpback whale,  3 sightings)
    └── Odontoceti (toothed whales)
        ├── Physeter_macrocephalus    (sperm whale,     90 sightings)
        └── Kogia_breviceps           (pygmy sperm,      1 sighting)
"""

import pytest
from taxonomy import Taxonomy

WHALE_NEWICK = (
    "("
    "(Eubalaena_glacialis:1,"
    "(Balaenoptera_musculus:1,Megaptera_novaeangliae:1)Balaenopteridae:1"
    ")Mysticeti:1,"
    "(Physeter_macrocephalus:1,Kogia_breviceps:1)Odontoceti:1"
    ")Cetacea;"
)

SPECIES_SIGHTINGS = {
    "Physeter_macrocephalus": 90,  # the white whale — dominates the novel
    "Eubalaena_glacialis": 5,
    "Balaenoptera_musculus": 2,
    "Megaptera_novaeangliae": 3,
    "Kogia_breviceps": 1,
}

TOTAL_SIGHTINGS = sum(SPECIES_SIGHTINGS.values())  # 101


@pytest.fixture
def whale_tax():
    """Whale taxonomy with sighting counts from Moby Dick."""
    tax = Taxonomy.from_newick(WHALE_NEWICK)
    for node_id in SPECIES_SIGHTINGS:
        tax.edit_node(node_id, rank="species")
    for node_id, count in SPECIES_SIGHTINGS.items():
        tax.set_data(node_id, "sightings", count)
    return tax


# ── set_data ─────────────────────────────────────────────────────────────────


def test_set_data(whale_tax):
    assert whale_tax["Physeter_macrocephalus"]["sightings"] == 90


def test_set_data_overwrite(whale_tax):
    whale_tax.set_data("Physeter_macrocephalus", "sightings", 999)
    assert whale_tax["Physeter_macrocephalus"]["sightings"] == 999


def test_set_data_new_key(whale_tax):
    whale_tax.set_data("Physeter_macrocephalus", "chapters_mentioned", 135)
    assert whale_tax["Physeter_macrocephalus"]["chapters_mentioned"] == 135


# ── node.get ─────────────────────────────────────────────────────────────────


def test_node_get_existing_key(whale_tax):
    assert whale_tax["Physeter_macrocephalus"].get("sightings") == 90


def test_node_get_missing_key_returns_none(whale_tax):
    # Internal nodes have no sightings set
    assert whale_tax["Mysticeti"].get("sightings") is None


def test_node_get_missing_key_custom_default(whale_tax):
    assert whale_tax["Mysticeti"].get("sightings", 0) == 0


# ── node.data ─────────────────────────────────────────────────────────────────


def test_node_data_returns_dict(whale_tax):
    data = whale_tax["Physeter_macrocephalus"].data
    assert isinstance(data, dict)
    assert data["sightings"] == 90


def test_node_data_empty_for_internal_node(whale_tax):
    assert whale_tax["Mysticeti"].data == {}


def test_node_data_is_snapshot(whale_tax):
    """TaxonomyNode is a snapshot — set_data does not update existing references."""
    node = whale_tax["Physeter_macrocephalus"]
    whale_tax.set_data("Physeter_macrocephalus", "sightings", 999)
    assert node["sightings"] == 90  # old snapshot
    assert whale_tax["Physeter_macrocephalus"]["sightings"] == 999  # re-fetched


# ── reduce_up ─────────────────────────────────────────────────────────────────


def clade_sightings(node, child_results):
    return node.get("sightings", 0) + sum(child_results)


def test_reduce_up_leaf(whale_tax):
    """Leaf clade sightings equals its own sightings."""
    annotated = whale_tax.reduce_up(
        "Physeter_macrocephalus", "clade_sightings", clade_sightings
    )
    assert annotated["Physeter_macrocephalus"]["clade_sightings"] == 90


def test_reduce_up_internal_node(whale_tax):
    """Odontoceti clade = sperm whale + pygmy sperm whale."""
    annotated = whale_tax.reduce_up(
        "Cetacea", "clade_sightings", clade_sightings
    )
    assert annotated["Odontoceti"]["clade_sightings"] == 90 + 1


def test_reduce_up_mysticeti(whale_tax):
    """Mysticeti clade = right whale + blue whale + humpback."""
    annotated = whale_tax.reduce_up(
        "Cetacea", "clade_sightings", clade_sightings
    )
    assert annotated["Mysticeti"]["clade_sightings"] == 5 + 2 + 3


def test_reduce_up_root(whale_tax):
    """Root clade sightings equals total across all species."""
    annotated = whale_tax.reduce_up(
        "Cetacea", "clade_sightings", clade_sightings
    )
    assert annotated["Cetacea"]["clade_sightings"] == TOTAL_SIGHTINGS


def test_reduce_up_preserves_original(whale_tax):
    """reduce_up returns a new taxonomy; original is unchanged."""
    whale_tax.reduce_up("Cetacea", "clade_sightings", clade_sightings)
    assert whale_tax["Odontoceti"].get("clade_sightings") is None


def test_reduce_up_count_species(whale_tax):
    """Count species with any sightings per clade."""
    annotated = whale_tax.reduce_up(
        "Cetacea",
        "detected_species",
        lambda node, child_results: sum(child_results)
        + (1 if node.rank == "species" and node.get("sightings", 0) > 0 else 0),
    )
    assert annotated["Mysticeti"]["detected_species"] == 3
    assert annotated["Odontoceti"]["detected_species"] == 2
    assert annotated["Cetacea"]["detected_species"] == 5


def test_reduce_up_max_subclade(whale_tax):
    """Max clade sightings among all subclades — uses child_results explicitly."""
    annotated = whale_tax.reduce_up(
        "Cetacea", "clade_sightings", clade_sightings
    )
    annotated = annotated.reduce_up(
        "Cetacea",
        "dominant_subclade_sightings",
        lambda node, child_results: max(child_results)
        if child_results
        else node.get("clade_sightings", 0),
    )
    # Physeter (90) is the dominant individual species; that propagates up through Odontoceti
    assert annotated["Cetacea"]["dominant_subclade_sightings"] == 90


# ── map_down ──────────────────────────────────────────────────────────────────


def test_map_down_depth_root(whale_tax):
    annotated = whale_tax.map_down(
        "Cetacea", "depth", 0, lambda parent_depth, node: parent_depth + 1
    )
    assert annotated["Cetacea"]["depth"] == 1


def test_map_down_depth_internal(whale_tax):
    annotated = whale_tax.map_down(
        "Cetacea", "depth", 0, lambda parent_depth, node: parent_depth + 1
    )
    assert annotated["Odontoceti"]["depth"] == 2
    assert annotated["Balaenopteridae"]["depth"] == 3


def test_map_down_depth_leaf(whale_tax):
    annotated = whale_tax.map_down(
        "Cetacea", "depth", 0, lambda parent_depth, node: parent_depth + 1
    )
    assert annotated["Physeter_macrocephalus"]["depth"] == 3
    assert annotated["Balaenoptera_musculus"]["depth"] == 4


def test_map_down_lineage(whale_tax):
    annotated = whale_tax.map_down(
        "Cetacea",
        "lineage",
        "",
        lambda parent, node: f"{parent};{node.id}" if parent else node.id,
    )
    assert annotated["Physeter_macrocephalus"]["lineage"] == (
        "Cetacea;Odontoceti;Physeter_macrocephalus"
    )
    assert annotated["Balaenoptera_musculus"]["lineage"] == (
        "Cetacea;Mysticeti;Balaenopteridae;Balaenoptera_musculus"
    )


def test_map_down_preserves_original(whale_tax):
    """map_down returns a new taxonomy; original is unchanged."""
    whale_tax.map_down("Cetacea", "depth", 0, lambda d, node: d + 1)
    assert whale_tax["Cetacea"].get("depth") is None


# ── chaining ──────────────────────────────────────────────────────────────────


def test_chain_reduce_up_then_map_down(whale_tax):
    """
    Compute clade sightings bottom-up, then propagate clade fraction top-down.
    Each node stores what fraction of the root's sightings belong to its clade.
    """
    annotated = whale_tax.reduce_up(
        "Cetacea", "clade_sightings", clade_sightings
    )
    annotated = annotated.map_down(
        "Cetacea",
        "clade_fraction",
        1.0,
        lambda _, node: node["clade_sightings"] / TOTAL_SIGHTINGS,
    )

    assert annotated["Cetacea"]["clade_fraction"] == pytest.approx(1.0)
    assert annotated["Odontoceti"]["clade_fraction"] == pytest.approx(91 / 101)
    assert annotated["Mysticeti"]["clade_fraction"] == pytest.approx(10 / 101)
