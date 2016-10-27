import pytest
import networkx

from taxonomy import Taxonomy


@pytest.fixture
def sample_taxonomy():
    tax_graph = networkx.DiGraph()
    tax_graph.add_edges_from([
        ('562', '561'),
        ('561', '543'),
        ('543', '91347'),
        ('91347', '1236'),
        ('1236', '1224'),
        ('1224', '2'),
        ('2', '131567'),
        ('131567', '1'),
        ('1219', '1218'),
        ('1218', '1213'),
        ('1213', '1890424'),
        ('1890424', '1117'),
        ('1117', '1798711'),
        ('1798711', '1783272'),
        ('1783272', '2'),
        ('1080230', '1148'),
        ('1080229', '1148'),
        ('1148', '1142'),
        ('1142', '1890428'),
        ('1890428', '1890424'),
    ])
    tax_graph.add_nodes_from([
        ('562', {'name': 'Escherichia coli', 'rank': 'species'}),
        ('561', {'name': 'Escherichia', 'rank': 'genus'}),
        ('543', {'name': 'Enterobacteriaceae', 'rank': 'family'}),
        ('91347', {'name': 'Enterobacterales', 'rank': 'order'}),
        ('1236', {'name': 'Gammaproteobacteria', 'rank': 'class'}),
        ('1224', {'name': 'Proteobacteria', 'rank': 'phylum'}),
        ('2', {'name': 'Bacteria', 'rank': 'superkingdom'}),
        ('131567', {'name': 'cellular organisms', 'rank': 'no rank'}),
        ('1', {'name': 'root', 'rank': 'no rank'}),
        ('1219', {'name': 'Prochlorococcus marinus', 'rank': 'species'}),
        ('1218', {'name': 'Prochlorococcus', 'rank': 'genus'}),
        ('1213', {'name': 'Prochloraceae', 'rank': 'family'}),
        ('1890424', {'name': 'Synechococcales', 'rank': 'order'}),
        ('1117', {'name': 'Cyanobacteria', 'rank': 'phylum'}),
        ('1798711', {'name': 'Cyanobacteria/Melainabacteria group', 'rank': 'no rank'}),
        ('1783272', {'name': 'Terrabacteria group', 'rank': 'no rank'}),
        ('1080229', {'name': 'Synechocystis sp. PCC 6803 substr. PCC-N', 'rank': 'no rank'}),
        ('1080230', {'name': 'Synechocystis sp. PCC 6803 substr. PCC-P', 'rank': 'no rank'}),
        ('1148', {'name': 'Synechocystis sp. PCC 6803', 'rank': 'species'}),
        ('1142', {'name': 'Synechocystis', 'rank': 'genus'}),
        ('1890428', {'name': 'Merismopediaceae', 'rank': 'family'}),
    ])
    return Taxonomy(tax_graph, {})


def test_parents(sample_taxonomy):
    # e. coli has the right parents
    assert sample_taxonomy.parents('562') == ['561', '543', '91347', '1236', '1224', '2', '131567',
                                              '1']

    # root has no parents
    assert sample_taxonomy.parents('1') == []

    # a non-existent tax_id will cause an error
    with pytest.raises(KeyError):
        sample_taxonomy.parents('27')


def test_children(sample_taxonomy):
    # root should have 1 immediate child (for now) and all other nodes as it's "all" children
    assert len(sample_taxonomy.children('1')) == 1
    n_nodes = len(sample_taxonomy.tax_graph.nodes())
    assert len(sample_taxonomy.children('1', all_children=True)) == n_nodes - 1


def test_lowest_common_ancestor(sample_taxonomy):
    # different child nodes resolve to the right parent node
    assert sample_taxonomy.lowest_common_ancestor('1080229', '1080230') == '1148'
    assert sample_taxonomy.lowest_common_ancestor('1080229', '1219') == '1890424'
    assert sample_taxonomy.lowest_common_ancestor('1080229', '562') == '2'

    assert sample_taxonomy.lowest_common_ancestor('1080229', '1080230', '1219', '562') == '2'

    # a child and a parent resolve to the parent
    assert sample_taxonomy.lowest_common_ancestor('562', '2') == '2'
    assert sample_taxonomy.lowest_common_ancestor('1080229', '1117') == '1117'

    # a node and a non-existent id will throw an error
    with pytest.raises(KeyError):
        sample_taxonomy.lowest_common_ancestor('1080229', '27')


def test_parent_at_rank(sample_taxonomy):
    # a species is itself
    assert sample_taxonomy.parent_at_rank('562') == '562'

    # map e. coli back to proteobacteria
    assert sample_taxonomy.parent_at_rank('562', 'phylum') == '1224'

    # strains find their species
    assert sample_taxonomy.parent_at_rank('1080229') == '1148'
