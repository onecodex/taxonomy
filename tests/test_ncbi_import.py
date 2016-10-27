import os

import taxonomy
from taxonomy import Taxonomy
from taxonomy.formats import taxonomy_from_ncbi


# Download taxdmp.zip once to local directory and uncompress if not already present
# we're storing this 1000 short list in the repo now so we don't have to do this every time
# (but we might still want to expose a "download NCBI" sequence function at some point
# import subprocess
# import urllib
#
# import pytest
# @pytest.fixture
# def taxonomy_setup():
#     if not os.path.exists(os.path.join(os.getcwd(), 'sample_data/')):
#         os.mkdir(os.path.join(os.getcwd(), 'sample_data/'))
#     full_tax_path = os.path.join(os.getcwd(), "sample_data/taxonomy.json.gz")
#     taxdmp_path = os.path.join(os.getcwd(), "sample_data/taxdmp.zip")
#     taxdmp_dir = os.path.join(os.getcwd(), "sample_data/taxdmp")
#     names_head = os.path.join(os.getcwd(), "sample_data/taxdmp/names.1000.dmp")
#     nodes_head = os.path.join(os.getcwd(), "sample_data/taxdmp/nodes.1000.dmp")
#     ncbi_ftp = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
#     if not os.path.exists(taxdmp_path):
#         urllib.urlretrieve(ncbi_ftp,
#                            taxdmp_path)
#     if not os.path.exists(taxdmp_dir):
#         res0 = subprocess.call(['unzip', taxdmp_path, '-d', taxdmp_dir])
#         if res0 != 0:
#             raise Exception('Failed to unzip necessary taxdmp files.')
#
#     if not os.path.exists(names_head) or not os.path.exists(nodes_head):
#         res1 = subprocess.call(['head", "-n", "1000',
#                                os.path.join(taxdmp_dir, 'names.dmp')],
#                                stdout=open(names_head, mode='w'))
#
#         res2 = subprocess.call(["head", "-n", "1000",
#                                os.path.join(taxdmp_dir, 'nodes.dmp')],
#                                stdout=open(nodes_head, mode='w'))
#
#         if res1 != 0 or res2 != 0:
#             raise Exception("Failed to take the head of the names.dmp or "
#                             "nodes.dmp files.")
#
#     if not os.path.exists(full_tax_path):
#         tax = taxonomy_from_ncbi(os.path.join(taxdmp_dir, 'names.dmp'),
#                                  os.path.join(taxdmp_dir, 'nodes.dmp'),
#                                  'FTP Revision - See Date')
#         tax.save(full_tax_path)


def test_small_tax_create():
    out_tax = os.path.join(os.getcwd(), 'sample_data/small_tax.json')
    names_head = os.path.join(os.getcwd(), 'tests/sample_ncbi_data/names.1000.dmp')
    nodes_head = os.path.join(os.getcwd(), 'tests/sample_ncbi_data/nodes.1000.dmp')
    tax = taxonomy_from_ncbi(names_head, nodes_head, 'Sample NCBI Revision')

    assert tax.__class__ is taxonomy.core.Taxonomy
    tax.save(out_tax)

    # Assertion on root
    assert tax.tax_graph.node['1'] == {'hidden': True, 'name': 'root', 'rank': 'no rank'}

    # Make sure the root doesn't link to itself (regression test)
    assert not tax.tax_graph['1']

    # Check that files are properly output
    assert not os.path.exists(out_tax)
    assert os.path.exists(out_tax + '.gz')  # auto-append works

    # Reread file in
    tax_b = Taxonomy.load(out_tax + '.gz')
    assert len(tax_b.tax_graph.nodes()) == len(tax.tax_graph.nodes())
    assert len(tax_b.tax_graph.edges()) == len(tax.tax_graph.edges())

    # Assert on new file
    assert tax_b.tax_graph.node['1'] == {'hidden': True, 'name': 'root', 'rank': 'no rank'}
    # Final cleanup
    os.remove(out_tax + '.gz')
