from nose.tools import *
import os
import subprocess
import urllib

import taxonomy
from taxonomy import Taxonomy


# Download taxdmp.zip once to local directory
# and uncompress if not already present
def taxonomy_setup():
    if not os.path.exists(os.path.join(os.getcwd(), 'sample_data/')):
        os.mkdir(os.path.join(os.getcwd(), 'sample_data/'))
    full_tax_path = os.path.join(os.getcwd(), "sample_data/taxonomy.json.gz")
    taxdmp_path = os.path.join(os.getcwd(), "sample_data/taxdmp.zip")
    taxdmp_dir = os.path.join(os.getcwd(), "sample_data/taxdmp")
    names_head = os.path.join(os.getcwd(), "sample_data/taxdmp/names.1000.dmp")
    nodes_head = os.path.join(os.getcwd(), "sample_data/taxdmp/nodes.1000.dmp")
    ncbi_ftp = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
    if not os.path.exists(taxdmp_path):
        urllib.urlretrieve(ncbi_ftp,
                           taxdmp_path)
    if not os.path.exists(taxdmp_dir):
        res0 = subprocess.call(["unzip", taxdmp_path,
                               "-d", taxdmp_dir])
        if res0 != 0:
            raise Exception("Failed to unzip necessary taxdmp files.")

    if not os.path.exists(names_head) or not os.path.exists(nodes_head):
        res1 = subprocess.call(["head", "-n", "1000",
                               os.path.join(taxdmp_dir, "names.dmp")],
                               stdout=open(names_head, mode='w'))

        res2 = subprocess.call(["head", "-n", "1000",
                               os.path.join(taxdmp_dir, "nodes.dmp")],
                               stdout=open(nodes_head, mode='w'))

        if res1 != 0 or res2 != 0:
            raise Exception("Failed to take the head of the names.dmp or "
                            "nodes.dmp files.")

    if not os.path.exists(full_tax_path):
        tax = Taxonomy.build_from_ncbi(os.path.join(taxdmp_dir, "names.dmp"),
                                       os.path.join(taxdmp_dir, "nodes.dmp"),
                                       ncbi_ftp, ncbi_ftp,
                                       "FTP Revision - See Date",
                                       verbose=True)
        tax.save(full_tax_path)


@with_setup(taxonomy_setup)
def test_small_tax_create():
    out_tax = os.path.join(os.getcwd(), "sample_data/small_tax.json")
    names_head = os.path.join(os.getcwd(), "sample_data/taxdmp/names.1000.dmp")
    nodes_head = os.path.join(os.getcwd(), "sample_data/taxdmp/nodes.1000.dmp")
    tax = Taxonomy.build_from_ncbi(names_head, nodes_head,
                                   "Test Names", "Test Nodes",
                                   "Sample NCBI Revision")

    assert tax.__class__ is taxonomy.taxonomy.Taxonomy
    print out_tax
    tax.save(out_tax)

    # Assertion on root
    print tax.G.node[1]
    assert tax.G.node[1] == {'hidden': True, 'name': 'root',
                             'rank': 'no rank'}

    # Check that files are properly output
    assert not os.path.exists(out_tax)
    assert os.path.exists(out_tax + ".gz")  # auto-append works

    # Reread file in
    tax_b = Taxonomy.load(out_tax + ".gz")
    assert len(tax_b.G.nodes()) == len(tax.G.nodes())
    assert len(tax_b.G.edges()) == len(tax_b.G.edges())

    # Assert on new file
    assert tax_b.G.node[1] == {'hidden': True, 'name': 'root',
                               'rank': 'no rank'}
    # Final cleanup
    os.remove(out_tax + ".gz")
