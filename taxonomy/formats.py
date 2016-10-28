import networkx
from networkx.readwrite import json_graph

from taxonomy.core import Taxonomy
from taxonomy.exceptions import TaxonomyException
from taxonomy.helpers import creation_time, logger, read_json


def taxonomy_from_build_json(f):
    """
    Load a Taxonomy from a One Codex build.json file

    Args:
        f (str, file): A filepath or file handle pointing to
                       a build.json file, which must have a
                       "taxonomy" key (or an already loaded dict).

    Returns:
        Taxonomy
    """
    input_json = read_json(f).get('taxonomy')
    node_link_data = input_json['node_link_data']
    # make sure all the tax ids have been coerced into strings
    for i in range(len(node_link_data['nodes'])):
        node_link_data['nodes'][i]['id'] = str(node_link_data['nodes'][i]['id'])

    try:
        tax_graph = json_graph.node_link_graph(node_link_data)
        metadata = input_json.get('metadata', {})
    except KeyError:
        raise TaxonomyException("Improper input. Expects a JSON blob "
                                "with node_link_graph and metadata "
                                "parent nodes")
    return Taxonomy(tax_graph, metadata)


def taxonomy_from_ncbi(names, nodes, ncbi_release=None, downloaded_date=None):
    """
    Parse an NCBI nodes.dmp and names.dmp file and build a taxonomy.
    Available from: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip

    Args:
        names (str): Path to nodes.dmp file
        nodes (str): Path to names.dmp file

    Kwargs:
        names_uri (str): URI for names.dmp file (e.g., S3 URI, NCBI FTP)
        nodes_uri (str): URI for nodes.dmp file
        ncbi_release (str): NCBI release # or string
        downloaded_date (datetime): Defaults to UTC now if not specified

    Returns:
        Taxonomy: object built from NCBI data
    """
    logger.debug('Loading NCBI taxonomy from %s', names)

    # 0. Create metadata (TODO: should we calculate the MD5's again?)
    metadata = {
        'names': names,
        'nodes': nodes,
        'ncbi_release': ncbi_release,
        'downloaded_data': downloaded_date if downloaded_date is not None else creation_time(names),
    }

    # 1. Extract authoritative (scientific name type) names and map
    #    tax IDs -> scientific names
    names_dict = {}
    with open(names, mode='r') as names_file:
        names_file.readline()  # Discard header
        for i, line in enumerate(names_file):
            line = line.strip().split("\t|\t")
            tax_id = line[0].strip()
            name_txt = line[1].strip()
            name_type = line[3].rstrip("\t|")
            if name_type == 'scientific name':
                names_dict[tax_id] = name_txt

    # 2. Create a networkx object
    tax_graph = networkx.DiGraph()

    # 2. Create
    with open(nodes, mode='r') as nodes_file:
        for i, line in enumerate(nodes_file):
            line = line.split("\t|\t")
            tax_id = line[0].strip()
            parent_tax_id = line[1].strip()
            rank = line[2].strip()
            names_txt = names_dict.get(tax_id, 'No name found.')
            hidden = line[10].strip()

            attr_dict = {
                'name': names_txt,
                'rank': rank,
                'hidden': bool(hidden)
            }

            tax_graph.add_node(tax_id, attr_dict=attr_dict)
            if tax_id == parent_tax_id:
                # this is probably the root node?
                continue

            tax_graph.add_edge(tax_id, parent_tax_id)
            if i % 1000 == 0:
                logger.debug('(Line %d) Adding %s (%s) (ID: %d. Parent: %d)', i, names_txt, rank,
                             tax_id, parent_tax_id)

    return Taxonomy(tax_graph, metadata)
