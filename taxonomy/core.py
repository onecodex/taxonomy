"""
Main Taxonomy object/class code.
"""
import gzip
import itertools
import networkx
import os
try:
    import simplejson as json
except ImportError:
    import json

from networkx.readwrite import json_graph

from taxonomy.exceptions import TaxonomyException
from taxonomy.helpers import read_json


class Taxonomy(object):
    """
    networkx-based taxonomy object. Main object type of `taxonomy` module.

    Args:
        tax_graph (networkx.DiGraph): Input Taxonomy graph
        metadata (dict): Associated metadata

    Returns:
        Taxonomy object
    """
    def __init__(self, tax_graph, metadata):
        self.tax_graph = tax_graph
        self.metadata = metadata

    def __repr__(self):
        return '<Taxonomy ({} nodes)>'.format(self.tax_graph.number_of_nodes())

    @staticmethod
    def load(f):
        """
        Load a Taxonomy from a .json file

        Args:
            f (str, file): A filepath or file handle pointing to
                           a taxonomy.json file (or an already loaded dict)

        Returns:
            Taxonomy
        """
        input_json = read_json(f)
        try:
            tax_graph = json_graph.node_link_graph(input_json['node_link_data'])
            metadata = input_json['metadata']
        except KeyError:
            raise TaxonomyException("Improper input. Expects a JSON blob "
                                    "with node_link_graph and metadata "
                                    "parent nodes")
        return Taxonomy(tax_graph, metadata)

    def save(self, f, compress=True):
        """
        Save a Taxonomy object out to a taxonomy.json file

        Args:
            f (str, file): A filepath or file handle

        Kwargs:
            compress (bool): Gzip the output? Default is True.
                             Note: Is only used if a filepath is
                             passed.

        Returns:
            None
        """
        out = {}
        out['node_link_data'] = json_graph.node_link_data(self.tax_graph)
        out['metadata'] = self.metadata
        if isinstance(f, (file, gzip.GzipFile)):
            json.dump(out, f)
        else:
            if gzip:
                if os.path.splitext(f)[1] != '.gz':
                    f = f + '.gz'
                json.dump(out, gzip.open(f, mode='w'))
            else:
                json.dump(out, open(f, mode='w'))

    def tax_info(self, tax_id):
        if not self.tax_graph.has_node(tax_id):
            raise KeyError('Tax ID {} not in taxonomy'.format(tax_id))

        return self.tax_graph.node[tax_id]

    def parent_at_rank(self, tax_id, rank='species'):
        if self.tax_graph.node.get(tax_id, {}).get('rank') == rank:
            return tax_id

        for parent_id in self.parents(tax_id):
            if self.tax_graph.node.get(parent_id, {}).get('rank') == rank:
                return parent_id

        # nothing at this level or above is the right rank; return None
        return None

    def children(self, tax_id, all_children=False):
        if not self.tax_graph.has_node(tax_id):
            raise KeyError('Tax ID {} not in taxonomy'.format(tax_id))

        if all_children:
            return networkx.ancestors(self.tax_graph, tax_id)
        else:
            return self.tax_graph.predecessors(tax_id)

    def parents(self, tax_id):
        if not self.tax_graph.has_node(tax_id):
            raise KeyError('Tax ID {} not in taxonomy'.format(tax_id))

        return list(networkx.dfs_preorder_nodes(self.tax_graph, tax_id))[1:]

    def lowest_common_ancestor(self, *tax_ids):
        len_tax_ids = len(tax_ids)
        if len_tax_ids == 0:
            return None
        elif len_tax_ids == 1:
            return tax_ids[0]

        if not all(self.tax_graph.has_node(tax_id) for tax_id in tax_ids):
            raise KeyError('Tax IDs {} not in taxonomy'.format(','.join(str(t) for t in tax_ids)))

        return reduce(self.lowest_common_ancestor_double, tax_ids)

    def lowest_common_ancestor_double(self, tax_id_1, tax_id_2):
        parent_gen_1 = networkx.dfs_postorder_nodes(self.tax_graph, tax_id_1)
        parent_gen_2 = networkx.dfs_postorder_nodes(self.tax_graph, tax_id_2)
        parent_id = None
        for ptax_id_1, ptax_id_2 in itertools.izip(parent_gen_1, parent_gen_2):
            if ptax_id_1 != ptax_id_2:
                break
            parent_id = ptax_id_1

        # the tax_ids don't share a common parent; return None
        return parent_id
