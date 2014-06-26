"""
Main Taxonomy object/class code.
"""
import datetime
import gzip
import networkx as nx
from networkx.readwrite import json_graph
import os
try:
    import simplejson as json
except ImportError:
    import json

from .helpers import md5_for_file


class TaxonomyException(Exception):
    pass


class TaxonomyMetadata(object):
    """
    A simple metadata object for the Taxonomy.

    Args:
        names_file (str): Filepath to local names.dmp file
        nodes_file (str): Filepath to local nodes.dmp file
        names_uri (str): URI for names.dmp file (e.g., S3 URI, NCBI FTP)
        nodes_uri (str): URI for nodes.dmp file
        ncbi_release (str): NCBI release # or string

    Kwargs:
        downloaded_date (datetime): Defaults to UTC now if not specified
        names_md5 (str): Precomputed MD5 checksum for input file.
                         Only use if known. Will otherwise be calculated.
        nodes_md5 (str): Same as names_md5 but for nodes.dmp
    """
    def __init__(self, names_file, nodes_file, names_uri, nodes_uri,
                 ncbi_release, downloaded_date=None,
                 names_md5=None, nodes_md5=None):
        self.names_uri = names_uri
        self.nodes_uri = nodes_uri
        if names_md5 is not None:
            self.names_md5 = names_md5
        elif names_file is not None:
            self.names_md5 = md5_for_file(names_file)
        else:
            raise TaxonomyException("Need to pass an md5 hash or a file.")

        if nodes_md5 is not None:
            self.nodes_md5 = nodes_md5
        elif nodes_file is not None:
            self.nodes_md5 = md5_for_file(nodes_file)
        else:
            raise TaxonomyException("Need to pass an md5 hash or a file.")

        self.ncbi_release = ncbi_release
        if downloaded_date is None:
            self.downloaded_date = datetime.datetime.utcnow().strftime(
                "%Y-%m-%d %H:%M")
        else:
            self.downloaded_date = downloaded_date

    def __repr__(self):
        return ("Taxonomy Metadata for names.dmp at %s and nodes.dmp at %s\n"
                "NCBI release: %s\n"
                "Downloaded: %s" %
                (self.names_uri, self.nodes_uri,
                 self.ncbi_release, self.downloaded_date))

    def to_dict(self):
        return {"names_uri": self.names_uri,
                "nodes_uri": self.nodes_uri,
                "names_md5": self.names_md5,
                "nodes_md5": self.nodes_md5,
                "ncbi_release": self.ncbi_release,
                "downloaded_date": self.downloaded_date
                }

    def to_json(self):
        return json.dumps(self.to_dict())


class Taxonomy(object):
    """
    networkx-based taxonomy object. Main object type of `taxonomy` module.

    Args:
        G (networkx.Graph or networkx.DiGraph): Input Taxonomy graph
        metadata (TaxonomyMetadata): Associated metadata

    Returns:
        Taxonomy object
    """
    def __init__(self, G, metadata):
        self.G = G
        self.metadata = metadata

    def __repr__(self):
        return ("Taxonomy object.\n"
                "Metadata: \n%s" % self.metadata)

    @staticmethod
    def load(f):
        """
        Load a Taxonomy from a .json file

        Args:
            f (str, file): A filepath or file handle pointing to
                           a taxonomy.json file

        Returns:
            Taxonomy
        """
        if isinstance(f, (file, gzip.GzipFile)):
            input_json = json.load(f)
        else:
            if os.path.splitext(f)[1] == ".gz":
                input_json = json.load(gzip.open(f, mode='r'))
            else:
                input_json = json.load(open(f, mode='r'))
        try:
            G = json_graph.node_link_graph(input_json["node_link_data"])
            metadata = TaxonomyMetadata(None, None, **input_json["metadata"])
        except KeyError:
            raise TaxonomyException("Improper input. Expects a JSON blob "
                                    "with node_link_graph and metadata "
                                    "parent nodes")
        return Taxonomy(G, metadata)

    @staticmethod
    def build_from_ncbi(names, nodes, names_uri, nodes_uri,
                        ncbi_release, downloaded_date=None, directed=True,
                        verbose=False):
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
            directed (bool): Build a directed graph? Defaults to True
            verbose (bool): Verbose build output. Defaults to False

        Returns:
            Taxonomy object built from NCBI data
        """
        # 0. Create metadata
        metadata = TaxonomyMetadata(names, nodes,
                                    names_uri, nodes_uri,
                                    ncbi_release,
                                    downloaded_date)

        # 1. Extract authoritative (scientific name type) names and map
        #    tax IDs -> scientific names
        names_dict = {}
        with open(names, mode='r') as names_file:
            _ = names_file.readline()  # Discard header
            for i, line in enumerate(names_file):
                line = line.strip().split("\t|\t")
                tax_id = int(line[0].strip())
                name_txt = line[1].strip()
                name_type = line[3].rstrip("\t|")
                if name_type == 'scientific name':
                    names_dict[tax_id] = name_txt

        # 2. Create a networkx object
        if directed:
            G = nx.DiGraph()
        else:
            G = nx.Graph()

        # 2. Create
        with open(nodes, mode='r') as nodes_file:
            for i, line in enumerate(nodes_file):
                line = line.split("\t|\t")
                tax_id = int(line[0].strip())
                parent_tax_id = int(line[1].strip())
                rank = line[2].strip()
                names_txt = names_dict.get(tax_id, "No name found.")
                hidden = line[10].strip()

                attr_dict = {
                    "name": names_txt,
                    "rank": rank,
                    "hidden": bool(hidden)
                }

                G.add_node(tax_id, attr_dict=attr_dict)
                G.add_edge(tax_id, parent_tax_id)
                if verbose and i % 1000 == 0:
                    print ("(Line %d) Adding %s (%s) (ID: %d. Parent: %d)" %
                           (i, names_txt, rank, tax_id, parent_tax_id))

        return Taxonomy(G, metadata)

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
        out["node_link_data"] = json_graph.node_link_data(self.G)
        out["metadata"] = self.metadata.to_dict()
        if isinstance(f, (file, gzip.GzipFile)):
            json.dump(out, f)
        else:
            if gzip:
                if os.path.splitext(f)[1] != ".gz":
                    f = f + ".gz"
                json.dump(out, gzip.open(f, mode='w'))
            else:
                json.dump(out, open(f, mode='w'))
