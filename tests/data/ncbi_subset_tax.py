""" Script used to generate ncbi_subset_tax.json

NCBI taxonomy downloaded in 2019 July 1st.
"""

from taxonomy import Taxonomy

tax = Taxonomy.from_ncbi("nodes.dmp", "names.dmp")
new_tax = tax.prune(remove=['28384', '12908', '10239', '2', '2157'] + [t for t in tax.children("2759") if t != "543769"])

with open("new_tax.json", "wb") as f:
    f.write(new_tax.to_json(True))
