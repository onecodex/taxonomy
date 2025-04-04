#!/usr/bin/env python3

import taxonomy
from downloads import download
import subprocess
import unittest


class LatestNCBITestCase(unittest.TestCase):
    def test_load_latest_ncbi_taxonomy(self):
        download("https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
        subprocess.check_output(["tar", "-zxvf", "taxdump.tar.gz"])
        taxonomy.Taxonomy.from_ncbi(".")


if __name__ == "__main__":
    unittest.main()
