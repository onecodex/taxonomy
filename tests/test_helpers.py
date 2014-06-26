from nose.tools import *
import os
from taxonomy.helpers import md5_for_file


def test_md5_for_file():
    init_md5 = md5_for_file(os.path.join(os.getcwd(), "tests/__init__.py"))
    assert ("d41d8cd98f00b204e9800998ecf8427e" == init_md5)
