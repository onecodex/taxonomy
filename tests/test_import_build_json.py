import os

from taxonomy.formats import taxonomy_from_build_json


def test_import_build_json():
    file_name = os.path.join(os.getcwd(), 'tests/build.json.gz')
    tax = taxonomy_from_build_json(file_name)

    assert tax.parents('9606') == ['3', '1']
    assert set(tax.children('1')) == {'2', '3'}
