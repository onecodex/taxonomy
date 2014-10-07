"""
``taxonomy``
---------

``taxonomy`` provides a set of convenient graph-based objects and functions
for understanding taxonomic/phylogenetic relationships.



Links
`````
* `Code on Github: <http://github.com/refgenomics/taxonomy/>`

"""
from setuptools import setup


setup(
    name='taxonomy',
    version='0.1.0',
    url='http://github.com/refgenomics/taxonomy/',
    license='All rights reserved',
    author='Nick Boyd Greenfield',
    author_email='nick@refgenomics.com',
    description='Graph-based taxonomy objects for genomics',
    long_description=__doc__,
    packages=['taxonomy'],
    zip_safe=True,
    platforms='any',
    install_requires=[
        'networkx>=1.9',
        'pandas>=0.14',
        'nose>=1.3.3'
    ],
    test_suite='nose.collector'
)
