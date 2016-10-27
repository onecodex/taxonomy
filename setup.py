#!/usr/bin/env python
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
    version='0.2.0',
    url='http://github.com/onecodex/taxonomy/',
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
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)
