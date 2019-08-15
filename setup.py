#!/usr/bin/env python3
import os
from setuptools import setup

try:
    from setuptools_rust import Binding, RustExtension
except ImportError:
    import subprocess

    errno = subprocess.call([
        '/usr/bin/env', 'python3', '-m', 'pip', 'install', 'setuptools-rust'
    ])
    if errno:
        print('Please install setuptools-rust package')
        raise SystemExit(errno)
    else:
        from setuptools_rust import Binding, RustExtension


with open('README.md', 'r') as fh:
    long_description = fh.read()


with open('Cargo.toml', 'r') as fh:
    # fast and dirty TOML parsing to get the version
    for line in fh:
        if '=' not in line:
            continue
        key, value = line.split('=')
        if key.strip() == 'version':
            version = value.strip().strip('"')
            break
    else:
        raise KeyError('version not found in Cargo.toml')


setup(
    name='taxonomy',
    version=version,
    author='One Codex',
    description='Load, save, and manipulate taxonomic trees',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/onecodex/taxonomy',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    rust_extensions=[
        RustExtension(
            'taxonomy',
            'Cargo.toml',
            features=['python'],
            debug=False,
            binding=Binding.PyO3
        ),
    ],
    packages=[],
    setup_requires=['setuptools-rust>=0.10.1', 'wheel'],
    install_requires=[],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False
)
