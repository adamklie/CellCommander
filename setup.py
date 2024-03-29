#!/usr/bin/env python

import os
import setuptools
import codecs
from typing import List


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version() -> str:
    for line in read('cellcommander/__init__.py').splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


def readme() -> str:
    with open('README.md') as f:
        return f.read()


def _readlines(filename, filebase=''):
    with open(os.path.join(filebase, filename)) as f:
        lines = f.readlines()
    return lines


def get_requirements() -> List[str]:
    requirements = _readlines('requirements.txt')
    if 'READTHEDOCS' in os.environ:
        requirements.extend(get_rtd_requirements())
    return requirements


def get_rtd_requirements() -> List[str]:
    requirements = _readlines('requirements-rtd.txt')
    return requirements


setuptools.setup(
    name='cellcommander',
    version=get_version(),
    description='Common single cell analysis tasks made runnable with the command line',
    long_description=readme(),
    long_description_content_type='text/markdown',
    classifiers=[
      'Development Status :: 4 - Beta',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: BSD License',
      'Programming Language :: Python :: 3.7',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='single-cell bioinformatics',
    url='https://github.com/adamklie/CellCommander',
    author='Adam Klie',
    license='BSD (3-Clause)',
    packages=setuptools.find_packages(),
    install_requires=get_requirements(),
    extras_require={
        "dev": ["pytest", "scikit-learn"],
        "docs": get_rtd_requirements(),
    },
    entry_points={
        'console_scripts': ['cellcommander=cellcommander.base_cli:main'],
    },
    include_package_data=True,
    package_data={'': ['*.ipynb']},  # include the report template
    zip_safe=False
)
