#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import, print_function

import io
import os
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import relpath
from os.path import splitext

from setuptools import find_packages
from setuptools import setup

def read(*names, **kwargs):
    return io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ).read()

tests_require=['nose']

setup(
    name='absorbing_centrality',
    version='0.1.0',
    license='ISC',
    description='An implementation of the absorbing random-walk centrality measure for graphs.',
    long_description='%s\n%s' % (read('README.rst'), re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst'))),
    author='Charalampos Mavroforakis',
    author_email='cmav@bu.edu',
    url='https://github.com/harrymvr/absorbing-centrality',
    packages=['absorbing_centrality'],
    # package_dir={'': 'absorbing_centrality'},
    # py_modules=[splitext(basename(path))[0] for path in glob('absorbing_centrality/*.py')],
    include_package_data=True,
    zip_safe=False,
    test_suite='nose.collector',
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Natural Language :: English'
    ],
    keywords=['graph mining', 'node centrality', 'random walks', 'algorithms',
              'data mining'
    ],
    install_requires=[
        'networkx>=1.9.1',
        'numpy==1.9.2',
        'scipy==0.16'
    ],
    extras_require={
        'tests': tests_require,
    },
)
