================================
Absorbing Random-Walk Centrality
================================

|docs| |travis| |coveralls|
    
This is an implementation of the *absorbing random-walk centrality* measure for 
nodes in graphs. For the definition of the measure, as well as a study of the
related optimization problem and algorithmic techniques, please see the pre-print
publication on arXiv_. A short version of this paper will appear in the
`ICDM 2015`__.

.. _arXiv: http://arxiv.org/abs/1509.02533
__ http://icdm2015.stonybrook.edu/

To cite this work, please use

::

   Mavroforakis, Charalampos, Michael Mathioudakis, and Aristides Gionis.
   "Absorbing random-walk centrality: Theory and algorithms"
   Data Mining (ICDM), 2015 IEEE International Conference on. IEEE, 2015.


Installation
------------

You can install the *absorbing_centrality* package by executing the following command in a terminal.

::

   pip install git+https://github.com/harrymvr/absorbing-centrality#Egg=absorbing_centrality

Documentation
-------------

For instructions on how to use the package, consult `its documentation`__.

__ https://absorbing-centrality.readthedocs.org/

Development
-----------

To run all the tests for the code, you will need tox_ -- check its webpage for instructions on how to install it.

.. _tox: https://testrun.org/tox/latest/

Once tox_ is installed, use your terminal to enter the directory with the local copy of the code (here it's named '*absorbing-centrality*') and simply type the following command.

::

    absorbing-centrality $ tox

If everything goes well, you'll receive a congratulatory message. 


Note that the code is distributed under the Open Source Initiative (ISC) license.
For the exact terms of distribution, see the LICENSE_.

.. _LICENSE: ./LICENSE

::

   Copyright (c) 2015, absorbing-centrality contributors,
   Charalampos Mavroforakis <cmav@bu.edu>,
   Michael Mathioudakis <michael.mathioudakis@aalto.fi>,
   Aristides Gionis <aristides.gionis@aalto.fi>

    
.. |docs| image:: https://readthedocs.org/projects/absorbing-centrality/badge/?version=latest
    :target: https://absorbing-centrality.readthedocs.org/en/latest/
    :alt: Documentation Statu

.. |travis| image:: https://travis-ci.org/harrymvr/absorbing-centrality.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/harrymvr/absorbing-centrality

.. |requires| image:: https://requires.io/github/harrymvr/absorbing-centrality/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/harrymvr/absorbing-centrality/requirements/?branch=master


.. |coveralls| image:: https://coveralls.io/repos/harrymvr/absorbing-centrality/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/github/harrymvr/absorbing-centrality?branch=master


.. |version| image:: https://img.shields.io/pypi/v/absorbing_centrality.svg?style=flat
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/absorbing_centrality

.. |downloads| image:: https://img.shields.io/pypi/dm/absorbing_centrality.svg?style=flat
    :alt: PyPI Package monthly downloads
    :target: https://pypi.python.org/pypi/absorbing_centrality

.. |wheel| image:: https://img.shields.io/pypi/wheel/absorbing_centrality.svg?style=flat
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/absorbing_centrality

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/absorbing_centrality.svg?style=flat
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/absorbing_centrality

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/absorbing_centrality.svg?style=flat
    :alt: Supported imlementations
    :target: https://pypi.python.org/pypi/absorbing_centrality

