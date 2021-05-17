**pyPhase** is an open-source Python package for phase retrieval from phase contrast images in the Fresnel regime. For an overview, check out the pyPhase manuscript: https://arxiv.org/abs/2012.07942

Features
========

* Phase retrieval algorithms
* Wave propagation.
* Handling of different image sources and formats
* Tools for pre-processing such as registration of phase contrast images and motion estimation

Installation
============

Installation is currently through PyPI. Create a virtual environment using your favourite virtual environment manager, verify that pip is installed, then:

.. code-block::

   pip install pyphase


PyPhase currently requires Elastix 4.9 installed for registration. 
To manually install elastix 4.9 go to https://elastix.lumc.nl/download.php

* Unzip the archive.
* Add the path for elastix/bin to your .bashrc: add YOUR_PATH_TO_elastix/bin to your environment variable PATH.
* Add the path for elastix/lib to your .bashrc: add YOUR_PATH_TO_elastix/lib to your environment variable LD_LIBRARY_PATH.

Test your installation:

.. code-block::

    python3 import pyphase







