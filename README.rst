galaxia
-------

`galaxia <http://galaxia.sourceforge.net>`_ is a code for generating a
synthetic model of the galaxy. This repository is a modified version of
galaxia that enables the user to modify the parameters of the galaxy model.
This README outlines a slightly different installation and usage than the
one described on the galaxia documentation page.


Installation
------------

To install galaxia, first clone the GitHub repository.

.. code-block:: bash

    git clone git@github.com:jluastro/galaxia.git

Next download and extract the GalaxiaData folder
from our online server at :code:`http://bhs.astro.berkeley.edu/GalaxiaData.tar.gz`.

.. code-block:: bash

    wget http://bhs.astro.berkeley.edu/GalaxiaData.tar.gz
    tar xvf GalaxiaData.tar.gz

This folder will need to be unique to each galaxy model that you
run on galaxia, so you may need to download multiple versions for
different galaxy models parameters.

Next compile, make, and install the galaxia installation. We recommend following the
non root installation outlined in
the `galaxia documentation <http://galaxia.sourceforge.net/Galaxia3pub.html>`_.
The prefix argument specifies the folder where the executable will be written to.
As stated in the original documentation, if you want to install in :code:`home/user/sw/bin/` then run:

.. code-block:: bash

    ./configure --prefix=home/user/sw
    make
    make install

Note that our version does not require specifying a data directory. Doing so
has no effect on the execution of the code. Instead the user is asked to
specify the location of the GalaxiaData folder in a galaxyModel parameter file
as described below.

You will now have a galaxia executable located in :code:`home/user/sw/bin/`. If this
directory does not already exist in your PATH, then make sure to prepend it.

.. code-block:: bash

    export PATH=home/user/sw/bin/:$PATH

Running galaxia
----------------

The functions and features of galaxia are outlined on the
`galaxia documentation page <http://galaxia.sourceforge.net/Galaxia3pub.html>`_.
We provide an example of the required galaxia parameter file
at `example_galaxiaParams.txt <docs/example_galaxiaParams.txt>`_.

Our version requires an additional parameter file that states
the location of the GalaxiaData directory and the galaxy model parameters.
An example galaxyModel parameter file is located
at `example_galaxyModelParams.txt <docs/example_galaxyModelParams.txt>`_.

To run galaxia with this parameter file, place it as the next argument after the
regular galaxia parameter file.

.. code-block:: bash

    galaxia -r example_galaxiaParams.txt example_galaxyModelParams.txt

Make sure that the GalaxiaData directory specified in your galaxyModel parameter file
points to a unique directory for each different set of galaxy model parameters.

Our version will automatically generate the BHTree files required to run galaxia
if they are not present in the GalaxiaData directory. This will occur the first time
you run a galaxia simulation. However you can choose to
force the calculation of the BHTree files by running:

.. code-block:: bash

    galaxia -s [warp or nowarp] example_galaxyModelParams.txt

BHTree files will now exist in the GalaxyData folder specified
in :code:`example_galaxyModelParams.txt`.