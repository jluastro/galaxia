galaxia
-------

`galaxia <http://galaxia.sourceforge.net>`_ is a code for generating a
synthetic model of the galaxy. This repository is a modified version of
galaxia that enables the user to modify the parameters of the galaxy model.
This repository outlines a slightly different installation than the
one described on the galaxia documentation page.


Installation
------------

To install galaxia, first clone the GitHub repository.

.. code-block:: bash

    git clone git@github.com:jluastro/galaxia.git

Next download and extract the GalaxiaData folder
from our online server at astro.berkeley.edu/path/to/GalaxiaData_ .

.. code-block:: bash

    wget astro.berkeley.edu/path/to/GalaxiaData.tar.gz
    tar xvf GalaxiaData.tar.gz

This folder will need to be unique to each galaxy model that you
run on galaxia, so you may need to copy multiple versions for different
galaxy models.

Next compile and make the galaxia installation. We recommend following the
non root installation outlined in the galaxia documentation. The prefix argument
specifies the folder where the executable will be written to. As stated in the
original documentation, if you want to install in home/user/sw/bin/ the run:

.. code-block:: bash

    ./configure --prefix=home/user/sw
    make
    make install

Note that our version does not require specifying a data directory. Doing so
has no effect on the execution of the code. Instead the user is asked to
specify the location of the GalaxiaData folder in a galaxyModel parameter file.

You will now have a galaxia executable located on your computer.

Running galaxia
----------------

How to run galaxia is outlined on the `galaxia documentation page <http://galaxia.sourceforge.net/Galaxia3pub.html>`_.
We provide an example galaxia parameter file located
at `example_galaxiaParams.txt <example_galaxiaParams.txt>`_.

Our version requires an additional parameter file that states
the location of the GalaxiaData directory and the galaxy model parameters.
An example galaxyModel parameter file is located
at `example_galaxyModelParams.txt <example_galaxyModelParams.txt>`_.

To run galaxia with this parameter file, place it as the next argument after the
regular galaxia parameter file.

.. code-block:: bash

    galaxia -r example_galaxiaParams.txt example_galaxyModelParams.txt

Make sure that the GalaxiaData directory specified in your galaxyModel parameter file
points to a unique directory for each different set of galaxy model parameters.

Our version will automatically generate the BHTree files required to run galaxia
if they are not present in the GalaxiaData directory. However you can choose to
force the calculation of the BHTree files by running:

.. code-block:: bash

    galaxia -s [warp or nowarp] example_galaxyModelParams.txt
