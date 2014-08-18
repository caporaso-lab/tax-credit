short-read-tax-assignment
=========================

A repository for storing code and data related to a systematic comparison of
short read taxonomy assignment tools. Includes a framework to facilitate the
evaluation and comparison of additional taxonomy assigners in the context of
the results presented in (Bokulich, Rideout, et al. (Nature Methods, under review)).

Install
-------

To run the code in this repository, run the following commands:

    git clone https://github.com/gregcaporaso/short-read-tax-assignment.git
    cd short-read-tax-assignment/code
    pip install numpy
    pip install .

Testing
-------

To run the unit tests, you should install [nose](http://nose.readthedocs.org/en/latest/). Then change to the ```short-read-tax-assignment``` directory and run:

    nosetests code

Using the IPython Notebooks included in this repository
-------------------------------------------------------

To view and interact with an [IPython Notebook](http://ipython.org/notebook.html), change into the ``short-read-tax-assignment`` directory and launch the IPython Notebook server:

    ipython notebook

This will launch the IPython Notebook interface in a new web browser window. You can then navigate to the ``ipynbs`` directory, and open the notebook of interest.

Evaluation workflows, for testing new taxonomic assignment methods
------------------------------------------------------------------

Several [IPython Notebooks](http://ipython.org/notebook.html) are provided to illustrate how to run a parameter sweep on a concept design for a new taxonomic assigner, and then evaluate that in the context of pre-computed results (the evaluation data from this study). The notebooks can be found under ``short-read-tax-assignment/ipynbs``.

The notebooks whose names begin with ``0`` illustrate how to generate results that can be analyzed using the evaluation framework presented here. These will likely have software requirements in addition to those  that are installed with ``short-read-tax-assignment``. Those will be listed in the top of each individual notebook.

The notebooks whose names begin with ``1`` are used to run the analyses performed in our study, and which you can re-run to include your data. These do not have requirements beyond those included in the ``short-read-tax-assignment`` installation, though you will need to update filepaths in one cell in that notebook to refer to locations on the system where you are executing the notebook.
