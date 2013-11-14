short-read-tax-assignment
=========================

A repository for storing code and data related to a systematic comparison of
short read taxonomy assignment tools. Includes a framework to facilitate the
evaluation and comparison of additional taxonomy assigners in the context of
the results presented in (Bokulich, Rideout, et al. (in preparation)).

Installing
----------

To run the code in this repository, you will need the development version of
[QIIME](http://www.qiime.org) (tested on 1.7.0-dev, master@dc1ca5e) and
[IPython](http://ipython.org/) >= 1.0.0 installed, with
```short-read-tax-assignment/code``` added to your ```PYTHONPATH``` and
```short-read-tax-assignment/code/scripts``` added to your ```PATH```.

Testing your installation
-------------------------

To test your installation/setup, you can run the unit tests with the following
commands (assuming you are in the ```short-read-tax-assignment``` directory and
you have [nose](http://nose.readthedocs.org/en/latest/) installed):

    nosetests code

Using the IPython Notebooks included in this repository
-------------------------------------------------------

To view and interact with an [IPython Notebook](http://ipython.org/notebook.html), change into the directory containing the IPython Notebooks of interest, and start an IPython Notebook server. For example:

```
    cd short-read-tax-assignment/demo/eval-demo
    ipython notebook
```

This will launch a web browser. You can then open the notebook by clicking on its name, and execute the notebook.

Evaluation workflows, for testing new taxonomic assignment methods
------------------------------------------------------------------

Two [IPython Notebooks](http://ipython.org/notebook.html) are provided to illustrate how to run a parameter sweep on a concept design for a new taxonomic assigner, and then evaluate that in the context of pre-computed results (the evaluation data from this study). The notebooks can be found under ```short-read-tax-assignment/demo/eval-demo```

Requirements for running the usearch example parameter sweep notebook, in addition to the code in this repository:

* Python 2.7
* [IPython Notebook](http://ipython.org/notebook.html)
* [usearch](http://www.drive5.com/usearch/) v5.2.236. 
* [pyqi](http://bipy.github.io/pyqi/doc/index.html) 0.2.0
* [biom-format](http://www.biom-format.org) 1.2.0
* [uc_to_assignments](https://gist.github.com/gregcaporaso/6083538) (A gist written for use in this notebook)

Requirements for running the evaluation framework notebook, in addition to the code in this repository:

* Python 2.7
* [IPython Notebook](http://ipython.org/notebook.html)
* [QIIME](http://qiime.org/)
* [PyCogent](https://github.com/pycogent/pycogent)
* [matplotlib](http://matplotlib.org/)
