short-read-tax-assignment
=========================

A repository for storing code and data related to a systematic comparison of
short read taxonomy assignment tools.

To run the code in this repository, you will need the development version of
[QIIME](http://www.qiime.org) (tested on 1.7.0-dev, master@7b0920b) and
[IPython](http://ipython.org/) 1.0.dev installed, with
```short-read-tax-assignment/code``` added to your ```PYTHONPATH``` and
```short-read-tax-assignment/code/scripts``` added to your ```PATH```.

To test your installation/setup, you can run the unit tests with the following
commands (assuming you are in the ```short-read-tax-assignment``` directory):

    python code/tests/test_multiple_assign_taxonomy.py
    python code/tests/test_generate_taxa_compare_table.py

Two [IPython Notebooks](http://ipython.org/notebook.html) are provided to
demonstrate the workflows that were used to compare taxonomy assigners using
mock community data and natural community data. The notebooks and demonstration
material can be found under ```short-read-tax-assignment/demo/mock-demo``` and
```short-read-tax-assignment/demo/natural-demo```. To view and interact with a
notebook, change into the demo's directory and start an IPython Notebook
server. For example:

    cd short-read-tax-assignment/demo/mock-demo
    ipython notebook

This will launch a web browser. You can then open the notebook by clicking on
its name and run the demo workflow.
