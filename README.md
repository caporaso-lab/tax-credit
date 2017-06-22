# TAX CREdiT: TAXonomic ClassifieR Evaluation Tool

[![Build Status](https://travis-ci.org/caporaso-lab/tax-credit.svg?branch=master)](https://travis-ci.org/caporaso-lab/tax-credit)

### A standardized and extensible evaluation framework for taxonomic classifiers

To view static versions of the reports , [start here](https://github.com/caporaso-lab/tax-credit/blob/master/ipynb/Index.ipynb).


Environment
-----------------
This repository contains python-3 code and Jupyter notebooks, but some taxonomy assignment methods (e.g., using QIIME-1 legacy methods) may require different python or software versions. Hence, we use conda parallel environments to support comparison of myriad methods in a single framework.

The first step is to install conda and install QIIME2 following the instructions provided [here](https://docs.qiime2.org/2017.6/install/native/).

An example of how to load different environments to support other methods can be see in the [QIIME-1 taxonomy assignment notebook](https://github.com/caporaso-lab/tax-credit/blob/master/ipynb/mock-community/taxonomy-assignment-qiime1.ipynb).


Setup and install
-----------------
The library code and IPython Notebooks are then installed as follows:

```
git clone https://github.com/gregcaporaso/tax-credit.git
cd tax-credit/
pip install .
```

Finally, download and unzip the reference databases:

```
wget https://unite.ut.ee/sh_files/sh_qiime_release_20.11.2016.zip
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
unzip sh_qiime_release_20.11.2016.zip
tar -xzf gg_13_8_otus.tar.gz
```

Equipment
------------------
The analyses included here can all be run in standard, modern laptop, provided you don't mind waiting a few hours on the most memory-intensive step (taxonomy classification of millions of sequences). With the exception of the `q2-feature-classifier naive-bayes*` classifier sweeps, which were run on a high-performance cluster, all analyses presented in ``tax-credit`` were run in a single day using a MacBook Pro with the following specifications:
**OS** OS X 10.11.6 "El Capitan"
**Processor** 2.3 GHz Intel Core i7
**Memory** 8 GB 1600 MHz DDR3

If you intend to perform extensive parameter sweeps on a classifier (e.g., several hundred or more parameter combinations), you may want to consider running these analyses using cluster resources, if available.


Using the Jupyter Notebooks included in this repository
-------------------------------------------------------

To view and interact with [Jupyter Notebook](http://jupyter.org/), change into the ``/tax-credit/ipynb`` directory and run Jupyter Notebooks from the terminal with the command:

``jupyter notebook index.ipynb``

The notebooks menu should open in your browser. From the main index, you can follow the menus to browse different analyses, or use ``File --> Open`` from the notebook toolbar to access the full file tree.


Citing
------

A publication is on its way! For now, if you use any of the data or code included in this repository, please cite https://github.com/caporaso-lab/tax-credit
