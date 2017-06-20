# TAX CREdiT: TAXonomic ClassifieR Evaluation Tool

### A standardized and extensible evaluation framework for taxonomic classifiers

To view static versions of the reports presented in [Bokulich, et al., (Microbiome, under review)](https://peerj.com/preprints/934/), [start here](http://nbviewer.jupyter.org/github/nbokulich/short-read-tax-assignment/blob/master/ipynb/Index.ipynb).


Environment
-----------------
This repository contains python-3 code and Jupyter notebooks, but some taxonomy assignment methods (e.g., using QIIME-1 legacy methods) may require different python or software versions. Hence, we use conda parallel environments to support comparison of myriad methods in a single framework.

The first step is to create a conda environment with the necessary dependencies. This requires installing [miniconda 3](http://conda.pydata.org/miniconda.html) to manage parallel python environments. After miniconda (or another conda version) is installed, proceed with [installing QIIME 2](https://docs.qiime2.org/2.0.6/install/).

An example of how to load different environments to support other methods can be see in the [QIIME-1 taxonomy assignment notebook](https://github.com/nbokulich/short-read-tax-assignment/tree/master/ipynb/mock-community/generate-tax-assignments.ipynb).


Setup and install
-----------------
The library code and IPython Notebooks are then installed as follows:

```
cd $HOME/projects
git clone https://github.com/gregcaporaso/short-read-tax-assignment.git
cd $HOME/projects/short-read-tax-assignment/code
sudo pip install .
```

To run the unit tests, you should install run:

```
cd $HOME/projects/short-read-tax-assignment/code
nosetests .
```

Finally, download and unzip the reference databases:

```
cd $HOME/ref_dbs/
wget https://unite.ut.ee/sh_files/sh_qiime_release_20.11.2016.zip
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
unzip sh_qiime_release_20.11.2016.zip
tar -xzf gg_13_8_otus.tar.gz
```

Equipment
------------------
The analyses included here can all be run in standard, modern laptop, provided you don't mind waiting a few hours on the most memory-intensive step (taxonomy classification of millions of sequences). All analyses presented in ``tax-credit`` were run in a single afternoon using a MacBook Pro with the following specifications:
**OS** OS X 10.11.6 "El Capitan"
**Processor** 2.3 GHz Intel Core i7
**Memory** 8 GB 1600 MHz DDR3


Using the Jupyter Notebooks included in this repository
-------------------------------------------------------

To view and interact with [Jupyter Notebook](http://jupyter.org/), change into the ``/short-read-tax-assignment/ipynb`` directory and run Jupyter Notebooks from the terminal with the command:

``jupyter notebook index.ipynb``

The notebooks menu should open in your browser. From the main index, you can follow the menus to browse different analyses, or use ``File --> Open`` from the notebook toolbar to access the full file tree.


Citing
------

If you use any of the data or code included in this repository, please cite with:

Bokulich NA, Rideout JR, Kopylova E, Bolyen E, Patnode J, Ellett Z, McDonald D, Wolfe B, Maurice CF, Dutton RJ, Turnbaugh PJ, Knight R, Caporaso JG. (2015) A standardized, extensible framework for optimizing classification improves marker-gene taxonomic assignments. PeerJ PrePrints 3:e1156 https://dx.doi.org/10.7287/peerj.preprints.934v1
