A standardized and extensible evaluation framework for taxonomic classifiers
============================================================================

To view static versions of the reports presented in (Bokulich, Rideout, et al., (In preparation)), [start here](http://nbviewer.ipython.org/github/gregcaporaso/short-read-tax-assignment/blob/master/ipynb/Index.ipynb).

This repository contains python code and IPython Notebooks designed to run on the QIIME 1.9.0 AWS instance. You can find
that AMI on the [QIIME resources page](http://qiime.org/home_static/dataFiles.html).

Setup and install
-----------------

The analyses presented in our paper were run on an m2.4xlarge instance (some of the assignment procedures can require large amounts of memory, so this allows for use of 4 parallel engines). 

A 100GB EBS volume was mounted as ``$HOME/data`` (generic instructions for this are [here](http://qiime.org/tutorials/working_with_aws.html#creating-a-volume-for-persistent-storage-across-different-launches-of-an-instance-or-different-instances), but any instructions for mounting an EBS volume should suffice and others may be more up-to-date. 

The ``$HOME/tmp`` directory was updated to be a symbolic link to ``/mnt/``, and ownership of ``/mnt`` was transferred to the ``ubuntu`` user. This was done as follows:

```
rm -r $HOME/temp/
ln -s /mnt/ $HOME/temp
sudo chown ubuntu /mnt
sudo chgrp ubuntu /mnt
```

The library code and IPython Notebooks are then installed as follows:

```
cd $HOME/data
git clone https://github.com/gregcaporaso/short-read-tax-assignment.git
cd $HOME/data/short-read-tax-assignment/code
sudo pip install .
```
(If you are not running this on the QIIME 1.9.0 AWS instance, you may need to run ``sudo pip install numpy`` before ``sudo pip install .``.)

To run the unit tests, you should install run:

```
cd $HOME/data/short-read-tax-assignment/code
nosetests .
```

Finally, download and unzip the reference databases:

```
cd $HOME/data/
wget https://dl.dropboxusercontent.com/u/2868868/unite-97-rep-set.tgz
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
tar -xzf unite-97-rep-set.tgz
tar -xzf gg_13_8_otus.tar.gz
```

Using the IPython Notebooks included in this repository
-------------------------------------------------------

To view and interact with an [IPython Notebook](http://ipython.org/notebook.html), change into the ``$HOME/data/short-read-tax-assignment/ipynbs`` directory and [start the IPython Notebook server in a screen session](http://qiime.org/tutorials/working_with_aws.html#connecting-to-your-qiime-ec2-instance-using-the-ipython-notebook).


**Everything below here needs to be updated after the re-analyses are completed.**

Evaluation workflows, for testing new taxonomic assignment methods
------------------------------------------------------------------

Several [IPython Notebooks](http://ipython.org/notebook.html) are provided to illustrate how to run a parameter sweep on a concept design for a new taxonomic assigner, and then evaluate that in the context of pre-computed results (the evaluation data from this study). The notebooks can be found under ``short-read-tax-assignment/ipynbs``.

The notebooks whose names begin with ``0`` illustrate how to generate results that can be analyzed using the evaluation framework presented here. These will likely have software requirements in addition to those  that are installed with ``short-read-tax-assignment``. Those will be listed in the top of each individual notebook.

The notebooks whose names begin with ``1`` are used to run the analyses performed in our study, and which you can re-run to include your data. These do not have requirements beyond those included in the ``short-read-tax-assignment`` installation, though you will need to update filepaths in one cell in that notebook to refer to locations on the system where you are executing the notebook.

Citing
------

If you use any of the data or code included in this repository, please cite with the URL: https://github.com/gregcaporaso/short-read-tax-assignment.
