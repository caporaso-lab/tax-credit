# Computational runtime comparison

The goal of these notebooks is to benchmark computational runtimes for taxonomy assignment methods.

* [Compute runtimes](./compute-runtimes.ipynb): Creates and executes commands for timing taxonomic assignment processes as a function of reference database size and query database size (where size is number of sequences in both cases). If you're adding new methods to this comparison, you should re-run all of the commands as the runtimes are only comparable if they are computed on the same system.
* [Analysis](./analysis.ipynb): Analysis of the runtimes of different methods.