# Mock community evaluations

This notebook describes how to apply the mock community evaluations presented in (Bokulich, Rideout, et al. (in preparation)) to reproduce the analyses in that paper, or to extend them to other data sets. 

## Structuring new results for comparison to precomputed results
To prepare results from another classifier for analysis, you'll need to have [BIOM](http://www.biom-format.org) files with taxonomy assignments as an observation metadata category called ``taxonomy``. An example of how to generate these is presented in the [data generation notebook](./generate-tax-assignments.ipynb) in this directory, which was used to generated the precomputed data in the [tax-credit repository](https://github.com/caporaso-lab/short-read-tax-assignment/).

Your BIOM tables should be called ``table.biom``, and nested in the following directory structure:

```
results_dir/
 mock-community/
  dataset-id/ 
   reference-db-id/
    method-id/
     parameter-combination-id/
      table.biom
```

``results_dir`` is the name of the top level directory, and you will set this value in the first code cell of the analysis notebooks. You can name this directory whatever you want to. ``mock-community`` describes the specific analysis that is being run, and must be named ``mock-community`` for the framework to find your results.

This directory structure is identical to that for the [precomputed results](https://github.com/caporaso-lab/short-read-tax-assignment/tree/master/data/precomputed-results). You can review that directory structure for an example of how this should look.

## Contents
* [Data generation](./generate-tax-assignments.ipynb): Creates and executes commands for generating taxonomic assignments for the mock community contained in this package. The results of running this notebook are included in the repository, so it's not necessary to re-run this.
* [Analysis base](./base.ipynb): Template for mock community analysis at different taxonomic levels (all of the notebooks below are auto-generated from this notebook).
 * [Pre-computed phylum report](./phylum.ipynb)
 * [Pre-computed class report](./class.ipynb)
 * [Pre-computed order report](./order.ipynb)
 * [Pre-computed family report](./family.ipynb)
 * [Pre-computed genus report](./genus.ipynb)
 * [Pre-computed species report](./species.ipynb) 