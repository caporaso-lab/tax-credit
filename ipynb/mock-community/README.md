# Mock community evaluations

This notebook describes how to apply the mock community evaluations presented in (Bokulich, Kaehler, et al. (in preparation)) to reproduce the analyses in that paper, or to extend them to other data sets. 

## Structuring new results for comparison to precomputed results
To prepare results from another classifier for analysis, you'll need to have [BIOM](http://www.biom-format.org) files with taxonomy assignments as an observation metadata category called ``taxonomy``. An example of how to generate these is presented in the [data generation notebook](./mock-dataset-generation.ipynb) in this directory, which was used to generated the precomputed data in the [tax-credit repository](https://github.com/caporaso-lab/tax-credit/).

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

This directory structure is identical to that for the [precomputed results](https://github.com/caporaso-lab/tax-credit/tree/master/data/precomputed-results). You can review that directory structure for an example of how this should look.

## Contents
* [Data Generation](./dataset-generation.ipynb): retrieve and process mock community data sets for analysis. (These data are stored in tax-credit and will not need to be performed again, unless if preparing new mock community data sets for analysis.)
* [Taxonomy Assignment](./taxonomy-assignment-template.ipynb): Creates and executes commands for generating taxonomic assignments for the mock community contained in this package. The results of running this notebook are included in the repository, so it's not necessary to re-run this. Start [here](./taxonomy-assignment-template.ipynb) for testing new methods. Examples are provided for the following methods:
    * [Qiime 2](./taxonomy-assignment-q2-feature-classifier.ipynb): ``q2-feature-classifier``.
    * [Qiime 1](./taxonomy-assignment-qiime1.ipynb): RDP, sortmerna, uclust, and BLAST wrappers.
    * [Reference database comparison](./evaluate-classification-database-comparison.ipynb): and example of taxonomy assignment using different reference databases, this is the taxonomy generation step preceding the [comparative evaluation](./evaluate-classification-database-comparison.ipynb).
* [Analysis](./evaluate-classification-accuracy.ipynb): Template for mock community analysis at multiple taxonomic levels.
* Comparison of taxonomic assignment with different reference databases, for example, can be performed through modifications of the above notebooks to [format reference databases](./format-reference-databases.ipynb) to different specifications, [assign taxonomy](./generate-tax-assignments-trimmed-dbs.ipynb) with these databases, and [compare classification](./evaluate-classification-database-comparison.ipynb) of mock communities assigned with different reference databases.