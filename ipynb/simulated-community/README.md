# Simulated Communities

Simulated communities are synthetically composed collections of sequences, which are intended to resemble sequence collections assembled from natural microbial communities. Simulated communities share some similarities to the simulated reads used for [cross-validation of taxonomy assignment](../cross-validated/index.ipynb), in that they are assembled from reference sequences with known taxonomy annotations. Simulated communities have the added exceptions that these artifical sequences are present at known abundances, and simulate a natural community. While the composition of a simulated community is derived from the observed composition of a natural community, there is an essential difference: as the sequences in a simulated community are drawn from annotated reference sequences, the taxonomic composition of a simulated community is known in advance, whereas a natural community is not (and even after taxonomic classification, is still not known with certainty). This makes simulated communities useful for assessing the performance of taxonomy classification and other methods. Given the similarities with simulated reads used for cross-validated taxonomy assignment, method performance should perform fairly similar, but simulated communities give us some idea of how method and parameter configuration affect the observed composition of microbial communities.

Simulated communities have an additional function: to test prior probability training for machine-learning classifiers. Prior probabilities are used by some classification methods to predict the likelihood that a given class label will be observed in a population. Training prior probabilities may improve classification accuracy, particularly in situations when specific taxa are unique to an environment. This is opposed to the use of a uniform prior, which assumes that any class is equally likely to be observed in a given environment. While this is a prudent approach when handling unknown samples, it is a naive approach when handling samples from previously characterized sample types.

## Contents
The following notebooks generate, assign taxonomy to, and evaluate simulated communities:

1) **[Generate simulated communities](./dataset-generation.ipynb)** derived from natural community compositions

2) **[Assign taxonomy](./taxonomy-assignment.ipynb)** to simulated community sequences.

3) **[Evaluate classification accuracy](./evaluate-classification-accuracy.ipynb)** of simulated communities.
