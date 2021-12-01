---
layout: default
title: clustering
nav_order: 10
parent: Modules
---

10\. `clustering`: This module performs density-based, hierarchical clustering (HDBSCAN) on tSNE or UMAP embeddings of post-QC data using an approach similar to that described in the [`metaQC`]({{ site.baseurl }}/modules/metaQC). Clicking the "save" button after an optimal `min_cluster_size` has been selected causes the program to append the current cluster IDs to the post-QC dataframe and proceed to the next module.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `embeddingAlgorithm` | UMAP (str) | Embedding algorithm to use for clustering (options: TSNE or UMAP) |
| `channelExclusionsClustering` | [ ] (list) | Immunomarkers to exclude from clustering and all subsequent modules (strs) |
| `samplesToRemoveClustering` | [ ] (list) | Samples to exclude from clustering and all subsequent modules (strs). |
| `fracForEmbedding` | 1.0 (float) | Fraction of cells to be embedded (range: 0.0-1.0) limits the amount of data passed to downstream modules |
| `dimensionEmbedding` | 2 (int) | Dimension of the embedding (int, typically 2) |
| `metricQC` | euclidean (str) | Distance metric for computing embedding. Choose from valid metrics used by scipy.spatial.distance.pdist: ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’, ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘jensenshannon’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’. |
| `perplexity` | 50.0 (float) | This is a [tSNE-specific configuration](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html) related to the number of nearest neighbors used in other manifold learning algorithms. Larger datasets usually require larger perplexity. Different values can result in significantly different results. |
| `earlyExaggeration` | 12.0 (float) | This is a [tSNE-specific configuration](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html). For larger values, the space between natural clusters will be larger in the embedded space. |
| `learningRateTSNE` | 200.0 (float) | This is a [tSNE-specific configuration](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html). tSNE learning rate (typically between 10.0 and 1000.0) |
| `randomState` | 5 (int) | This is a [tSNE-specific configuration](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html). It determines the random number generator for reproducible results across multiple function calls. |
| `nNeighbors` | 5 (int) | This is a [UMAP-specific configuration](https://umap-learn.readthedocs.io/en/latest/api.html). It determines the size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100. |
| `learningRateUMAP` | 1.0 (float) | This is a [UMAP-specific configuration](https://umap-learn.readthedocs.io/en/latest/api.html). It Determines the initial learning rate for the embedding optimization. |
| `minDist` | 0.1 (float) | This is a [UMAP-specific configuration](https://umap-learn.readthedocs.io/en/latest/api.html). Determines the effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together, while larger values will result on a more even dispersal of points. The value should be set relative to the spread value, which determines the scale at which embedded points will be spread out (float). |
| `repulsionStrength` | 5.0 (float) | This is a [UMAP-specific configuration](https://umap-learn.readthedocs.io/en/latest/api.html). Determines the weighting applied to negative samples in low dimensional embedding optimization. Values higher than one will result in greater weight being given to negative samples. |
