---
layout: default-cylinter
title: clustering
nav_order: 10
parent: Modules
---

10\. `clustering`: This module performs density-based, hierarchical clustering with HDBSCAN on tSNE or UMAP embeddings of cells using an approach similar to that described in the [`metaQC`]({{ site.baseurl }}/modules/metaQC) module. Clicking the "save" button after arriving at an optimal clustering causes the program to append the current cluster IDs to the dataframe before proceeding to the next module.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `embeddingAlgorithm` | "UMAP" | (str) Embedding algorithm to use for clustering (options: "TSNE" or "UMAP"). |
| `channelExclusionsClustering` | [ ] | (list of strs) Immunomarkers to exclude from clustering. |
| `samplesToRemoveClustering` | [ ] | (list of strs) Samples to exclude from clustering. |
| `normalizeTissueCounts` | True | (bool) Make the number of cells per tissue for clustering more similar through sample-weighted random sampling. |
| `fracForEmbedding` | 1.0 | (float) Fraction of cells to be embedded (range: 0.0-1.0). Limits amount of data passed to downstream modules. |
| `dimensionEmbedding` | 2 | (int) Dimension of the embedding (options: 2 or 3). |
| `topMarkers` | "channels" | (str) Normalization axis ("channels" or "clusters") used to define highest expressed markers per cluster. |
| `colormapChannel` | null | (null or str) Channel to colormap to the embedding. |
| `metricQC` | "euclidean" | (str) Distance metric for computing embedding. Choose from valid metrics used by scipy.spatial.distance.pdist: "braycurtis", "canberra", "chebyshev", "cityblock", "correlation", "cosine", "dice", "euclidean", "hamming", "jaccard", "jensenshannon", "kulsinski", "mahalanobis", "matching", "minkowski", "rogerstanimoto", "russellrao", "seuclidean", "sokalmichener", "sokalsneath", "sqeuclidean", "yule". |
| `perplexity` | 50.0 | (float) This is a [tSNE-specific configuration](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html) related to the number of nearest neighbors used in other manifold learning algorithms. Larger datasets usually require larger perplexity. Different values can result in significantly different results. |
| `earlyExaggeration` | 12.0 | (float) This is a [tSNE-specific configuration](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html). For larger values, the space between natural clusters will be larger in the embedded space. |
| `learningRateTSNE` | 200.0 | (float) This is a [tSNE-specific configuration](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html). tSNE learning rate (typically between 10.0 and 1000.0). |
| `randomState` | 5 | (int) This is a [tSNE-specific configuration](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html). It determines the random number generator for reproducible results across multiple function calls. |
| `nNeighbors` | 5 | (int) This is a [UMAP-specific configuration](https://umap-learn.readthedocs.io/en/latest/api.html). It determines the size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100. |
| `learningRateUMAP` | 1.0 | (float) This is a [UMAP-specific configuration](https://umap-learn.readthedocs.io/en/latest/api.html). It Determines the initial learning rate for the embedding optimization. |
| `minDist` | 0.1 | (float) This is a [UMAP-specific configuration](https://umap-learn.readthedocs.io/en/latest/api.html). Determines the effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together, while larger values will result on a more even dispersal of points. The value should be set relative to the spread value, which determines the scale at which embedded points will be spread out. |
| `repulsionStrength` | 5.0 | (float) This is a [UMAP-specific configuration](https://umap-learn.readthedocs.io/en/latest/api.html). Determines the weighting applied to negative samples in low dimensional embedding optimization. Values higher than one will result in greater weight being given to negative samples. |
