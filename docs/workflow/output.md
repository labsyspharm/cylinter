---
layout: default
title: Output
nav_order: 2
parent: Workflow
---

# Output File Structure


In the below example, `<sample-name>` corresponds to the names given to the various tissue samples in a given analysis. `<test-name>` refers to the binary declarations specified as the fourth value elements in the `sampleMetadata` dictionary (see [general configurations](/input#yaml-configuration-file)). `<chunk>` refers to a slice of the combined single-cell feature table for QC status reclassification, and `<channel>` refers to the various immunomarkers used in the study.

``` bash
<output_dir>
├── area/
│   ├── <sample-name>.pdf
│   └── idxs_to_drop.csv
├── checkpoints/
│   ├── aggregateData.parquet
│   ├── areaFilter.parquet
│   ├── clustering.parquet
│   ├── clustermap.parquet
│   ├── curateThumbnails.parquet
│   ├── cycleCorrelation.parquet
│   ├── frequencyStats.parquet
│   ├── intensityFilter.parquet
│   ├── logTransform.parquet
│   ├── metaQC.parquet
│   ├── PCA.parquet
│   ├── pruneOutliers.parquet
│   ├── selectROIs.parquet
│   └── setContrast.parquet
├── clustering/
│   ├── clustermap_norm_channels.pdf
│   ├── clustermap_norm_clusters
│   ├── embedding.npy
│   frequency_stats/
│   ├── <test-name>/
│       ├── catplot.pdf
│       ├── plot.pdf
│       ├── stats_sig.csv
│       └── stats_total.csv
├── contrast/
│   └── contrast_limits.yml
├── cycles/
│   ├── cycle_correlation(logRatio).pdf
│   ├── cycle_correlation(perCycle).pdf
│   ├── cycle_correlation(perSample).png
│   └── idxs_to_drop.csv
├── intensity/
│   ├── <sample-name>.pdf
│   └── idxs_to_drop.csv
├── metaQC/
│   ├── <chunk>/
│   ├── chunk_index.txt
│   ├── QCData.pkl
│   └── reclass_storage_dict.pkl
├── PCA/
│   └── pcaScoresPlot.pdf
├── pruning/
│   ├── <channel>_pruned_rescaled.png
│   └── <channel>_raw.png
│   ├── data_copy1.parquet
│   └── pruning_dict.csv
└── ROIs/
    ├── idxs_to_drop.csv
    ├── images/
    │   └── <sample-name>.png
    └── polygon_dict.pkl
```
