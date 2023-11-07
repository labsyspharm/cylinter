---
layout: default-cylinter
title: Output
nav_order: 2
parent: Workflow
---

# Output File Structure

In the below example, `<sample-name>` corresponds to the names given to tissue samples for a given analysis. `<test-name>` refers to the binary declarations specified as the fourth value elements in the `sampleMetadata` dictionary (see [general configurations](input#yaml-configuration-file)). `<chunk>` refers to a slice of the combined single-cell feature table for QC status reclassification. `<cluster-ID>` refers to clusters identified in by the `clustering` module. `<thumbnail-size>` refers to the square dimension (in pixels) of the thumbnail images cropped from multiplex images by the `curateThumbnails` module, and `<channel>` refers to the immunomarkers used in the study.

``` bash
<OUTPUT DIR>
├── area/
│   ├── cutoffs.pkl
│   plots/
│       ├── <sample1>.pdf
│       └── <sample2>.pdf
├── checkpoints/
│   ├── aggregateData.parquet
│   ├── areaFilter.parquet
│   ├── clustering.csv
│   ├── clustering.parquet
│   ├── clustermap.parquet
│   ├── curateThumbnails.parquet
│   ├── cycleCorrelation.parquet
│   ├── frequencyStats.parquet
│   ├── gating.parquet
│   ├── intensityFilter.parquet
│   ├── logTransform.parquet
│   ├── metaQC.parquet
│   ├── PCA.parquet
│   ├── pruneOutliers.parquet
│   ├── selectROIs.parquet
│   └── setContrast.parquet
├── clustering/
│   2d/
│    ├── clustermap_cluster_2d_norm_channels.pdf
│    ├── clustermap_cluster_2d_norm_clusters.pdf
│    ├── emb_channels.png
│    ├── embedding.npy
│    frequency_stats/
│       ├── class/
│           ├── CANCER-TRUE_v_CANCER-FALSE/
│               ├── plot.pdf
│               ├── stats_sig.csv
│               └── stats_total.csv
│       └── cluster_2d
│           ├── CANCER-TRUE_v_CANCER-FALSE/
│               ├── plot.pdf
│               ├── stats_sig.csv
│               └── stats_total.csv
│    ├── MCS.txt
│    ├── ridgeplots.pdf
│    thumbnails/
│       ├── class/
│           ├── completed.txt
│           ├── Tumor.pdf
│           └── zarrs/
│                   ├── ...
│       └── cluster_2d
│           ├── 0.pdf
│           ├── 1.pdf
│           ├── 2.pdf
│           ├── 3.pdf
│           ├── 4.pdf
│           ├── completed.txt
│           └── zarrs/
│                   ├── ...
│    ├── UMAP_200_silho.png
│    └── UMAP_200.png
├── contrast/
│   └── contrast_limits.yml
├── cycles/
│   ├── cutoffs.pkl
│   ├── plots/
│       ├── 1.pdf
│       ├── 15.pdf
│       ├── 18.pdf
│       ├── 68.pdf
│       ├── correlation.png
│       └── correlation(sample).png
├── gating/
│   ├── cutoffs.pkl
│   ├── plots/
│       ├── 1.pdf
│       ├── 15.pdf
│       ├── 18.pdf
│       ├── 68.pdf
│       ├── correlation.png
│       └── correlation(sample).png
├── intensity/
│   ├── <sample-name>.pdf
│   └── idxs_to_drop.csv
├── metaQC/
│   ├── <chunk>/
│   ├── censored_by_stage.pdf
│   ├── chunk_index.txt
│   ├── MCS.txt
│   ├── QCData.pkl
│   ├── reclass_storage_dict.pkl
│   ├── RECLASS_TUPLE.txt
│   └── UMAP_<min_cluster_size>.png
├── PCA/
│   └── pcaScoresPlot.pdf
├── pruning/
│   ├── <channel>_pruned_rescaled.png
│   └── <channel>_raw.png
│   ├── data_copy1.parquet
│   └── pruning_dict.csv
└── ROIs/
    ├── images/
    │   └── <sample-name>.png
    └── polygon_dict.pkl
```
