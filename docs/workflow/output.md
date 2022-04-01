---
layout: default-cylinter
title: Output
nav_order: 2
parent: Workflow
---

# Output File Structure

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
│   thumbnails/
│   ├── cluster<cluster-ID>_thumbnails.pdf
│   ├── completed_clusters.txt
│   zarrs/
│   ├── clus<cluster-ID>_<sample-name>_win<thumbnail-size>.zarr/
│   └── clus<cluster-ID>_<sample-name>_win<thumbnail-size>_seg.zarr/
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

In the above example, `<sample-name>` corresponds to the names given to tissue samples for a given analysis. `<test-name>` refers to the binary declarations specified as the fourth value elements in the `sampleMetadata` dictionary (see [general configurations](input#yaml-configuration-file)). `<chunk>` refers to a slice of the combined single-cell feature table for QC status reclassification. `<cluster-ID>` refers to clusters identified in by the `clustering` module. `<thumbnail-size>` refers to the square dimension (in pixels) of the thumbnail images cropped from multiplex images by the `curateThumbnails` module, and `<channel>` refers to the immunomarkers used in the study.
