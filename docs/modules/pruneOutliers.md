---
layout: default
title: pruneOutliers
nav_order: 7
parent: Modules
---

7\. `pruneOutliers`: Cells affected by spurious visual artifacts such as antibody aggregates and illumination aberrations appear as outliers in the affected channel which can significantly impact the results of unsupervised cell clustering. In this module, users sensor channel outliers from each tissue by applying lower and upper percentile cutoffs on immunomarker signal intensity. Scatter plots (or `hexbins` which can reduce plot rendering times for larger datasets) are used to visualize intensity distributions versus cell segmentation area before and after cutoff selection (range: 0.0-100.0). Plots are re-rendered with outliers removed for confirmation or refinement. Closing the Napari window after cutoff selection causes the program to proceed to the next channel for cutoff assignment. Closing the Napari window without entering cutoffs causes all cells to be selected for that channel. After cutoff assignments have been made for all channels, the module applies the cutoffs to the data in the order in which they were curated then re-scales signal intensities for the remaining cells between 0 and 1 as a normalization procedure. Percentile cutoffs are stored as key:value pairs in `<output_dir>/pruning/pruning_dict.pkl`. Remove `pruning_dict.pkl` to re-define cutoffs.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `hexbins` | False (bool) | If True, use hexbins to plot single-cell signal intensities. If False, use scatter plots to plot single-cell signal intensities (allows for higher resolution, but requires longer rendering times).|
| `hexbinGridSize` | 20 (int) | Hexbin grid size (int); higher values increase bin resolution. Only applicable when hexbins parameter is True. |
