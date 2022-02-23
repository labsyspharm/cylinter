---
layout: default-cylinter
title: pruneOutliers
nav_order: 7
parent: Modules
---

7\. `pruneOutliers`: Cells affected by visual artifacts such as antibody aggregates and illumination aberrations appear as outliers in the affected channel and can significantly impact the results of unsupervised cell clustering. In this module, users sensor residual channel outliers from each tissue by applying lower and upper percentile cutoffs on immunomarker signal intensity. Scatter plots (or hexbins, see YAML configurations below) are used to visualize channel-specific intensity distributions versus cell segmentation area before and after cutoff selection (range: 0.0-100.0). Closing the Napari window after cutoff selection causes the program to proceed to the next channel for cutoff assignment. Closing the Napari window without entering cutoffs for a given channel causes all cells to be selected. Once cutoffs have been made for all channels, the module applies channel cutoffs in the order they were applied then re-scales signal intensities for the remaining cells between 0 and 1 as a normalization procedure. Percentile cutoffs are stored as key:value pairs in `<output_dir>/pruning/pruning_dict.pkl`. Remove `pruning_dict.pkl` to re-define cutoffs.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `hexbins` | False | (bool) Whether to use hexbins (True) or scatter plots (False) to plot single-cell signal intensities. Scatter plots allow for higher resolution, but may require longer rendering times.|
| `hexbinGridSize` | 20 | (int) Hexbin grid size when hexins is True. Higher values increase bin resolution. |
