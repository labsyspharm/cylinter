---
layout: default-cylinter
title: pruneOutliers
nav_order: 7
parent: Modules
---

7\. `pruneOutliers`: Cells affected by bright visual artifacts such as antibody aggregates and illumination aberrations or strong image background subtraction appear as outliers in the affected channel and can significantly impact the results of unsupervised cell clustering or manual gating in the identification of cell states. In this module, users remove any residual channel outliers from tissues by applying lower and upper percentile cutoffs  on immunomarker signal intensity (range: 0.0-100.0). Scatter plots (or hexbins, see YAML configurations below) are used to visualize channel-specific intensity distributions versus cell segmentation area before and after cutoff selection. Users can visualize dim and bright outliers as color-coded scatter points (dim= magenta, bright=cyan) in their respective image by entering the name of a sample in the `Sample Name` field and clicking the `view Outliers` button. Users move to the next channel in the series by clicking the `Apply Cutoffs and Move to Next Marker` button beneath the plots. If sliders are not adjusted for a given marker, all remaining cells will be carried forward in the analysis and shown in the distributions for the next marker. Note that in this module, cells are dropped in the order in whcih the markers are analyzed. Thus, users can elect to re-start the data curation process from a given marker by entering the name of a given channel in the `Re-start from Marker` field and clicking the enter button. Once cutoffs have been made for all markers, channel cutoffs are applied in the order in which they were applied and signals are re-scaled between 0 and 1 as a normalization procedure. Cutoffs can be re-adjusted at any time by re-running the pipeline from the `pruneOutliers` module.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `hexbins` | False | (bool) Whether to use hexbins (True) or scatter plots (False) to plot single-cell signal intensities. Scatter plots allow for higher resolution, but may require longer rendering times.|
| `hexbinGridSize` | 20 | (int) Hexbin grid size when hexins is True. Higher values increase bin resolution. |
