---
layout: default-cylinter
title: pruneOutliers
nav_order: 7
parent: Modules
---

7\. `pruneOutliers`: cells affected by visual artifacts such as antibody aggregates or illumination aberrations appear as bright outliers in the affected channel. On the other hand, strong image background subtraction (which can significantly enhance antibody signal-to-noise) can have the unintended consequence of creating cells with signal intensities near (or equal to) zero, that after log transformation, are far lower than values associated with the majority of cells in an image; both of these scenarios can significantly impact the results of unsupervised cell clustering or manual gating in the identification of cell states. Thus, in this module, users remove residual channel outliers from tissues not captured by the `selectROIs` module by applying lower and upper percentile cutoffs on marker intensity. Scatter plots (or hexbins, see YAML configurations below) are used to visualize channel-specific intensity distributions versus cell segmentation area before and after cutoffs are applied. Post-cutoff distributions are shown on a rescaled x-axis between 0-1. Users can visualize dim and bright outliers as color-coded scatter points (dim = magenta, bright = cyan) in their respective image by entering the name of a given sample in the `Sample Name` field and clicking the `view Outliers` button. Users move to the next channel in the series by clicking the `Apply Cutoffs and Move to Next Marker` button beneath the plots. Note that cells are dropped after cutoffs are applied to each channel in series. Users can elect to re-start the data curation process from a given marker by entering the name of the target channel in the `Re-start from Marker` field and clicking the enter button. If no cutoffs are applied for a given marker, all remaining cells will be carried forward in the analysis of the subsequent marker. Cutoffs can be re-adjusted at any time by re-running the pipeline from the `pruneOutliers` module.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `hexbins` | False | (bool) Whether to use hexbins (True) or scatter plots (False) to plot single-cell signal intensities. Scatter plots allow for higher resolution, but may lead to long rendering times.|
| `hexbinGridSize` | 20 | (int) Hexbin grid size when hexins is `True`. Higher values increase bin resolution. |
