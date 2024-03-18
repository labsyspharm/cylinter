---
layout: default-cylinter
title: pruneOutliers
nav_order: 7
parent: Modules
---

7\. `pruneOutliers`: cells affected by visual artifacts such as antibody aggregates or illumination aberrations appear as bright outliers in affected channels. Conversely, image background subtraction can have the unintended consequence of creating cells with signal intensities at or below zero that, on image clipping and log-transformation, are far lower than values associated with biologically relevant signals. Both of this scenarios can significantly impact data interpration. In this module, users remove any residual channel outliers from tissues not captured by the `selectROIs` module (e.g., small antibody aggregates) by applying lower and upper percentile cutoffs on marker intensity. Scatter plots (or hexbins, see YAML configurations below) are used to visualize channel-specific intensity distributions before and after cutoffs are applied. Marker intensites are plotted against cell segmentation area which is used as a dumby variable to create 2D plots so that small numbers of outliers can be easily detected. Post-cutoff distributions are shown on a normalized (0-1) x-axis. By entering the name of a given sample in the `Sample Name` field and clicking the `view Outliers` button, users can visualize dim and bright outliers as scatter points (dim = magenta, bright = cyan) in their respective image channels. Users will move to the next channel in the series by clicking the `Apply Cutoffs and Move to Next Marker` button beneath the plots. Note that cells are dropped from the marker channels in an ordered series. Thus, users can elect to re-start outlier removal from a given marker by entering the name of the target channel in the `Re-start from Marker` field and clicking the enter button, but must re-curate outliers in all subsequent channels as well. If no cutoffs are applied for a given marker, all cells in the plots will be carried forward into the analysis of the subsequent marker. To re-define percentile cutoffs, remove the metadata associated with the target channel(s) from `cylinter_report.yml` located in the CyLinter output directory specified in `cylinter_config.yml` and re-run the `pruneOutliers` module with `cylinter --module pruneOutliers cylinter_config.yml`.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `hexbins` | False | (bool) Whether to use hexbins (True) or scatter plots (False) to plot single-cell signal intensities. Scatter plots allow for higher resolution, but may lead to long rendering times with large datasets.|
| `hexbinGridSize` | 20 | (int) The number of hexagons in the x-direction; higher values increase bin resolution. The number of hexagons in the y-direction is chosen such that the hexagons are approximately regular. |
