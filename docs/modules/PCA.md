---
layout: default
title: PCA
nav_order: 9
parent: Modules
---

9\. `PCA`: This module performs principle component analysis (PCA) on mean immunomarker intensities of cells for each tissue sample in the analysis. The PCA scores plot is saved to `<output_dir/PCA/pcaScoresPlot.pdf>`.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `channelExclusionsPCA` | [ ] (list) | Immunomarkers to exclude from PCA analysis and all subsequent modules (strs) |
| `samplesToRemovePCA` | [ ] (list) | Samples to exclude from PCA analysis and all subsequent modules (strs) |
| `dimensionPCA` | 2 (int) | Number of PCs to compute (int, typically 2) |
| `pointSize` | 90.0 (float) | Scatter point size for sample scores plot |
| `labelPoints` | True (bool) | Annotate scatter points with condition abbreviations specified in the sample_metadata dict |
| `distanceCutoff` | 0.15 (bool) | Maximum distance between data points in PCA scores plot to be annotated with a common label. Useful for increasing visual clarity of PCA plots containing many data points. Only applicable when labelPoints is True. |
| `conditionsToSilhouette` | [ ] (list) | List of abbreviated condition names whose corresponding scores plot points will be greyed out, left unannotated, and sent to the back of the plot (zorder). Useful for increasing visual clarity of PCA plots containing many data points. |
