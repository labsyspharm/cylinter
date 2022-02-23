---
layout: default-cylinter
title: PCA
nav_order: 9
parent: Modules
---

9\. `PCA`: This module performs principle component analysis (PCA) on mean immunomarker intensities of cells in each tissue. The PCA scores plot is saved to `<output_dir/PCA/pcaScoresPlot.pdf>`.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `channelExclusionsPCA` | [ ] | (list of strs) Immunomarkers to exclude from PCA analysis. |
| `samplesToRemovePCA` | [ ] | (list of strs) Samples to exclude from PCA analysis. |
| `dimensionPCA` | 2 | (int) Number of PCs to compute. |
| `pointSize` | 90.0 | (float) Scatter point size for sample scores plot. |
| `labelPoints` | True | (bool) Annotate scatter points with condition abbreviations from sampleMetadata configuration. |
| `distanceCutoff` | 0.15 | (float) Maximum distance between data points in PCA scores plot to be annotated with a common label. Useful for increasing visual clarity of PCA plots containing many data points. Applicable when labelPoints is True. |
| `conditionsToSilhouette` | [ ] | (list of strs) Abbreviated condition names whose corresponding scores plot points will be greyed out, left unannotated, and sent to the back of the plot (zorder). Useful for increasing visual clarity of PCA plots containing many data points. |
