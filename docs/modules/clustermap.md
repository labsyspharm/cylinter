---
layout: default-cylinter
title: clustermap
nav_order: 13
parent: Modules
---

13\. `clustermap`: this module computes hierarchically clustered heatmaps (i.e. clustermaps) of mean immunomarker signals of clusters identified in the [`clustering`]({{ site.baseurl }}/modules/clustering) module. One is normalized across clusters (row-wise), the other is normalized across channels (column-wise). Plots are saved in `<output_dir>/clustering/2d/` in the case of 2D clustering or `<output_dir>/clustering/3d/` in the case of 3D clustering.

### No YAML configurations
