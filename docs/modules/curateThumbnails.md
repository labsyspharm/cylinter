---
layout: default-cylinter
title: curateThumbnails
nav_order: 15
parent: Modules
---

15\. `curateThumbnails`: this module is fully automated. It programmatically generates image patchs of cells drawn at random from each cluster identified in the [clustering module]({{ site.baseurl}}/modules/clustering) and each cell type defined in the [gating module]({{ site.baseurl}}/modules/gating) for visual review. The number of examples shown per cluster is adjusted using the `numThumbnails` parameter in `cylinter_config.yml`. The size of the image window areound the reference cells is controlled by the `squareWindowDimension` parameter in the same configuration file. A white pixel corresponding to the nuclear centroid of the example cell is shown in each image as a reference. Images can be saved with or without segmentation outlines superimposed by toggling the `segOutlines` parameter in configuration file. To facilitate interpretation, only the three most highly expressed protein markers are shown per cluster (based on channel z-scores. Image contrast settings defined in the [setContrast module]({{ site.baseurl }}/modules/setContrast) are applied to improve image appearance. Image galleries for each cluster and gate-based cell type class are saved to the `thumbnails` directory in the `clustering` subdirectory of the main CyLinter output directory. This path is `thumbnails/2d/frequency_stats` in the case of 2D clusterings and `thumbnails/3d/frequency_stats` in the case of 3D clusterings.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `numThumbnails` | 25 | (int) Number of example cells per cluster to be curated. |
| `windowSize` | 30 | (int) Number of pixels from the centroid of the reference cell in x and y dimensions. |
| `segOutlines` | True | (bool) Whether to overlay cell segmentation outlines on thumbnail images. |
