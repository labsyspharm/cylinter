---
layout: default-cylinter
title: curateThumbnails
nav_order: 15
parent: Modules
---

15\. `curateThumbnails`: this module programmatically generates cropped regions of tissue showing reference cells drawn at random from each cluster identified by the [`clustering`]({{ site.baseurl }}/modules/clustering) module. Image contrast settings defined in the [`setContrast`]({{ site.baseurl }}/modules/setContrast) are applied to images for visualization purposes. A single white pixel corresponding to the nuclear centroid of the reference cell is shown in each image as a fiducial reference. Optionally, images can be saved with segmentation outlines superimposed by toggling the `segOutlines` parameter in YAML configuration file (`config.yml`). The number and size of thumbnail examples may be adjusted using the `numThumbnails` and `squareWindowDimension` parameters in `config.yml`. Image galleries for each cluster and gate-based cell type class are saved in `<output_dir>/clustering/2d/thumbnails/` in the case of 2D clustering or `<output_dir>/clustering/3d/thumbnails/` in the case of 3D clustering.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `numThumbnails` | 10 | (int) Number of examples per cluster to be curated. |
| `topMarkersThumbnails` | "channels" | (str) Normalization axis ("channels" or "clusters") used to define highest expressed markers per cluster. |
| `windowDimension` | 35 | (int) Number of pixels from the centroid of the reference cell in x and y dimensions. |
| `segOutlines` | False | (bool) Whether to overlay cell segmentation outlines on thumbnail images. |
