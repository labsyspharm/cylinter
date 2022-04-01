---
layout: default-cylinter
title: curateThumbnails
nav_order: 14
parent: Modules
---

14\. `curateThumbnails`: This module programmatically generates cropped regions of tissue showing cells drawn at random from each cluster identified by the [`clustering`]({{ site.baseurl }}/modules/clustering) module. Image contrast settings defined in the [`setContrast`]({{ site.baseurl }}/modules/setContrast) are applied to images prior to thumbnail curation. A single white pixel corresponding to the nuclear centroid of the reference cell is shown in each image as a fiducial reference. Optionally, images can be saved with segmentation outlines superimposed by toggling the `segOutlines` parameter in `config.yml`. The number and size of thumbnail examples can also be specified using the `numThumbnails` and `squareWindowDimension` parameters in `config.yml`. Image galleries per cluster are saved in `<output_dir>/clustering/thumbnails/`.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `numThumbnails` | 10 | (int) Number of examples per cluster to be curated. |
| `topMarkersThumbnails` | "channels" | (str) Normalization axis ("channels" or "clusters") used to define highest expressed markers per cluster. |
| `squareWindowDimension` | 35 | (int) Number of pixels from the centroid of the reference cell in x and y dimensions. |
| `segOutlines` | False | (bool) Whether to overlay cell segmentation outlines on thumbnail images. |
