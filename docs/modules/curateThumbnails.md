---
layout: default
title: curateThumbnails
nav_order: 13
parent: Modules
---

13\. `curateThumbnails`: This module programmatically generates cropped regions of tissue showing cells drawn at random from each cluster identified by the [`clustering`]({{ site.baseurl }}/modules/clustering) module. Image contrast settings defined in the [`setContrast`]({{ site.baseurl }}/modules/setContrast) are applied to images prior to thumbnail curation. A single white pixel corresponding to the nuclear centroid of the reference cell is shown in each image as a fiducial reference. Optionally, images can be saved with segmentation outlines superimposed by toggling the `segOutlines` parameter in `config.yml`. The number and size of thumbnail examples can also be specified using the `numThumbnails` and `squareWindowDimension` parameters in `config.yml`. Image galleries per cluster are saved in `<output_dir>/clustering/thumbnails/`.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `numThumbnails` | 10 (int) | Number of examples per cluster to be curated |
| `squareWindowDimension` | 35 (int) | Number of pixels from the reference cell centroid in the x and y directions |
| `segOutlines` | False (bool) | If True, overlay cell segmentation outlines on thumbnail images |
