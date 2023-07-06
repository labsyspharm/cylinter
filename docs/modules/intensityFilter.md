---
layout: default-cylinter
title: intensityFilter
nav_order: 3
parent: Modules
---

3\. `intensityFilter`: Out-of-focus and counterstain-oversaturated cell nuclei tend to have altered signal intensity values relative to in-focus and properly counterstained cells. In this module, users interact with histogram widgets of mean counterstain signal intensity to assign upper and lower cutoffs on DNA intensity. Users can then visualize cells between lower and upper cutoffs as scatter points colored by nuclear signal intensity in their respective image by clicking the `Plot Points` button beneath the histogram. Segmentation outlines are provided as a fiducial reference for evaluating cutoff assignments. Selected data points are carried forward into downstream analyses. Closing the Napari window causes the program to apply the current cutoffs and proceed to the next tissue for nuclear intensity cutoff assignment.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `numBinsIntensity` | 50 | (int) Number of bins for DNA intensity histograms. |
