---
layout: default
title: selectROIs
nav_order: 2
parent: Modules
---

### Description
2\. `selectROIs`: Polygon selection tools built into the Napari image viewer (i.e. triangle, circle, and square icons) are used to gate on tissue regions of interest (ROIs). Positive and negative selection modes are available (see YAML configurations below). In positive selection, cells within ROI boundaries are carried forward into downstream analysis. In negative selection, selected cells are dropped from analysis. Negative selection is preferred when datasets exhibit diffuse artifacts. Closing the Napari window after drawing ROIs for a given tissue causes the program to advance to the next tissue in the batch. Closing the window without drawing ROIs causes all cells to be selected for that tissue. ROI polygon vertices are stored as key:value pairs in `<output_dir>/ROI/polygon_dict.pkl`. This file must be removed (or individual key:value pairs must be deleted) in order to re-define ROIs.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `delintMode` | False (bool) | Whether to retain (False) or drop (True) cells in selected ROIs. False = positive ROI selection; True = negative ROI selection |
| `showAbChannels` | True (bool) | True = show all immunomarker channels when Napari is open (may be memory limiting). False = show cycle1 DNA channel only |
| `samplesForROISelection` | [ ] (list) | List a sample names (strs) for ROI selection |
