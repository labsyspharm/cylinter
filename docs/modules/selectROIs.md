---
layout: default-cylinter
title: selectROIs
nav_order: 2
parent: Modules
---

### Description
2\. `selectROIs`: Polygon selection tools built into the Napari image viewer (i.e. triangle, circle, and square icons) are used to gate on tissue regions of interest (ROIs). Positive and negative selection modes are available. In positive selection, cells within ROI boundaries are carried forward into downstream analysis; in negative selection, selected cells are dropped from the analysis. Negative selection is preferred when datasets exhibit diffuse artifacts. Closing the Napari window after drawing ROIs for a given tissue causes the program to advance to the next tissue for ROI selection. Closing the window without drawing ROIs causes all cells in that tissue to be selected for downstream analysis. ROI polygon vertices are stored as key:value pairs in `<output_dir>/ROIs/polygon_dict.pkl`. This file must be removed (or individual key:value pairs must be deleted) to re-define ROIs.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `delintMode` | False | (bool) Whether to drop (True; negative selection) or retain (False; positive selection) cells selected by ROIs. |
| `showAbChannels` | True | (bool) Whether to show all immunomarker channels (True) when Napari is open (may be memory limiting) or show only cycle 1 DNA (False). |
| `samplesForROISelection` | [ ] | (list of strs) Sample names for ROI selection specified according to the first elements of [sampleMetadata]({{ site.baseurl }}/workflow/input#general-configurations) configuration. |
