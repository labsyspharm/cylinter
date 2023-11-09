---
layout: default-cylinter
title: intensityFilter
nav_order: 3
parent: Modules
---

3\. `intensityFilter`: out-of-focus or counterstain-oversaturated cell nuclei introduce noise into tissue-derived single-cell data, as out-of-focus cells tend to have unreliable immunomarker signals and oversaturated nuclei tend to be poorly segmented. In this module, users interact with histogram widgets of per-cell counterstain signal intensities to assign upper and lower bounds on DNA signals. Users then visualize selected cells falling between lower and upper cutoffs as scatter points in their respective image colored by DNA signal intensity by clicking the `Plot Points` button. Selected data points are then carried forward into downstream analysis. Users can move to the next sample in the series by clicking the `Apply Gates and Move to Next Sample` button beneath the histogram. If no gate selections are made, all cells in the current sample will be carried forward into downstream analysis. Users may jump between tissues in the series by entering the name of a given sample in the `Sample Name` field of the `Arbitrary Sample Selection` widget at the right of the Napari viewer to adjust gates of previously analyzed tissues or refer to alternative tissues in the curation of cutoffs for a given tissue. Gates may be re-adjusted at any time by re-running the `intensityFilter` module.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `numBinsIntensity` | 50 | (int) Number of bins used to construct DNA intensity histograms. |
