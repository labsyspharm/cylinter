---
layout: default-cylinter
title: intensityFilter
nav_order: 3
parent: Modules
---

3\. `intensityFilter`: Cell nuclei that are either out-of-focus or oversaturated with counterstain can introduce noise into single-cell data, as out-of-focus cells tend to have unreliable immunomarker signals and those oversaturated with counterstain are generally poorly segmented. In this module, users interact with histogram widgets of per-cell counterstain signal intensities to assign upper and lower bounds on DNA intensity. Users can then visualize cells between lower and upper cutoffs as scatter points colored by counterstain signal intensity in their respective image by clicking the `Plot Points` button beneath the histogram. Segmentation outlines are provided in the `layer list` at the left of the Napari viewer as a fiducial reference for evaluating cutoff assignments. Data points between the lower and upper cutoffs are carried forward into downstream analyses. Users can move to the next sample in the batch series by clicking the `Apply Gates and Move to Next Sample` button at the top right of the Napari window. If no ROIs are drawn for a given sample, all cells in that tissue will be carried forward into downstream analysis. Alternatively, users can jump between samples by entering their name in the `Sample Name` field at the right of the Napri viewer to adjust gates of previously analyzed tissues or refer to tissues further in the batch series when considering gate selections for another tissue. Gates can be re-adjusted at any time by re-running the pipeline from the `intensityFilter` module.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `numBinsIntensity` | 50 | (int) Number of bins in DNA intensity histograms. |
