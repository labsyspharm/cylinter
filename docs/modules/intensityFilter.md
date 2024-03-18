---
layout: default-cylinter
title: intensityFilter
nav_order: 3
parent: Modules
---

3\. `intensityFilter`: out-of-focus cells and those oversaturated with nuclear counterstain introduce noise into image-derived, single-cell data. This is because out-of-focus cells tend to have unreliable immunomarker signals and oversaturated nuclei tend to be poorly segmented. In this module, users interact with histogram widgets of per-cell counterstain signal intensities to assign upper and lower bounds on DNA signal intensity. Gaussian mixture modeling (GMM) is used to identify default cutoffs that can be manually refined. Users can visualize cells falling between lower and upper cutoffs as scatter points in their respective image colored by DNA signal intensity by clicking the `Plot Points` button. Selected data points are then carried forward into downstream analysis. Users will move to the next sample in the series by clicking the `Apply Gates and Move to Next Sample` button beneath the histogram. Users may jump between tissues in the series by entering the name of a given sample in the `Sample Name` field of the `Arbitrary Sample Selection` widget at the right of the Napari viewer to adjust thresholds of previously analyzed tissues. To re-define thresholds, remove the metadata associated with the target sample(s) from `cylinter_report.yml` located in the CyLinter output directory specified in `cylinter_config.yml` and re-run the `intensityFilter` module with `cylinter --module intensityFilter cylinter_config.yml`.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `numBinsIntensity` | 50 | (int) Number of bins used to construct DNA intensity histograms. |
