---
layout: default-cylinter
title: areaFilter
nav_order: 4
parent: Modules
---

4\. `areaFilter`: Cell segmentation errors introduce significant noise into tissue-derived, single-cell data. In this module, users assign lower and upper cutoffs on the segmentation area (measured in pixels) of cells comprising each tissue as calculated by the quantification module in [MCMICRO](https://mcmicro.org/) (or other feature extraction method) to remove under- and over-segmented cells. The module functions in the same way as the `intensityFilter` module through assigmnet of lower and upper cutoffs on per-cell segmentation area. After cutoffs have been set, users can visualize selected cells in their corresponding tissue by clicking the `Plot Points` button beneath the histogram. Segmentation outlines are provided in the `layer list` at the left of the Napari viewer as a fiducial reference for evaluating cutoff assignments. Data points between the lower and upper cutoffs are carried forward into downstream analyses. Users can move to the next sample in the batch series by clicking the `Apply Gates and Move to Next Sample` button at the top right of the Napari window. If no ROIs are drawn for a given sample, all cells in that tissue will be carried forward into downstream analysis. Alternatively, users can jump between samples by entering their name in the `Sample Name` field at the right of the Napri viewer to adjust gates of previously analyzed tissues or refer to tissues further in the batch series when considering gate selections for another tissue. Gates can be re-adjusted at any time by re-running the pipeline from the `areaFilter` module.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `numBinsArea` | 50 | (int) Number of bins in DNA area histograms. |
