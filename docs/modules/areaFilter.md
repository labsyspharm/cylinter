---
layout: default-cylinter
title: areaFilter
nav_order: 4
parent: Modules
---

4\. `areaFilter`: cell segmentation errors can be a significant source of noise in image-derived, single-cell data. In this module, users assign lower and upper bounds on cell segmentation area (pixel count) to remove severely under- and over-segmented cells. Cell segmentation area is a standard component of [MCMICRO](https://mcmicro.org/parameters/core.html#mcquant) feature table output and is calculated using [skimage.measure.regionprops()](https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.regionprops). This module functions similar to the `intensityFilter` module in that it allows users to assign lower and upper thresholds on interactive histogram widgets of per-cell data. Gaussian mixture modeling (GMM) assist in identifying default cutoffs that can be manually refined. Once thresholds have been adjusted for a given sample, users can visualize selected cells in their corresponding image by clicking the `Plot Points` button. Segmentation outlines are provided in the `layer list` at the left of the Napari viewer as a reference for evaluating segmentation quality. Data points falling between lower and upper sliders are carried forward into downstream analysis. Users will move to the next sample in the series by clicking the `Apply Gates and Move to Next Sample` button beneath the histogram. Users may jump between tissues in the series by entering the name of a given sample in the `Sample Name` field of the `Arbitrary Sample Selection` widget at the right of the Napari viewer to adjust thresholds of previously analyzed tissues. To re-define thresholds, remove the metadata associated with the target sample(s) from `cylinter_report.yml` located in the CyLinter output directory specified in `cylinter_config.yml` and re-run the `areaFilter` module with `cylinter --module areaFilter cylinter_config.yml`.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `numBinsArea` | 50 | (int) Number of bins used to construct DNA area histograms. |
