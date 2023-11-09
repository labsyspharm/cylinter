---
layout: default-cylinter
title: areaFilter
nav_order: 4
parent: Modules
---

4\. `areaFilter`: cell segmentation errors introduce significant noise into tissue-derived, single-cell data. In this module, users assign lower and upper cutoffs on the segmentation area (pixel count; as calculated by the quantification module in [MCMICRO](https://mcmicro.org/parameters/core.html#mcquant) or other feature extraction method) of individual cell segmentation instances in a given tissue to remove under- and over-segmented cells. The function of this module is analogous to the `intensityFilter` module in that it allows users to assign lower and upper cutoffs on per-cell segmentation area by interacting with histogram widgets. Once cutoffs have been adjusted for a given sample, users can visualize selected cells in their corresponding image by clicking the `Plot Points` button. Segmentation outlines are provided in the `layer list` at the left of the Napari viewer as a reference for evaluating segmentation quality. Data points falling between lower and upper sliders are carried forward into downstream analysis. Users may move to the next sample in the series by clicking the `Apply Gates and Move to Next Sample` button beneath the histogram. If no gate selections are made, all cells in the current sample will be carried forward into downstream analysis. Users may jump between tissues in the series by entering the name of a given sample in the `Sample Name` field of the `Arbitrary Sample Selection` widget at the right of the Napari viewer to adjust gates of previously analyzed tissues or refer to alternative tissues in the curation of cutoffs for a given tissue. Gates may be re-adjusted at any time by re-running the `areaFilter` module.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `numBinsArea` | 50 | (int) Number of bins used to construct DNA area histograms. |
