---
layout: default-cylinter
title: cycleCorrelation
nav_order: 5
parent: Modules
---

5\. `cycleCorrelation`: this module is relevant to cyclic imaging technologies (i.e. CyCIF, CODEX, mIHC) and is designed to remove cells that have shifted or become detached from the microscope slide over multi-cycle imaging, as these cells will have false-negative signals for all subsequent markers. Similar to the `intensityFilter` and `areaFilter` modules, users gate on interactive histogram widgets of per-cell signals, only this time histograms represent the log<sub>10</sub>-transformed ratio of DNA intensities between the first and last imaging cycles (log<sub>10</sub>[cycle<sub>1</sub>/cycle<sub>n</sub>]). Lower and upper cutoff sliders are adjusted to select cells with highly-correlated signals at or around zero (log<sub>10</sub>[1/1] = 0). Once lower and upper cutoffs are adjusted, users can visualize selected cells in their corresponding image by clicking the `Plot Points` button. DNA channels for the first and last imaging cycles are shown for reference to visualize cells that have shifted or become detached from the slide between the first and last imaging cycles. Data points between lower and upper cutoffs are carried forward into downstream analysis. Users can move to the next sample in the series by clicking the `Apply Gates and Move to Next Sample` button beneath the histogram. If no gate selections are made, all cells in the current sample will be carried forward into downstream analysis. Users may jump between tissues in the series by entering the name of a given sample in the `Sample Name` field of the `Arbitrary Sample Selection` widget at the right of the Napari viewer to adjust gates of previously analyzed tissues or refer to alternative tissues in the curation of cutoffs for a given tissue. Gates may be re-adjusted at any time by re-running the `cycleCorrelation` module.


### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `numBinsCorrelation` | 50 | (int) Number of bins used to construct DNA<sub>1</sub>/DNA<sub>n</sub> histograms. |
