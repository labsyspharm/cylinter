---
layout: default-cylinter
title: cycleCorrelation
nav_order: 5
parent: Modules
---

5\. `cycleCorrelation`: this module is relevant to cyclic imaging technologies (e.g., CyCIF, CODEX, mIHC) and is designed to remove cells that have shifted or become detached from the microscope slide over multi-cycle imaging studies, as these cells appear negative for all markers after the movement or loss event. Similar to the `intensityFilter` and `areaFilter` modules, users will gate on interactive histogram widgets of per-cell signals. However, the histograms in this module represent the log<sub>10</sub>-transformed ratio of DNA intensities between the first and last imaging cycles (log<sub>10</sub>[cycle<sub>1</sub>/cycle<sub>n</sub>]). Lower and upper cutoff sliders are adjusted to select cells with highly-correlated signals (typically at or around zero, as log<sub>10</sub>[1/1] = 0). Like in the `intensityFilter` and `areaFilter` modules, Gaussian mixture modeling (GMM) is used to identify initial default cutoffs that can be manually refined. Once lower and upper cutoffs are adjusted, users can visualize selected cells in their corresponding image by clicking the `Plot Points` button. DNA channels for the first and last imaging cycles are shown for reference to visualize cells that have shifted or become detached from the slide between the first and last imaging cycles. Data points between lower and upper cutoffs are carried forward into downstream analysis. Users will move to the next sample in the series by clicking the `Apply Gates and Move to Next Sample` button beneath the histogram. Users may jump between tissues in the series by entering the name of a given sample in the `Sample Name` field of the `Arbitrary Sample Selection` widget at the right of the Napari viewer to adjust thresholds of previously analyzed tissues. To re-define thresholds, remove the metadata associated with the target sample(s) from `cylinter_report.yml` located in the CyLinter output directory specified in `cylinter_config.yml` and re-run the `cycleCorrelation` module with `cylinter --module cycleCorrelation cylinter_config.yml`.


### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `numBinsCorrelation` | 50 | (int) Number of bins used to construct DNA<sub>1</sub>/DNA<sub>n</sub> histograms. |
