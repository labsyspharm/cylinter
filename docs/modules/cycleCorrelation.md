---
layout: default-cylinter
title: cycleCorrelation
nav_order: 5
parent: Modules
---

5\. `cycleCorrelation`: This module pertains to cyclic imaging technologies (i.e. CyCIF, CODEX, mIHC) and is designed to remove cells that have shifted or became detached from the microscope slide over the course of imaging. This phenomenon, referred to as "cell dropout", leads to false-negative signals in cells for all immunomarkers used subsequent to the dropout event. In this module, users apply lower and upper cutoffs on the log<sub>10</sub>-transformed ratio of DNA intensities between cells at the first and last imaging cycles (log<sub>10</sub>[cycle<sub>1</sub>/cycle<sub>n</sub>]) for each tissue. Sliders are adjusted to select cells with highly-correlated signals at and around zero (log<sub>10</sub>[1/1] = 0). After lower and upper cutoffs are set, users can visualize selected cells in their corresponding tissue by clicking the `Plot Points` button beneath the histogram. DNA channels for the first and last imaging cycles are shown for reference. Data points between the lower and upper cutoffs are carried forward into downstream analyses. Users can move to the next sample in the batch series by clicking the `Apply Gates and Move to Next Sample` button at the top right of the Napari window. If no ROIs are drawn for a given sample, all cells in that tissue will be carried forward into downstream analysis. Alternatively, users can jump between samples by entering their name in the `Sample Name` field at the right of the Napri viewer to adjust gates of previously analyzed tissues or refer to tissues further in the batch series when considering gate selections for another tissue. Gates can be re-adjusted at any time by re-running the pipeline from the `cycleCorrelation` module.


### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `numBinsCorrelation` | 50 | (int) Number of bins in DNA1/DNAn histograms. |
