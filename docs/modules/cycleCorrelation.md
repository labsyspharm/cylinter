---
layout: default-cylinter
title: cycleCorrelation
nav_order: 5
parent: Modules
---

5\. `cycleCorrelation`: This module pertains to cyclic imaging technologies (i.e. CyCIF) and is designed to remove cells that have shifted or become detached from the microscope slide over the course of imaging. This phenomenon, referred to as "cell dropout", leads to false-negative signals for immunomarkers probed for subsequent to the dropout event. Users apply lower and upper cutoffs on the log<sub>10</sub>-transformed ratio of DNA intensities between cells at the first and last imaging cycles (log<sub>10</sub>[cycle<sub>1</sub>/cycle<sub>n</sub>]) for each tissue. Users adjust sliders to select cells with highly-correlated signals centered near zero (log<sub>10</sub>[1/1] = 0). After gate placement, users can visualize selected cells in their corresponding tissue by clicking the `Plot Points` button beneath the histogram. DNA channels for the first and last imaging cycles are shown together with scatter points at the nuclear centroids of selected cells to be carried forward into downstream analyses.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `numBinsCorrelation` | 50 | (int) Number of bins for DNA1/DNAn histograms. |
