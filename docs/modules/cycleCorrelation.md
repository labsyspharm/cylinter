---
layout: default-cylinter
title: cycleCorrelation
nav_order: 5
parent: Modules
---

5\. `cycleCorrelation`: This module pertains to cyclic imaging technologies (i.e. CyCIF) and is designed to remove cells that have shifted or become detached from the microscope slide over the course of imaging. This phenomenon, referred to as "cell dropout", leads to false-negative signals for immunomarkers probed for subsequent to the dropout event. Users apply lower and upper cutoffs on the log<sub>10</sub>-transformed ratio of DNA intensities between cells at the first and last imaging cycles (log<sub>10</sub>[cycle<sub>1</sub>/cycle<sub>n</sub>]) for each tissue. Users adjust sliders to select cells with highly-correlated signals centered near zero (log<sub>10</sub>[1/1] = 0). After gate placement, users can visualize selected cells in their corresponding tissue by clicking the `Plot Points` button beneath the histogram. DNA channels for the first and last imaging cycles are shown together with scatter points at the nuclear centroids of selected cells to be carried forward into downstream analyses. Tissues with little cell loss may be gated by applying a single Y-axis cutoff on cell count (see YAML configuration below). In this method, users are first presented with histograms for all tissues in the analysis to determine a single cutoff to be applied to all samples. Histograms are re-rendered with a horizontal line at the assigned cutoff for confirmation or refinement. While y-axis gating can expedite the removal of cell dropout, it may be unsuitable for tissue samples exhibiting large amounts of cell loss.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `yAxisGating` | False | (bool) Whether to apply a single y-axis counts cutoff (True) or lower and upper bins cutoffs of log10(cycle1/cycleN) histograms of cells per tissue. |
