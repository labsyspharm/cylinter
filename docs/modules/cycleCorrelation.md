---
layout: default
title: cycleCorrelation
nav_order: 5
parent: Modules
---

5\. `cycleCorrelation`: This module pertains to cyclic imaging technologies such as t-CyCIF and aids in the removal of cells that have shifted or become detached from the microscope slide over the course of imaging. This phenomenon, referred to as "cell dropout", leads to false-negative signals for all immunomarkers probed for after the dropout event. To remove cell dropout, users apply lower and upper cutoffs on the log<sub>10</sub>-transformed ratio of DNA intensities between cells at the first and last imaging cycles (log<sub>10</sub>[cycle<sub>1</sub>/cycle<sub>n</sub>]) for each tissue. Users adjust sliders to select cells with highly-correlated signals centered around zero (log<sub>10</sub>[1/1] = 0). After gate placement, users visualize where selected cells reside in a given tissue using the Napari image viewer. DNA channels for the first and last imaging cycles are shown together with scatter points at the nuclear centroids of selected cells that will be carried forward into downstream analyses. Tissues with little cell loss may be gated by applying a single Y-axis cutoff on cell count. In this method, users first visualize histograms for all tissues in the analysis, then select a single cutoff. Histograms are re-rendered with a horizontal line at the assigned cutoff for confirmation or refinement. Y-axis gating may be toggled using the `yAxisGating` parameter in `config.yml` (see below). While y-axis gating can expedite the removal of cell dropout, it may be unsuitable for tissue samples exhibiting large amounts of cell loss.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `yAxisGating` | False (bool) | If True, apply a single y-axis gate on counts of log10[cycle1/cycleN] histograms of cells per tissue. If False, apply lower and upper gates on log10[cycle1/cycleN] histogram bins. |
