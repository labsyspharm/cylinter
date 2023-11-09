---
layout: default-cylinter
title: gating
nav_order: 11
parent: Modules
---

11\. `gating` (optional): this module allows users to perform manual gating on per-cell marker instensities using the [SYLARAS](https://www.sciencedirect.com/science/article/pii/S2405471220302854) approach to high-dimensional, single-cell gating (see the "Cell Subset Identification" section [here](https://www.sylaras.org/#details) for details). In doing so, users interact with immunomarker x cell segmentation area scatter plots to apply one-dimensional gates along the immunomarker axis (x-axis) on a per-sample basis. Users can visualize gated cells (i.e. those falling to the right of the gate) as scatter points in their respective image channels by clicking the `Plot Points` button. After an optimal gate has been identified, the gate is applied and the program moves to the next sample in the series by clicking the `Apply Gate and Move to Next Sample` button beneath the scatter plot. If no gate selection is made, all cells in the current sample will be carried forward into downstream analysis. Users may jump between channels and tissues in the series by entering their names in respective fields in the `Arbitrary Sample/Marker Selection` widget at the bottom right of the Napari viewer and clicking the `Enter` button to adjust gates of previously analyzed tissues or refer to alternative tissues in the curation of gates for a target tissue. PDF pages are stored in the `gating` output folder to document gate placement and can be updated at any time by entering the name of a specific marker in the `Marker Name` field and clicking the `Refresh PDF(s)` button at the bottom right of the Napari viewer; typing "ALL" into the `Marker Name` field will render gated scatter plots for all markers. Gates may be re-adjusted at any time by re-running the `gating` module.

After gates have been applied to all combinations of markers and samples in the series, signal intensities are programmatically binarized according to the gates such that cells falling to the right of the gate are considered immunopositive and take the value 1, and those falling to the left of the gate are considered immunonegative and take the value 0. Unique Boolean vectors emerging from this process are then mapped to biologically-meaningful cell types previously defined in the YAML configuration file (`config.yml`). This module can be bypassed by toggling the `gating` parameter to `False` (see YAML configurations below).

### YAML configurations


| Parameter | Default | Description |
| --- | --- | --- |
| `gating` | "True" | (bool) Whether to perform SYLARAS-style gating on single-cell data (Cell Syst. 2020 Sep 23;11(3):272-285.e9 PMID: 32898474) |
| `channelExclusionsGating` | [ ] | (list of strs) Immunomarkers to exclude from gating. |
| `samplesToRemoveGating` | [ ] | (list of strs) Samples to exclude from gating. |
| `vectorThreshold` | 100 | (int) vizualize Boolean vectors with cell counts >= vectorThreshold |
| `classes` | Tumor: definition: [+pan-CK, +KI67] subsets: [] | (dict) Boolean immunophenotype signatures. +marker = immunopositive , -marker = immunonegative, marker = don't care |
