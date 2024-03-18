---
layout: default-cylinter
title: gating
nav_order: 11
parent: Modules
---

11\. `gating` (optional): this module allows users to classify cell types present in the datset using the [SYLARAS](https://www.sylaras.org/#details) approach to high-dimensional, single-cell gating[[1]](#1). In doing so, users assign a set of manual gating thresholds on a per-marker and per-sample basis using interactive scatter plots of marker (x-axis) x cell segmentation area (y-axis). Gated cells (i.e. those falling to the right of the gate) can be visualized as scatter points in their respective image channel by clicking the `Plot Points` button to confirm accurate gate placement. After an optimal gate has been identified, users will move to the next marker/sample combination in the series by clicking the `Apply Gate and Move to Next Sample` button beneath the scatter plot. If no gate selection is made, all cells in the current plot will be carried forward into downstream analysis. Users may jump between markers and tissues in the series by entering their names into respective fields in the `Arbitrary Sample/Marker Selection` widget at the bottom right of the Napari viewer and clicking the `Enter` button. This can allow for the adjustment of previously defined gates. PDFs showing scatter plots with superimposed gates are stored in the `gating` output directory as a reference which can be updated any time by entering the name of a specific marker in the `Marker Name` field and clicking the `Refresh PDF(s)` button at the bottom right of the Napari viewer; typing "ALL" into the `Marker Name` field will render gated scatter plots for all markers in the analysis. Gates may be re-defined, by removing the metadata associated with particular marker/sample combinations in `cylinter_report.yml` located in the CyLinter output directory specified in `cylinter_config.yml` and re-running the `gating` module with `cylinter --module gating cylinter_config.yml`.

After all gates have been applied, signal intensities are automatically binarized according to the defined gating thresholds such that cells falling to the right of the gate are considered immunopositive, and those falling to the left of the gate are considered immunonegative. Unique Boolean vectors (i.e., binary phenotype profiles) emerging from this procedure are then mapped to biologically-meaningful cell types previously defined in the YAML configuration file (`cylinter_config.yml`). This module can be bypassed by toggling the `gating` parameter to `False` (see YAML configurations below).

### YAML configurations


| Parameter | Default | Description |
| --- | --- | --- |
| `gating` | "False" | (bool) Whether to perform SYLARAS-style gating on single-cell data |
| `channelExclusionsGating` | [ ] | (list of strs) Immunomarkers to exclude from gating. |
| `samplesToRemoveGating` | [ ] | (list of strs) Samples to exclude from gating. |
| `vectorThreshold` | 100 | (int) vizualize Boolean vectors with cell counts >= vectorThreshold |
| `classes` | Tumor: definition: [+pan-CK, +KI67, -aSMA, -CD45] subsets: [CDKN1A] | (dict) Boolean immunophenotype signatures. +marker = immunopositive , -marker = immunonegative |

## References

<a id="1">[1]</a>
Baker GJ. et al. [SYLARAS: A Platform for the Statistical Analysis and Visual Display of Systemic Immunoprofiling Data and Its Application to Glioblastoma](https://www.sciencedirect.com/science/article/pii/S2405471220302854). **Cell Systems** (2020)

