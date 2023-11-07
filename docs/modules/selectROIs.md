---
layout: default-cylinter
title: selectROIs
nav_order: 2
parent: Modules
---

{: .no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

2\. `selectROIs`: A combination of [manual](#manual-roi-selection) and [automated](#automated-artifact-detection) tools are used to highlight regions of tissue affected by microscopy artifacts (e.g. illumination aberrations, slide debris, out-of-focus image tiles, mis-registered regions of tissue, etc.).

### Manual ROI selection
Regions of interest (ROIs) can be drawn manually around artifacts in specific image channels by first selecting the `Manual ROI Selections (neg.)` image layer in the `layer list` at the left of the Napari viewer, selecting one of the built-in polygon selection tools from the `layer controls` (i.e. circle, square, triangle triangle, or lasso icon), then clicking and holding the mouse or track pad and outlining the artifact. Both positive and negative ROI selection methods are available (see `delint` configuration in `config.yml` for details). In the case of negative selection (`delint=True`), cells in ROI boundries are dropped from the analysis; negative selection is the preferred method for tissues exhibiting diffuse artifacts. By contrast, positive selection, works best for tissues exhibiting large regions of artifact-free tissue that can be highlighted by one or a few ROIs whose corresponding cells are carried forward in the analysis. 

### Automated Artifact Detection
To faciliate artifact curation, users can also run an automated artifact detection algorithm on different antibody channels listed in the pulldown menu at the right of the Napari window in the `Automated Artifact Detection` widget by selecting the channel of interest and clicking the `Compute Artifact Mask` button. Translucent artifact masks will then appear over regions of tissue that the model flags as putative artifacts. By default the model is run using a reasonable default sensitivity parameter (when `auto` is checked), but the  sensitivity of the algorithm can be adjusted by changing the value in the `Sensitivity` spinbox. The algorithm produces two additional layers per channel that the algorithm is run on which are added to the `layers list`. The first shows the seed points of artifacts flagging by the algorithm; the second shows their corresponding artifact masks. Highlighting individual seed points for a particular channel after selecting the target `Artifacts Seeds` and selecting the `Select shapes`(open arrow icon) in the `layer controls`, allows the user to fine-tune the artifact mask associated with that seed point by changing the `Tolerance` value in '`Fine-tuning` widget at the right of the Napari viewer.   

Users can jump between samples by entering the name of an arbitrary sample in the `Sample Name` field in the `Arbitrary Sample Selection` widget to add, delete, or modify ROIs of previously curated tissues or refer to other tissues as reference for curating ROIs in another ROIs in another. After all manual and automated curations have been made for a given sample, users can move to the next sample in the batch by clicking the `Apply ROI(s) and Move to Net Sample` at the top right of the Napari window. If no ROIs are drawn for a given sample, all cells in that tissue will be carried forward into downstream modules. ROIs can be added, removed, or modified at any time by rerunning the this module.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `delintMode` | False | (bool) Whether to drop (True; negative selection) or retain (False; positive selection) cells selected by ROIs. |
| `showAbChannels` | True | (bool) Whether to show all immunomarker channels (True) when Napari is open (may be memory limiting) or show only cycle 1 DNA (False). |
| `samplesForROISelection` | [ ] | (list of strs) Sample names for ROI selection specified according to the first elements of [sampleMetadata]({{ site.baseurl }}/workflow/input#general-configurations) configuration.
| `autoArtifactDetection` | True | (bool) Whether to display tools for automated artifact detection in Napari window. |
| `artifactDetectionMethod` | "classical" | (str) Algorithm used for automated artifact detection (current option: "classical"). Multi-layer perceptron method ("MLP") currently under development.