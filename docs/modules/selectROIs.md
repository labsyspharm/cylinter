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

2\. `selectROIs`: [manual](#manual-roi-selection) and [automated](#automated-artifact-detection) tools are used to highlight regions of tissue affected by microscopy artifacts (e.g. illumination aberrations, slide debris, out-of-focus image tiles, mis-registered regions of tissue, etc.).

### Manual ROI selection
Regions of interest (ROIs) are manually drawn around artifacts by clicking on the `Manual ROI Selections (neg.)` image layer in the `layer list` at the left of the Napari viewer then clicking on one of the built-in polygon selection tools from the `layer controls` dock (i.e. circle, square, triangle, or lasso icons above the `layers list`). The mouse button (or track pad) is then clicked and held to outline an artifact in the image window. Clicking the escape key allows for additional ROIs to be drawn. Both positive and negative ROI selection methods are available (see `delint` configuration in `cylinter_config.yml` for details). In the case of negative selection (i.e. `delint=True`, default), cells in ROI boundaries are dropped from the analysis; negative selection is the preferred method for tissues exhibiting diffuse artifacts. Positive selection works best on samples exhibiting large regions of artifact-free tissue that can be highlighted by one or a few ROIs. Cells selected in this case are carried forward into downstream analysis. 

### Automated Artifact Detection
To supplement manual artifact curation, users can choose to run an automated artifact detection (AAD) algorithm on individual image channels by selecting the target channel from the pulldown menu in the `Automated Artifact Detection` widget at the right of the Napari window and clicking the `Compute Artifact Mask` button. Translucent white artifact masks will then appear over regions of tissue that the model flags as putative artifacts. When the `auto` box is checked, the model is run using a reasonable default sensitivity parameter. Sensitivity of the algorithm can be adjusted manually by changing the value in the spinbox labeled `Sensitivity`. Each time the algorithm is run on a given channel it adds two layers to the `layers list` at the left of the Napari viewer. The first layer shows the artifact masks. The second layer shows the seed points corresponding to the different artifacts in the image. Seed points are not visible by default, but can be toggled on by clicking the eye icon shown in the `Artifact Seeds` layer. Individual seed points (and their corresponding artifact masks) can be modified or removed from a given channel by highlighting the `Artifact Seeds` layer, selecting the `arrow icon` in the `layers control` dock to enable point selection mode, and pressing and holding the mouse button to drag over the target seed point to highlight it. Once highlighted, users can fine-tune the artifact mask associated with the seed by changing the `Tolerance` value in the '`Fine-tuning` widget at the right of the Napari viewer or delete the seed entirely by clicking the `x` button in the `layer controls` dock. These AAD tailoring features are designed to give users flexibility over automated artifact masks without the need to re-run the AAD algorithm.

### Workflow
Once all ROIs for a given sample have been generated, users will move to the next sample in the series by clicking the `Apply ROI(s) and Move to Next Sample` button at the top right of the Napari window. If no ROIs are drawn for a given sample, all cells in that tissue will be carried forward into downstream analysis. Users may also jump between samples by entering the name of a given sample in the `Sample Name` field at the right of the Napri viewer to add, delete, or modify manual or automated ROIs of previously analyzed samples or refer to arbitrary tissues in the curation of ROIs for a given samples. ROIs can be added, removed, or modifiedby re-running the `selectROIs` module.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `delintMode` | False | (bool) Whether to drop (True; negative selection) or retain (False; positive selection) cells selected by ROIs. |
| `showAbChannels` | True | (bool) Whether to show all immunomarker channels (True) when Napari is open (may be memory limiting) or show only cycle 1 DNA (False). |
| `samplesForROISelection` | [ ] | (list of strs) Sample names for ROI selection specified according to the first elements of [sampleMetadata]({{ site.baseurl }}/structure/#general-configurations) configuration.
| `autoArtifactDetection` | True | (bool) Whether to display tools for automated artifact detection in Napari window. |
| `artifactDetectionMethod` | "classical" | (str) Algorithm used for automated artifact detection (current option: "classical"). Deep learning method currently under development.