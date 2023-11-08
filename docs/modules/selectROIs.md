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
Regions of interest (ROIs) can be manually drawn around channel-specific artifacts by selecting the `Manual ROI Selections (neg.)` image layer in the `layer list` at the left of the Napari viewer, selecting one of the built-in polygon selection tools from the `layer controls` dock (i.e. circle, square, triangle, or lasso icons above the `layers list`), then pressing and holding the mouse (or track pad) and outlining the artifact in the image window. Both positive and negative ROI selection methods are available (see `delint` configuration in `config.yml` for details). In the case of negative selection (i.e. `delint=True`), cells in ROI boundaries are dropped from the analysis; negative selection is the preferred method for tissues exhibiting diffuse artifacts. By contrast, positive selection works best with tissues exhibiting large regions of artifact-free tissue that can be highlighted by one or a few ROIs. Selected cells in this case are the only ones carried forward into downstream analysis. 

### Automated Artifact Detection
To facilitate artifact curation, users can also run an automated artifact detection (AAD) algorithm on individual antibody channels by first selecting a target channel from the pulldown menu in the `Automated Artifact Detection` widget at the right of the Napari window then clicking the `Compute Artifact Mask` button. Translucent white artifact masks will appear over regions of tissue that the model flags as putative artifacts. When the `auto` box is selected, the model is run using a reasonable default sensitivity parameter. However, the sensitivity of the algorithm can be adjusted manually by changing the value in the spinbox labeled `Sensitivity`. Each time the algorithm is run on a given channel it adds two layers to the `layers list`. The first layer shows the seed points of artifacts flagged by the algorithm; these are not visible in the image by default, but can be visualized by toggling the eye icon in the `Artifact Seeds` layer. The second layer shows the artifact masks themselves. Individual seed points can be modified or removed from a given channel by highlighting the corresponding `Artifact Seeds` layer, selecting the `arrow icon` in the `layers control` dock to enable point selection mode, then pressing and holding the mouse button and dragging over the seed point to highlight it. Once highlighted, users can fine-tune the artifact mask associated with the highlighted seed by changing the `Tolerance` value in the '`Fine-tuning` widget at the right of the Napari viewer or deleting the seed point by clicking the `x` button in the `layer controls` dock. These AAD tailoring features are designed to give users flexibility over automated artifact masks without having to re-run the AAD algorithm.

### Workflow
After all ROI curations have been made on a given sample, users can move to the next sample in the batch series by clicking the `Apply ROI(s) and Move to Next Sample` button at the top right of the Napari window. If no ROIs are drawn for a given sample, all cells in that tissue will be carried forward into downstream analysis. Alternatively, users can jump between samples by entering their name in the `Sample Name` field at the right of the Napri viewer to add, delete, or modify manual or automated ROIs of previously analyzed tissues or refer to tissues further in the tissue batch series when considering ROI selections for another tissue. ROIs can be added, removed, or modified at any time by re-running the `selectROIs` module.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `delintMode` | False | (bool) Whether to drop (True; negative selection) or retain (False; positive selection) cells selected by ROIs. |
| `showAbChannels` | True | (bool) Whether to show all immunomarker channels (True) when Napari is open (may be memory limiting) or show only cycle 1 DNA (False). |
| `samplesForROISelection` | [ ] | (list of strs) Sample names for ROI selection specified according to the first elements of [sampleMetadata]({{ site.baseurl }}/workflow/input#general-configurations) configuration.
| `autoArtifactDetection` | True | (bool) Whether to display tools for automated artifact detection in Napari window. |
| `artifactDetectionMethod` | "classical" | (str) Algorithm used for automated artifact detection (current option: "classical"). Multi-layer perceptron method ("MLP") currently under development.