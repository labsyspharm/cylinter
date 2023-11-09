---
layout: default-cylinter
title: Modules
nav_order: 5
has_children: true
---

# Module list

| Name | Purpose | Description/YAML Configs |
| :-- | :-- | :-- |
| `aggregateData` | Combine feature tables | [Details]({{ site.baseurl }}/modules/aggregateData) |
| `selectROIs` | Define tissue ROIs | [Details]({{ site.baseurl }}/modules/selectROIs) |
| `intensityFilter` | Filter out-of-focus and counterstain oversaturated cells | [Details]({{ site.baseurl }}/modules/intensityFilter) |
| `areaFilter` | Filter over- and under-segmented cells | [Details]({{ site.baseurl }}/modules/areaFilter) |
| `cycleCorrelation` | Filter unstable cells | [Details]({{ site.baseurl }}/modules/cycleCorrelation) |
| `logTransform` | Log10-transform immunomarker signals | [Details]({{ site.baseurl }}/modules/logTransform)
| `pruneOutliers` | Filter channel outliers | [Details]({{ site.baseurl }}/modules/pruneOutliers) |
| `metaQC` |  Reclassify cells according to QC status  | [Details]({{ site.baseurl }}/modules/metaQC)
| `PCA` | Run principle component analysis | [Details]({{ site.baseurl }}/modules/PCA)
| `setContrast` | Adjust image contrast settings | [Details]({{ site.baseurl }}/modules/setContrast)
| `clustering` | Identify cell states | [Details]({{ site.baseurl }}/modules/clustering)
| `clustermap` | Visualize cell state protein expression | [Details]({{ site.baseurl }}/modules/clustermap)

| `frequencyStats` | Compute cluster frequency statistics | [Details]({{ site.baseurl }}/modules/frequencyStats) |
| `curateThumbnails` | Visualize example cells from each cluster | [Details]({{ site.baseurl }}/modules/curateThumbnails)

<br/>

# Suggest a module
The CyLinter team is collaborating with NCI-sponsored consortia (CSBC and PS-ON) to host hackathons to improve and automate existing methods for microscopy quality control like those instantiated by the CyLinter pipeline. CyLinter modules are also being added incrementally by a diverse developer community seeded by the NCI [Human Tissue Atlas Network](https://humantumoratlas.org/). See what modules are currently available [here]({{ site.baseurl }}/modules/index). Module suggestions can be made by posting to [https://forum.image.sc/](https://forum.image.sc/) and tagging your post with the `cylinter` tag.
