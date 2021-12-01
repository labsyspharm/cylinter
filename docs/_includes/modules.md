# Module list

| Name | Purpose | Description/YAML Configs |
| :-- | :-- | :-- |
| `aggregateData` | Combines feature tables | [Details]({{ site.baseurl }}modules/aggregateData) |
| `selectROIs` | User-defined ROI selection | [Details]({{ site.baseurl }}modules/selectROIs) |
| `intensityFilter` | Filter on DNA intensity | [Details]({{ site.baseurl }}modules/intensityFilter) |
| `areaFilter` | Filter on cell segmentation area | [Details]({{ site.baseurl }}modules/areaFilter) |
| `cycleCorrelation` | Filter on detached cells | [Details]({{ site.baseurl }}modules/cycleCorrelation) |
| `logTransform` | Log10-transformation of immunomarker signals | [Details]({{ site.baseurl }}modules/logTransform)
| `pruneOutliers` | Filter on channel outliers | [Details]({{ site.baseurl }}modules/pruneOutliers) |
| `metaQC` |  Unsupervised QC status reclassification | [Details]({{ site.baseurl }}modules/metaQC)
| `PCA` | Principle component analysis | [Details]({{ site.baseurl }}modules/PCA)
| `clustering` | Identify cell states | [Details]({{ site.baseurl }}modules/clustering)
| `clustermap` | Visualize cell state protein expression | [Details]({{ site.baseurl }}modules/clustermap)
| `setContrast` | Adjust image contrast settings | [Details]({{ site.baseurl }}modules/setContrast)
| `curateThumbnails` | Visualize examples of cell states | [Details]({{ site.baseurl }}modules/curateThumbnails)
| `frequencyStats` | Compute cell state frequency statistics | [Details]({{ site.baseurl }}modules/frequencyStats) |

<br/>

# Suggest a module
The CyLinter team is collaborating with the NCI to run hackathons and challenges to improve existing modules and increasingly integrate automation into the pipeline. Modules are also being added to CyLinter incrementally by a diverse developer community seeded by the NCI [Human Tissue Atlas Network](https://humantumoratlas.org/). See what [modules]({{ site.baseurl }}modules/index) we are currently using. Module suggestions can be made by posting to [https://forum.image.sc/](https://forum.image.sc/) and tagging your post with the `cylinter` tag.
