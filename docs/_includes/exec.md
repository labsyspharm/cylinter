# Configuration settings

Specify experimental metadata before running CyLinter by editing the YAML configuration file (`config.yml`) in the corresponding CyLinter input directory.

``` yaml
in_dir: <cylinter_input_dir>
out_dir: <cylinter_output_dir>  # location is arbitrary
random_sample_size: 1.0  # floating point (0.0-1.0); 1.0 is full dataset
mask_object: cellMask  # cellMask, nucleiMask, cellRingMask; see csv table column headers
sample_metadata:
  unmicst-<sampleString>: [<fullConditionString>, <abbrConditionString>, <replicateInteger>]
samples_to_exclude: [<sampleString1>, <sampleString2>, ...]
markers_to_exclude: [<markerString1>, <markerString2>, ...]
```

# Pipeline execution

``` bash
# Activate the CyLinter virtual environment
source $HOME/cylinter/bin/activate

# Pass configuration file to CyLinter for analysis
cylinter --module (optional) <input_dir>/config.yml
```

* By default, the pipeline starts at the first module: `getSingleCellData`. Passing the name of specific modules with `--module <module_name>` causes the program to retrieve the cached parquet file returned by the preceding module.

# Current modules
`getSingleCellData`— combine tissue sample (or TMA core) data into single dataframe.
`setContrast`— adjust upper and lower image contrast settings.
`selectROIs`— gate tissue regions of interest. Cells within region bounds are carried forward into downstream analysis. Alternatively, gate areas to be excluded from downstream analysis.
`dnaIntensityCutoff`— assign upper and lower cutoffs on DNA intensity.
`dnaAreaCutoff`— assign upper and lower cutoffs on DNA area (nucleus size).
`crossCycleCorrelation`— remove cells poorly-correlated in DNA intensity across imaging cycles.
`log10transform`— log-transform and rescale data (per channel).
`pruneOutliers`— assign upper and lower percentile cutoffs on immunomarker signal intensity.  
`performPCA`— perform PCA on tissue samples (or TMA cores) using mean immunomarker intensity values.
`performTSNE`— perform TSNE embedding and density-based clustering (HDBSCAN) on high-quality cells.
`getClustermap`— compute and plot clusters vs. immunomarkers clustermap
`lassoClusters`— select cells of interest from specific clusters for further analysis.
`curateThumbnails`— curate images of cells from specific clusters for validation.
