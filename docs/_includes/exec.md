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
1. `getSingleCellData:` combine whole tissue or TMA core data into a single dataframe.
2. `setContrast:` adjust upper and lower image contrast bounds.
3. `selectROIs:` gate regions of interest. Positive (inclusive) and negative (exclusive) selection options.
4. `dnaIntensityCutoff:` assign upper and lower cutoffs on DNA counter stain (Hoechst, DAPI) intensity.
5. `dnaAreaCutoff:` assign upper and lower cutoffs on DNA area (nucleus size).
6. `crossCycleCorrelation:` filter poorly-correlated cells across imaging cycles according to DNA intensity.
7. `log10transform:` log-transform and rescale data (per channel).
8. `pruneOutliers:` assign upper and lower percentile cutoffs on immunomarker signal intensity.  
9. `performPCA:` perform PCA on whole tissue or TMA core mean immunomarker intensities.
10. `performTSNE:` perform TSNE embedding and density-based clustering (HDBSCAN) on cleaned single-cell data.
11. `getClustermap:` compute and plot clusters vs. immunomarkers clustermap
12. `lassoClusters:` lasso cells of interest from the embedding for further analysis.
13. `curateThumbnails:` curate thumbnail images of cells from each cluster for visual inspection.
