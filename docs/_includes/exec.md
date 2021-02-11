# Configuration settings

Update the template configuration file (`config.yml`) in <cylinter_input_dir> with experiment-specific metadata.

``` yaml
in_dir: "<cylinter_input_dir>"
out_dir: "<cylinter_output_dir>"  # location is arbitrary
random_sample_size: 1.0  # floating point (0.0-1.0); 1.0 is full dataset
mask_object: cellMask  # cellMask, nucleiMask, cellRingMask; see csv table column headers
sample_metadata:
  unmicst-<sampleString>: ["<fullConditionString>", "<abbrConditionString>", <replicateInteger>]
samples_to_exclude: ["<sampleString1>", "<sampleString2>", ...]
markers_to_exclude: ["<markerString1>", "<markerString2>", ...]
.
.
.
```

# Pipeline execution

``` bash
# Activate the CyLinter virtual environment.
source ~/cylinter/bin/activate

# Perform CyLinter analysis by passing the YAML configuration file.
cylinter --module (optional) <cylinter_input_dir>/config.yml
```

* By default, the pipeline starts at the first module: `getSingleCellData`. Passing the name of a specific module with `--module <module_name>` causes the program to read the cached parquet file returned by the preceding module.

# User guide
CyLinter users are presented a series of graphical-user-interface (GUI) windows for image inspection, parameter selection, and data visualization to filter of image-derived single-cell data. The program proceeds by passing a dataframe of single-cell data (cells x feature) containing the combined data from all tissue samples evaluated in a particular experiment. QC modules are run in the following ordered series:

1. `getSingleCellData`: Combines tabular data from whole tissue sections or TMA cores into one dataframe. This module is automated.

2. `setContrast`: Adjust upper and lower image contrast bounds. Visualize immunomarker channels by clicking on their corresponding buttons at the left-hand side of the Napari viewer. Use the contrast limit slider at the upper-left of the Napari viewer to increase image gain by sliding the right side of the slider to the left; remove background signal intensities by sliding the left side of the slider to the right. Once thresholds for all channels have been set, close the Napari viewer by clicking on the button at the upper-left of the window.

3. `selectROIs`: Select tissue regions-of-interest (ROIs). Positive and negative selection modes are available. Set configuration parameter `delint_mode` to `false` to positive selection: single-cell data corresponding to selected regions will be included in downstream analysis. Set `delint_mode` to `true` for negative selection: single-cell data corresponding to selected regions will be excluded from downstream analysis. Draw gates using a mouse or track pad after clicking the triangle icon in the upper-left of the Napari viewer. Multiple ROIs are permissible per tissue section. Click the up-arrow key to begin drawing a new ROI. Close the Napari viewer after ROI(s) for a given tissue have been drawn to save the ROIs and advance to the next tissue in the experiment. Saved ROIs are stored as key:value pairs (sample_name: [polygon vertices]) in the `<cylinter_output_dir>/ROI` subdirectory. The ROI selection step will be skipped in all future runs of the `selectROIs` module and the existing ROIs will be applied. Remove `<cylinter_output_dir>/ROI/polygon_dict.pkl` to redefine ROIs.

4. `dnaIntensityCutoff`: Assign upper and lower cutoffs on Hoechst or DAPI signal intensity to remove out-of-focus cells and saturated cells prone to nuclei mis-segmentation. Use a mouse or track pad to drag slider bars for upper and lower cutoffs on per cell mean DNA signal intensity. Data points falling between upper and lower bounds are carried forward to downstream modules. Visualize included data by entering the name of a tissue sample of interest and select return. A Napari viewer will come up showing the cycle 1 DNA channel with scatter points superimposed on cell centroids falling between upper and lower bounds. Scatter points are colored according to mean DNA signal intensity. Enter a new sample name then press return to inspect additional samples. Note: ensure that at least one Napari window is open while inspecting images. Closing all Napari windows will cause the current threshold settings to be saved and the program to move to the next module.

5. `dnaAreaCutoff`: Assign upper and lower cutoffs on DNA area (nuclei size) to remove under- and over- segmented cells. See module 4 (`dnaIntensityCutoff`) for usage.

6. `crossCycleCorrelation`: Filter cells whose DNA signal intensity is not correlated across all imaging cycles. First inspect the log10(cycle1/cycleN) histograms of DNA signal intensities (columns) that are automatically presented for all tissues in the experiment (rows). Decide on a common y-axis cutoff value. Close the window and enter the chosen value in the provided text window. Enter the name of a tissue of interest and evaluate the DNA images and cell scatter points colorized according to log10(cycleN/cycleN+1) ratio using a diverging color map. Signal intensities higher in cycle N (pink); signal intensities higher in in cycle N+1 (green). After closing the Napari window the histogram plots will come up again, this time with a black horizontal line at the selected count cutoff. Chose a different cutoff and repeat the process if needed. Apply the final cutoff by closing the Napari and text box windows.

7. `log10transform`: Log-transform and rescale data per channel. This module is automated; may become configurable in a future release.

8. `pruneOutliers`: Assign upper and lower percentile cutoffs on immunomarker signal intensity. Users are presented with hexbins or scatter plots of immunomarker signal intensity (x-axes) vs. nuclear area (y-axes) for every tissue in the experiment (see `hexbins` config parameter). After inspecting the plots, type reasonable lower and upper percentile cutoffs (scale: 0.0-100.0) separated by a comma into the text box window and press the return key. New plots (post-pruning) will come up for inspection; enter new cutoffs and repeat if necessary. Once the optimal cutoffs have been evaluated (assumed to be the last applied) and the text box window is closed, a copy of the log10transform.parquet file will be filtered with the selected cutoffs. If nothing is entered into the text box window, and the window is closed, lower and upper cutoff values of 0.1 and 99.9 will be automatically applied. Subsequent channel cutoffs will be applied to the dataframe copy. Cutoff values are stored as key:value pairs (i.e. sampleName: (upperCutoff, lowerCutoff)) in `<cylinter_output_dir>/hexbins`. Once cutoffs for all channels have been curated, the program will loop through the cutoff dictionary and drop cells after applying the cutoffs in the same order in which they were curated. The cutoff selection step will be skipped in future runs of this module and existing cutoffs will be applied. Remove `<cylinter_output_dir>/hexbins/hexbin_dict.pkl` to redefine percentile cutoffs.

9. `performPCA`: Perform PCA on mean immunomarker signal intensities from whole tissue sections or TMA cores. This module is automated and configurable; for parameters and their descriptions, see `config.yml`.

10. `performTSNE`: Perform TSNE embedding and density-based clustering (HDBSCAN) on single-cell data from whole tissue sections or TMA cores. This module is configurable; for parameters and their descriptions, see `config.yml`. After the t-SNE embedding has been computed, users will be presented with a text box for entering a minimum cluster size to visualize the resulting clustering; a range of minimum cluster sizes can also be passed: `startValue-endValue`. In this case, the number of clusters identified per minimum cluster size is printed to the terminal window for identifying a stable clustering solution. Once the optimal minimum cluster size has been identified, close the text box window, and the cluster IDs will be added to the current dataframe.

11. `getClustermap`: Compute and plot clusters vs. immunomarkers clustermap. This module is automated.

12. `lassoClusters`: Lasso cells of interest from the embedding for further analysis. Using a mouse or track pad, lasso cells of interest in the t-SNE embedding to print their 3 highest expressed immunomarkers to the terminal window. The module doesn't filter any data, but may become useful for performing off-hand, exploratory analysis in future releases.

13. `curateThumbnails`: Curate thumbnail images of cells from each cluster for visualization. This module is automated but configurable; see `config.yml`. Thumbnail images are saved to the `<cylinter_output_dir>/thumbnails` subdirectory.

Users can peruse intermediate data output files in `<cylinter_output_dir>` at any point in the pipeline. Parquet files containing the dataframe filtrate returned by each module can be found in the `<cylinter_output_dir>/checkpoints` subdirectory and be used in further analysis.
