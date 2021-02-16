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
CyLinter presents users with a series of graphical-user-interface (GUI) windows for image inspection, parameter selection, and data visualization in order to filter image-derived, single-cell data corrupted by microscopy artifacts. The program passes a dataframe of single-cell data from all tissues in the analysis through the following series of QC modules:

1. `getSingleCellData`: Combines single-cell data from tissue sections or TMA cores into a single dataframe. This module is automated.

2. `setContrast`: Adjust upper and lower image contrast settings. This module is for image visualization purposes only and will not affect the single-cell data. Visualize immunomarkers by clicking on their named buttons at the left of the Napari viewer. Use the contrast limit slider in the upper-left corner of the viewer to increase image gain (by sliding the right side of the slider left), or remove background signal intensities (by sliding the left side of the slider right). Once contrast settings have been adjusted for all immunomarkers, close the Napari viewer to save the image settings and advance to the next module. Image contrast settings are reflected in the images saved by the `curateThumbnails` module below. Current contrast settings can be updated by re-running the module.

3. `selectROIs`: Gate on tissue regions-of-interest (ROIs). Draw one or more polygon gates using a mouse or track pad after clicking the triangle icon in the upper-left of the Napari viewer. Set `delint_mode` to `true` in `<cylinter_input_dir>/config.yml` for negative selection. This will cause single-cell data corresponding to the selected regions to be excluded from downstream analysis. Close the Napari viewer after ROI(s) for a given tissue have been drawn to save the ROIs and advance to the next tissue in the analysis. To select all cells from a given tissue, simply close the Napari viewer without drawing any gates. ROIs are stored as a a dictionary of key:value pairs in the `<cylinter_output_dir>/ROI` subdirectory (i.e. sample: [polygon vertices]). The ROI selection step will be skipped and existing ROIs will be applied in further runs of the module. Remove `<cylinter_output_dir>/ROI/polygon_dict.pkl` to redefine ROIs.

4. `dnaIntensityCutoff`: Assign upper and lower cutoffs on Hoechst or DAPI signal intensity. This module is designed to remove dim or out-of-focus cells and hard-to-segment cells saturated with nuclear counterstain. Users are presented with a histogram of integrated DNA signal intensity for individual cells among all tissues in the analysis. Use a mouse or track pad to drag the sliders to assign upper and lower cutoffs. Cells falling between the lower and upper bounds are passed to downstream modules. Visualize where selected cells reside in a particular tissue by typing the tissue name into the text box window and pressing return. A Napari viewer will then show the first DNA channel for the specified tissue. Scatter points appear over cell centroids included by the current cutoffs. Scatter points are colored according to the integrated DNA signal intensity for individual cells scaled between the lower and upper cutoffs. Ensure that at least one Napari viewer remains open while inspecting images, closing all the viewers will cause the current threshold settings to be saved and the program to proceed to the next module.

5. `dnaAreaCutoff`: Assign upper and lower cutoffs on DNA area (nuclei size). This module is helpful for removing under- and over- segmented cells. Module usage is the same as module 4: `dnaIntensityCutoff`.

6. `crossCycleCorrelation`: Filter cells whose DNA signal intensity is not correlated across all imaging cycles. This  module guards against false immunomaker signatures incurred by tissue movement and. Users are presented with histograms of the log10 ratio of DNA signal intensity from cycle 1 relative to all additional cycles (columns) for each tissue in the analysis (rows). Enter a y-axis value into the text box window to apply a count cutoff to all tissues. Enter the name of a tissue of interest separated by a comma to visualize how DNA signal intensity changes in different areas of the tissue across imaging cycles. Scatter points of nuclear centroids colorized by their log10(cycleN/cycleN+1) ratio will be shown. Cells whose signal intensity is higher in cycle N will be increasingly pink; cells whose signal intensity is higher in cycle N+1 will be increasingly green. This step is currently for data exploration only and does not filter the single-cell data. Closing the Napari viewer will cause the histogram plots to come up again, this time with a black horizontal line at the selected count cutoff. Confirm the cutoff is appropriate for all tissues, or chose a different cutoff and repeat the process if needed. The final cutoff will be saved after closing all windows.

7. `log10transform`: Log-transform and rescale data per channel. This module is fully automated, but may become configurable in future release.

8. `pruneOutliers`: Assign upper and lower percentile cutoffs on immunomarker signal intensities. Users are presented with hexbins or scatter plots of immunomarker signal intensity (x-axes) vs. nuclear area (y-axes) for all tissues in the analysis (see the `hexbins` parameter in the configuration file for details). After inspecting the plots, type lower and upper percentile cutoffs (scale: 0.0-100.0) separated by a comma into the text box window and press return. Pruned plots will come up for inspection; enter new cutoffs and repeat as necessary. Once the optimal cutoffs have been identified, close the text box window and the data will be filtered with the selected cutoffs. Repeat for all immunomarkers in the analysis. If the text box window is closed with nothing entered into it, lower and upper cutoff values of 0.1 and 99.9 will be applied by default. Cutoff values are stored as key:value pairs (i.e. sample: (upperCutoff, lowerCutoff)) in `<cylinter_output_dir>/hexbins`. Once cutoffs for all channels have been made, the program will loop through the cutoff dictionary and drop cells from the dataframe in the same order in which they were curated. Cutoff selection is skipped and existing cutoffs are automatically applied in further runs of this module. Remove `<cylinter_output_dir>/hexbins/hexbin_dict.pkl` to chose different cutoffs.

9. `performPCA`: Perform PCA on mean immunomarker signal intensities from whole tissue sections or TMA cores. This module is automated and configurable; see `config.yml` for parameters and their descriptions.

10. `performTSNE`: Perform TSNE embedding and density-based clustering (HDBSCAN) on single-cell data from whole tissue sections or TMA cores. This module is configurable; see `config.yml` for parameters and their descriptions. After the t-SNE embedding has been computed, users will be presented with a text box for entering a `minimum cluster size`. This parameter significantly influences the resulting clustering solution. A range of minimum cluster size values may be passed as `startValue-endValue` and the number of clusters per minimum cluster size value will be printed to the terminal window. This feature aids in identify a stable clustering solution. Once the optimal minimum cluster size has been identified, add ".save" to the end of the value (e.g. 45.save) and press return. This will add the cluster IDs for the current clustering solution to the dataframe returned by the module.

11. `getClustermap`: Compute and plot clusters vs. immunomarkers clustermap. This module is automated. Close the clustermap after visual inspection.

12. `lassoClusters`: Lasso cells of interest from the embedding for further analysis. Using a mouse or track pad, lasso cells of interest in the t-SNE embedding to print their 3 most-highly expressed immunomarkers to the terminal window. Currently, this module does not filter the single-cell data, but may become useful for performing off-hand, exploratory analysis in future releases.

13. `curateThumbnails`: Curate thumbnail images of cells from each cluster for visualization. This module is automated but configurable; see `config.yml` for parameters and their descriptions. Thumbnail images are saved to the `<cylinter_output_dir>/thumbnails` subdirectory.

Parquet files containing the dataframe filtrate returned by each module are found in the `<cylinter_output_dir>/checkpoints` subdirectory and available for use in downstream analyses.
