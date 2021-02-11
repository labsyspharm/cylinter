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
CyLinter presents users with a series of graphical-user-interface (GUI) windows for image inspection, parameter selection, and data visualization in the filtering of image-derived single-cell data. The program proceeds by passing a dataframe of single-cell data (cells x features) containing the combined data from all tissues evaluated in the particular experiment. QC modules are run in the following order:

1. `getSingleCellData`: Combines single-cell data from tissue sections or TMA cores into a single dataframe. This module is automated.

2. `setContrast`: Adjust upper and lower image contrast settings. Visualize immunomarkers by clicking on their respective button at the left side of the Napari viewer. Use the contrast limit slider in the upper-left corner of the Napari viewer to increase image gain by sliding the right side of the slider left or remove background signal intensities by sliding the left side of the slider right. Once contrast settings have been adjusted for all immunomarkers, close the Napari viewer to save the image settings and advance to the next module.

3. `selectROIs`: Select tissue regions-of-interest (ROIs). Draw one or more gates using a mouse or track pad after clicking the triangle icon in the upper-left of the Napari viewer. Set `delint_mode` to `true` in `<cylinter_input_dir>/config.yml` for negative selection. In this case, single-cell data corresponding to selected regions will be excluded from downstream analysis. Close the Napari viewer after ROI(s) for a given tissue have been drawn to save the ROIs and advance to the next tissue in the experiment. To select all data points for a given tissue, simply close the Napari viewer without drawing any gates. ROIs are stored as key:value pairs of sample names and lists of polygon vertices in the `<cylinter_output_dir>/ROI` subdirectory. The ROI selection step will be skipped in future runs of the module and existing ROIs will be applied. Remove `<cylinter_output_dir>/ROI/polygon_dict.pkl` to redefine ROIs.

4. `dnaIntensityCutoff`: Assign upper and lower cutoffs on Hoechst or DAPI signal intensity. This module is designed to remove dim or out-of-focus cells and har-to-segment cells saturated with nuclear counterstain. Use a mouse or track pad to drag the sliders for upper and lower mean intensity cutoffs. Points falling between the upper and lower bounds will be passed to downstream modules. Visualize the selected cells in the context of a particular tissue by typing a tissue name in the provided text box and pressing return. A Napari viewer will show cycle 1 DNA for the selected image with scatter points colored by mean DNA signal intensity for cells falling between the upper and lower bounds. Ensure that at least one Napari viewer remains open while inspecting images, closing all the viewers will cause the current threshold settings to be saved and the program to proceed to the next module.

5. `dnaAreaCutoff`: Assign upper and lower cutoffs on DNA area (nuclei size). This module is designed to remove under- and over- segmented cells. See the description for the `dnaIntensityCutoff` module for usage.

6. `crossCycleCorrelation`: Filter cells whose DNA signal intensity is not correlated across all imaging cycles. Users are presented with histograms showing log10 ratio of DNA signal intensity from cycle 1 compared to all other cycles (columns) for each tissue in the analysis (rows): log10(cycle1/cycleN). Users then enter a y-axis (cell count) cutoff value into the text box window to be applied to all tissues and the name of a tissue of interest separated by a comma. Cycle 1 DNA for the chosen image, along with scatter points of nuclear centroids colorized by their log10(cycleN/cycleN+1) ratio will be shown. Cells whose signal intensity is higher in cycle N will be increasingly pink; cells whose signal intensity is higher in cycle N+1 will be increasingly green. This step is currently for data exploration only and does not filter the data. Close the Napari viewer and the histogram plots will come up again, this time with a black horizontal line at the selected count cutoff. Chose a different cutoff if unsatisfactory and repeat the process as needed. The final cutoff will be saved after closing the Napari viewer, histogram plots, and the text box window.

7. `log10transform`: Log-transform and rescale data per channel. This module is fully automated, but may become configurable in future release.

8. `pruneOutliers`: Assign upper and lower percentile cutoffs on immunomarker signal intensity. Users will be presented with hexbins or scatter plots of immunomarker signal intensity (x-axes) vs. nuclear area (y-axes) (see the `hexbins` parameter in the configuration file for details) for all tissues in the analysis. After inspecting the plots, users will type lower and upper percentile cutoffs (scale: 0.0-100.0) separated by a comma into the text box window and press return. Pruned plots will come up for inspection; enter new cutoffs and repeat as necessary. Once the optimal cutoffs have been identified, close the text box window and a copy of the data will be filtered with the selected cutoffs. Repeat for all immunomarker channels. If the text box window is closed with nothing entered into it, lower and upper cutoff values of 0.1 and 99.9 will be applied by default. Cutoff values are stored as key:value pairs (i.e. sampleName: (upperCutoff, lowerCutoff)) in `<cylinter_output_dir>/hexbins`. Once cutoffs for all channels have been made, the program will loop through the cutoff dictionary and drop cells from the actual dataframe in the same order in which they were curated. Cutoff selection will be skipped and existing cutoffs will be automatically applied in further runs of this module until `<cylinter_output_dir>/hexbins/hexbin_dict.pkl` is removed.

9. `performPCA`: Perform PCA on mean immunomarker signal intensities from whole tissue sections or TMA cores. This module is automated and configurable; see `config.yml` for parameters and their descriptions.

10. `performTSNE`: Perform TSNE embedding and density-based clustering (HDBSCAN) on single-cell data from whole tissue sections or TMA cores. This module is configurable; see `config.yml` for parameters and their descriptions. After the t-SNE embedding has been computed, users will be presented with a text box for entering a minimum cluster size. This parameter significantly influences the resulting clustering solution. A range of minimum cluster size values may be passed (`startValue-endValue`) and the number of clusters per minimum cluster size will be printed to the terminal window. This feature is designed to rapidly identify stable clustering solutions. Once the optimal minimum cluster size has been identified, close the text box window and the current clustering IDs will be added to the dataframe.

11. `getClustermap`: Compute and plot clusters vs. immunomarkers clustermap. This module is automated. Simply close the clustermap after visual inspection.

12. `lassoClusters`: Lasso cells of interest from the embedding for further analysis. Using a mouse or track pad, lasso cells of interest in the t-SNE embedding to print their 3 most-highly expressed immunomarkers to the terminal window. Currently, this module doesn't filter the data, but may become useful for performing off-hand, exploratory analysis in future releases.

13. `curateThumbnails`: Curate thumbnail images of cells from each cluster for visualization. This module is automated but configurable; see `config.yml` for parameters and their descriptions. Thumbnail images are saved to the `<cylinter_output_dir>/thumbnails` subdirectory.

Parquet files containing the dataframe filtrate returned by each module can be found in the `<cylinter_output_dir>/checkpoints` subdirectory and be used in downstream analysis.
