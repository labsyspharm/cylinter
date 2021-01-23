# Configuration settings

Before running CyLinter, the user must specify parameter values for a series of configurations settings by editing an input directory's YAML configuration file.

``` yaml
in_dir: path/to/input/dir
out_dir: path/to/chosen/output/dir
random_sample_size: 1.0 # floating point (0.0-1.0; 1.0 = full dataset)
mask_object: cellMask  # cellMask, nucleiMask, cellRingMask; indicated by column headers in csv files; determined by mcmicro run parameters.
sample_metadata:
  unmicst-<sampleString>: [<fullConditionString>, <abbrConditionString>, <replicateInteger>]
samples_to_exclude: [<sampleString1>, <sampleString2>, ...]
markers_to_exclude: [<markerString1>, <markerString2>, ...]
```

# Pipeline execution

To run CyLinter, point the tool at a previously-formatted input directory using the following command:

``` bash
cylinter --module (optional) <input_dir_path>
```

* Without passing --module, the pipeline starts by default at the first module `getSingleCellData`. Passing a module name will cause the pipeline to start from the specified module. In that case, the program will look for cached parquet files (compressed csv files) from the immediately preceding module stored in the `checkpoints` subdirectory of the top-level output directory specified in the corresponding input directory's YAML configuration file.
  * For example, running the pipeline with `cylinter --module selectROIs <input_dir_path>` will start the pipeline at the region of interest (ROI) section module and run all subsequent modules until the user terminates the program.
