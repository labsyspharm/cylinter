# Pipeline execution

After installation, Cylinter can be executed as a console script at the command using the following incantation:

``` bash
cylinter --module (optional) path/to/input/dir
```

By default, the pipeline starts from the "getSingleCellData" module. Use the `--module` flag to execute from any module in the pipeline instead.

By default CyLinter writes module-specific intermediate parquet files to a `checkpoints/` subdirectory inside whatever parent output directory specified in the config.yml file. (See the next section for more information on pipeline configuration settings.)


### Specifying configuration settings

Run-specific parameters are passed to a YAML configuration (`input/config.yml`), which is provided to CyLinter at the point of execution.

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
