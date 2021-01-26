# Configuration settings

Before running CyLinter, specify experimental metadata by editing the YAML configuration file in the corresponding CyLinter input directory.

``` yaml
in_dir: /cylinter/input/dir
out_dir: /arbitrary/output/dir
random_sample_size: 1.0 # floating point (0.0-1.0; 1.0 = full dataset)
mask_object: cellMask  # cellMask, nucleiMask, cellRingMask; determined by mcmicro run parameters (see csv column headers)
sample_metadata:
  unmicst-<sampleString>: [<fullConditionString>, <abbrConditionString>, <replicateInteger>]
samples_to_exclude: [<sampleString1>, <sampleString2>, ...]
markers_to_exclude: [<markerString1>, <markerString2>, ...]
```

# Pipeline execution

Run CyLinter by pointing the tool at a formatted input directory using the following command:

``` bash
source $HOME/cylinter/bin/activate  # activates the CyLinter virtual environment
cylinter --module (optional) <input_dir>/config.yml
```

By default, the pipeline starts at the first module (`getSingleCellData`). Passing the name of a specific module  with `--module <module_name>` causes the pipeline to start from the specified module. The program will look for cached parquet files from the preceding module located in the `checkpoints` subdirectory of the specified output directory.

  ``` bash
  cylinter --module selectROIs </input/dir>/config.yml
  ```
