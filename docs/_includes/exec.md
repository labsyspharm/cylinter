# Program execution

The CyLinter quality control pipeline proceeds with users interacting with the program through a series of graphical windows for image inspection, cutoff selection, and data visualization. A YAML configuration file (`<cylinter_input_path>/config.yml`) is used to define experiment-specific configuration settings.

The program is executed by typing the `cylinter` command followed by the path to a configuration file:  

``` bash
# Activate virtual environment
source ~/cylinter/bin/activate

# Run pipeline beginning from the first module:
cylinter <cylinter_input_path>/config.yml
```

CyLinter saves partially-redacted feature tables as `.parquet` files in `cylinter_output_path/checkpoints/` to allow for dynamic restarts. In addition to passing a configuration file path to `cylinter`, specifying the name of a particular module using the `--module` flag allows users to perform analysis starting from a particular QC module:

``` bash
cylinter --module <module_name> <cylinter_input_path>/config.yml
```
