# Transferring and organizing input data
CyLinter uses a support console script to programmatically format standard output files from the [MCMICRO](https://github.com/labsyspharm/mcmicro) image-processing pipeline as standard input for CyLinter. CyLinter supports analysis of whole tissue sections and tissue microarray (TMA) data pre-processed by MCMICRO.

To organize **WHOLE TISSUE DATA** as CyLinter input, run the following command:

``` bash
prep <mcmicro_output_path> <cylinter_input_path>  # < > indicates a variable.
```

To organize **TMA DATA** as CyLinter input, run the following command:

``` bash
# prepend a "-t" flag before the source and destination directory file paths.
prep -t <mcmicro_output_path> <cylinter_input_path>
```

SSH keys can be used to transfer MCMICRO output from the o2 compute cluster at HMS (requires access permission).

``` bash
prep <userID>@transfer.rc.hms.harvard.edu:<mcmicro_output_path> <cylinter_input_path>
```

CyLinter input directories are organized as follows:

``` bash
<cylinter_input_path>
├── config.yml
├── csv
│   ├── unmicst-<tissue/core1>.csv
│   └── unmicst-<tissue/core2>.csv
├── markers.csv
├── seg
│   ├── <tissue/core1>.tif (standard TIFF format)
│   └── <tissue/core2>.tif
└── tif
    ├── <tissue/core1>.tif (OME-TIFF file format)
    └── <tissue/core2>.tif
```

* In the above example CyLinter input directory, `config.yml` is a template configuration file for specifying program configurations and module parameters. The template is pre-formatted for use with an exemplar dataset (click "Exemplar data" tab at left for details).

* The `markers.csv` contains immunomarker metadata for the particular experiment. The file should be formatted as follows:

```
channel_number,cycle_number,marker_name
1,1,<DNA1>
2,1,<abx1>
3,1,<abx2>
4,1,<abx3>
5,2,<DNA2>
6,2,<abx4>
7,2,<abx5>
8,2,<abx6>
.
.
.
```
* Additional metadata columns are permissible, but are not currently used by CyLinter.
