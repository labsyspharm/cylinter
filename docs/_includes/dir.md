# Transferring and organizing input data
CyLinter uses a support console script to programmatically format standard output files from the [MCMICRO](https://github.com/labsyspharm/mcmicro) image-processing pipeline as standard input for CyLinter. CyLinter supports analysis of whole tissue sections and tissue microarray (TMA) data pre-processed by MCMICRO. To organize **WHOLE TISSUE DATA** as CyLinter input, run the following command:

``` bash
prep <mcmicro_output_dir> <cylinter_input_dir>  # < > indicates a variable.
```

To organize **TMA DATA** as CyLinter input, run the following command:

``` bash
# prepend a "-t" flag before the source and destination directory file paths.
prep -t <mcmicro_output_dir> <cylinter_input_dir>
```

SSH keys can be used to transfer MCMICRO output from the o2 compute cluster at HMS (requires access permission).

``` bash
prep <userID>@transfer.rc.hms.harvard.edu:<mcmicro_output_dir> <cylinter_input_dir>
```

CyLinter input directories are organized as follows:

``` bash
<cylinter_input_dir>
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

* `config.yml` is a template configuration file for specifying program configurations and module parameters.
The template is pre-formatted for use with an exemplar dataset (`emit22_demo`, see "Exemplar data" tab).

* `markers.csv` contains immunomarker metadata for the particular experiment. The file must be formatted as follows:

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
* Additional metadata columns are permissible, but not used by the current version of CyLinter.
