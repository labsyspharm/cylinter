# Transferring and organizing input data
CyLinter uses a support console script for programmatic formatting of standard output files from the [MCMICRO](https://github.com/labsyspharm/mcmicro) image-processing pipeline as CyLinter input. Whole tissue and tissue microarray (TMA) data are processed differently by MCMICRO. Organize MCMICRO output from **WHOLE TISSUE DATA** into a CyLinter input run the following command:

``` bash
prep <mcmicro_output_dir> <cylinter_input_dir>  # < > indicates a variable.
```

Organize MCMICRO output from **TMA DATA** into a CyLinter input prepend a "-t" flag before the source and destination directory file paths:

``` bash
prep -t <mcmicro_output_dir> <cylinter_input_dir>
```

SSH keys may be used to transfer MCMICRO output from the o2 compute cluster at HMS (requires access permission).

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

* `config.yml` is a template CyLinter configuration file pre-formatted for use with CyLinter's example dataset (`emit22_demo`, see exemplar tab).

* `markers.csv` contains experiment-specific immunomarker metadata. The file must be organized as follows:

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
* Additional metadata columns are permissible, but not used by the current version of the algorithm.
