# Input Directory Structure
CyLinter takes standard output from the [mcmicro](https://github.com/labsyspharm/mcmicro) multiplex image preprocessing pipeline as input. CyLinter comes with a support console script called `prep` that prepares mcmicro output as CyLinter input.

``` bash
# Transfer mcmicro output to CyLinter input directory and organize files for analysis
prep <mcmicro_output_dir> <cylinter_input_dir>

# mcmicro stores TMA data differently that whole tissue data. A "-t" flag must be passed before the source and destination paths to indicate that the mcmicro output are TMA data.
prep -t <mcmicro_output_dir> <cylinter_input_dir>
```

``` bash
# SSH keys maybe used to transfer mcmicro output from remote sources such as HMS's o2 compute cluster.
prep <userID>@transfer.rc.hms.harvard.edu:</path/to/top-level/mcmicro/output <cylinter_input_dir>
```

CyLinter input directories take the following form:

``` bash
<cylinter_input_dir>
├── config.yml
├── csv
│   ├── unmicst-<sample/core1>.csv
│   └── unmicst-<sample/core2>.csv
├── markers.csv
└── tif
    ├── <sample/core1>.ome.tif
    └── <sample/core2>.ome.tif
```

Notes:

* The CyLinter input directories should be named according to their particular experiments.
* The `markers.csv` file contains metadata about immunomarkers used in the experiment and must be present for analysis.
* Preprocessed multiplex imaging files (i.e. ome.tif files) are located in the `tif/` subdirectory.
* Corresponding tabular data are located in the `csv/` subdirectory.
