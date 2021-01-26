# Input Directory Structure
CyLinter uses standard output from the [mcmicro](https://github.com/labsyspharm/mcmicro) multiplex image preprocessing pipeline. A CyLinter support script (`prep`) programmatically transfers this output into a target destination file then organizes the new input directory for CyLinter analysis.

``` bash
# Transfer and organize mcmicro output for CyLinter analysis
prep <mcmicro_output_dir> <cylinter_input_dir>
```

``` bash
# Data may be transferred from the o2 cluster using an ssh key
prep <userID>@transfer.rc.hms.harvard.edu:</path/to/top-level/mcmicro/output <cylinter_input_dir>
```

For TMA data, pass the "-t" flag prior to source and destination paths. This is required, as mcmicro stores TMA data differently that whole tissue data.  

``` bash
prep -t <mcmicro_output_dir> <cylinter_input_dir>
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

* The CyLinter input directory should be named according to the particular experiment.
* `markers.csv` contains metadata about immunomarkers used in the experiment.
* Tabular data (i.e. csv files containing single-cell data) are located in the `csv/` subdirectory.
* Preprocessed multiplex imaging data (i.e. ome.tif files) are located in the `tif/` subdirectory.
