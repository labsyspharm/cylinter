# Input Directory Structure
CyLinter uses standard output files from the [mcmicro](https://github.com/labsyspharm/mcmicro) multiplex image preprocessing pipeline as input into the program.

mcmicro is typically deployed on compute clusters (e.g. Harvard Medical School's o2 cluster), but can also be run locally. Transfer and reorganize mcmicro output data files as CyLinter input by running the following command:

``` bash
prep <source_path> <destination_path>
```
To transfer mcmicro output from o2 <source_path> my be written as <userID>@transfer.rc.hms.harvard.edu:<path/to/top-level/mcmicro/output/dir>. The destination path is chosen by the user and serves as the CyLinter input directory.

For TMA data, an additional "-t" flag must be passed prior to the source and destination paths, as mcmicro stores TMA data differently that whole tissue data.  

``` bash
prep -t <source_path> <destination_path>.
```

Resulting CyLinter input directories will take the following form:

``` bash
<input_dir>  # eqivalent to <destination_path>
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
