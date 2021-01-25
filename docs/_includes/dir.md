# Input Directory Structure
CyLinter uses standard output files from the [mcmicro](https://github.com/labsyspharm/mcmicro) multiplex image preprocessing pipeline as input.

Before running CyLinter, transfer and organize mcmicro output data files (typically on the o2 cluster) as CyLinter input using the following command:

``` bash
prep <source_dir_path> <dest_dir_path>
```

In the case of TMA data, an additional "-t" flag must be passed prior to the source and destination file path variables, as mcmicro handles TMA data differently that whole tissue data.  

``` bash
prep -t <source_dir_path> <dest_dir_path>.
```

The resulting input directory will then be organized as follows:

```
</input/dir>
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

* The basename of the input directory (i.e. </input/dir>) is the name of the experiment.
* The file `markers.csv` contains metadata about experimental immunomarkers.
* Tabular data (i.e. `csv files`) are located in `csv/` subdirectories.
* Preprocessed multiplex imaging data (i.e. `ome.tif files`) are located in `tif/` subdirectories.
