# Input Directory Structure
CyLinter uses standard output files from the [mcmicro](https://github.com/labsyspharm/mcmicro) multiplex image preprocessing pipeline as input.

Before running CyLinter, transfer and organize mcmicro output data files (typically on the o2 cluster) according to the below tree diagram using the following command:

``` bash
prep <source_dir_path> <dest_dir_path>
```

```
<input_dir>
├── config.yml
├── csv
│   ├── unmicst-<sample1>.csv
│   └── unmicst-<sample2>.csv
├── markers.csv
└── tif
    ├── <sample1>.ome.tif
    └── <sample2>.ome.tif
```

In the case of TMA data, an additional "-t" flag must be passed prior to the source and destination file path variables

``` bash
prep -t <source_dir_path> <dest_dir_path>.
```

The resulting input directory will then be organized with TMA-specific subdirectories:

```
<input_dir>
├── <tma1>
│   ├── config.yml
│   ├── csv
│   │   ├── unmicst-<core1>.csv
│   │   └── unmicst-<core2>.csv
│   ├── markers.csv
│   └── tif
│       ├── <core1>.ome.tif
│       └── <core2>.ome.tif
└── <tma2>
    ├── config.yml
    ├── csv
    │   ├── unmicst-<core1>.csv
    │   └── unmicst-<core2>.csv
    ├── markers.csv
    └── tif
        ├── <core1>.ome.tif
        └── <core2>.ome.tif
```

Notes:

* The basename of the destination directory path (i.e. <dest_dir_path>) is assumed to be the experiment name.
* The file `markers.csv` contains metadata about experimental immunomarkers.
* Tabular data (i.e. `csv files`) are located in `csv/` subdirectories.
* Preprocessed multiplex imaging data (i.e. `ome.tif files`) are located in `tif/` subdirectories.
