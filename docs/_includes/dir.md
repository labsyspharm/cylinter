# Input Directory Structure
CyLinter uses standard output files from the [mcmicro multiplex image preprocessing pipeline](https://github.com/labsyspharm/mcmicro) as input. Input directory structure is shown below:  

Before running CyLinter, transfer and organize the preprocessed data from o2 according to tree diagram (shown below) using the following command: `prep <source_dir_path> <dest_dir_path>`

```
<SAMPLE_dir>
├── config.yml
├── csv
│   ├── unmicst-<sample1>.csv
│   └── unmicst-<sample2>.csv
├── markers.csv
└── tif
    ├── <sample1>.ome.tif
    └── <sample2>.ome.tif
```

In the case of TMA data, an additional "-t" flag must be passed prior to the source and destination file path variables: `prep -t <source_dir_path> <dest_dir_path>`. The file will then be organized as shown below, with subdirectories containing files from different TMA experiments:

```
<TMA_dir>
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
* In the case of TMA data, data from different TMA slides are located in experiment-specific subdirectories.
