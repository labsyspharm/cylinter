# Input Directory Structure
CyLinter can take standard output files from the [mcmicro multiplex image preprocessing pipeline](https://github.com/labsyspharm/mcmicro) as input into the pipeline.  

Run `prep.sh` to configure mcmicro output files into CyLinter-compatible input directory structure:

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

* The name of the parent input directory (e.g., `<input>`) is the experiment name.
* The file `markers.csv` contains metadata about experimental immunomarkers and must be placed in the parent directory.
* Tabular data (i.e. `csv files`) must be placed inside the `csv/` subdirectory.
* Preprocessed multiplex imaging data (i.e. `ome.tif files`) must be placed inside the `tif/` subdirectory.
* For TMA data (as opposed to whole tissue samples), subdirectories labeled by TMA name will be present in the parent input directory (`<input>`). Otherwise, the input directory structure is repeated per TMA subdirectory.
* The pipeline can be started from any module (see below), in which case it will look for the cached output of the immediately preceding module in the auto-generated `output/checkpoints` subdirectory.
  * For example, running the pipeline with `--module selectROIs` will start the pipeline at the region of interest (ROI) section module.
