# Transferring and organizing input data
CyLinter can use standard output from the [mcmicro](https://github.com/labsyspharm/mcmicro) multiplex image processing pipeline as input and comes with a support console script called `prep` that programmatically organizes mcmicro output as CyLinter input.

 using the following command:

``` bash
# Transfer and organize mcmicro output into a CyLinter input directory.
prep <mcmicro_output_dir> <cylinter_input_dir>  # < > indicates a variable.
```

The `<mcmicro_output_dir>` is expected to contain sample-specific subdirectories and be organized in the following way:

``` bash
<mcmicro_output_dir>
├── <sample1/TMA1>
│   ├── markers.csv
│   ├── quantification
│   │   └── unmicst-<sample1/TMA1>.csv
│   └── registration
│       └── <sample1/TMA1>.ome.tif
└── <sample2/TMA2>
    ├── markers.csv
    ├── quantification
    │   └── unmicst-<sample2/TMA2>.csv
    └── registration
        └── <sample2/TMA2>.ome.tif
```

Because mcmicro stores TMA data differently than whole tissue data, pass a "-t" flag before the source and destination directory paths to prepare TMA data.

``` bash
prep -t <mcmicro_output_dir> <cylinter_input_dir>
```

SSH keys can used to transfer mcmicro output from remote sources (e.g. HMS o2 cluster).

``` bash
prep <userID>@transfer.rc.hms.harvard.edu:<mcmicro_output_dir> <cylinter_input_dir>
```

A correctly formatted CyLinter input directory will have the following form:

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

* Single-cell data files from a given TMA or batch of whole tissue sections will be combined into the `csv/` subdirectory.
* Corresponding multiplex imaging files (i.e. ome.tif files) will be combined in the `tif/` subdirectory.

# Immunomarker metadata
The `markers.csv` file contains metadata about immunomarkers used in the study and must have the following form:.

```
channel_number,cycle_number,marker_name,Filter,excitation_wavelength,emission_wavelength
1,1,<DNA1>,Hoecsht,395,431
2,1,<abx1>,FITC,485,525
3,1,<abx2>,Sytox,555,590
4,1,<abx3>,Cy5,640,690
5,2,<DNA2>,Hoecsht,395,431
6,2,<abx4>,FITC,485,525
7,2,<abx5>,Sytox,555,590
8,2,<abx6>,Cy5,640,690
.
.
.
```
