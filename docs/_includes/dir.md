# Transferring and organizing CyLinter input data
CyLinter takes standard output from the [mcmicro](https://github.com/labsyspharm/mcmicro) multiplex image preprocessing pipeline as input and comes with a support console script called `prep` that programmatically organizes this output as CyLinter input.

``` bash
# Transfer mcmicro output to CyLinter input directory and organize files for analysis
prep <mcmicro_output_dir> <cylinter_input_dir>

# mcmicro stores TMA data differently than whole tissue data. A "-t" flag must be passed before the source and destination paths to indicate that the mcmicro output are TMA data.
prep -t <mcmicro_output_dir> <cylinter_input_dir>
```

The <mcmicro_output_dir> should contain subdirectories containing all samples (or TMAs) from a given experiment:

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

* Multiplex imaging files (i.e. ome.tif files) for a given TMA or set of whole tissue sections are combined in the `tif/` subdirectory.
* Corresponding tabular data are combined in the `csv/` subdirectory.
* The `markers.csv` file contains metadata about immunomarkers used in the experiment and should take the following form:.

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
