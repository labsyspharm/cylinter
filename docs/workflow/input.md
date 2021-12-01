---
layout: default-cylinter
title: Input
nav_order: 1
parent: Workflow
---

{: .no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

# Transferring and organizing input
CyLinter supports analysis of **whole slide image (WSI)** and **tissue microarray (TMA)** data processed by the [MCMICRO](https://mcmicro.org) image-processing pipeline. The program comes with a secondary command called `prep` which programmatically transfers and organizes standard MCMICRO output as CyLinter input (currently only available for Mac users). However, CyLinter can analyze any multiplex imaging data so long as they are in the correct format: TIFF/OME-TIFF images, OME-TIFF cell segmentation outlines, and CSV single-cell feature tables.

<br/>

Mac users, run the following command to organize **WSI DATA** as CyLinter input:  


``` bash
prep <mcmicro_output_dir> <cylinter_input_dir>  # < > indicates a variable.
```

Use the **"-t"** flag to organize **TMA DATA** as CyLinter input:

``` bash
prep -t <mcmicro_output_dir> <cylinter_input_dir>
```

 * *PC and linux users should organize CyLinter input data according to the below [input directory structure](#input-directory-structure).*

<br/>

For those with access to the o2 compute cluster at HMS, SSH keys may be used to transfer MCMICRO output from o2 compute to a local workstation using the following command:

``` bash
prep <userID>@transfer.rc.hms.harvard.edu:<mcmicro_output_dir> <cylinter_input_dir>
```

<br/>

# Input directory structure

CyLinter assumes input files take the standard MCMICRO output format. In the below example, `<sample-name>` corresponds to the names of the different samples used in a particular study.

``` bash
<input_dir>
├── config.yml
├── csv/
│   └── unmicst-<sample-name>.csv
├── markers.csv
├── seg/
│   └── <sample-name>.ome.tif (or tif)
└── tif/
    └── <sample-name>.ome.tif (or tif)
```

<br/>

# YAML configuration file

`config.yml` is the YAML configuration file passed to the `cylinter` command when the program is executed. It specifies module-specific parameters and general program configuration settings for a given analysis and should be included in the top-level CyLinter [input directory](#input-directory-structure). The template is pre-configured for use with [Example Data]({{ site.baseurl }}/exemplar) used to demonstrate CyLinter.

## General configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `inDir` | /Users/user/Desktop/cylinter_demo | Path to CyLinter input directory containing multi-channel image files (TIFF or OME-TIFF), segmentation outlines (OME-TIFF), and corresponding single-cell feature tables (CSV) |
| `outDir` | /Users/user/Desktop/cylinter_demo/output | Path to arbitrary CyLinter output directory (created automatically) |
| `randomSampleSize` | 1.0 (float) | Analyze a fraction of single-cell data selected at random (range: 0.0-1.0) |
| `sampleMetadata` | unmicst\<sample-name>\_cellMask: <br />  ["13", "Lung squamous cell carcinoma", "LSC", "CANCER-TRUE", 1] (dict) | Sample metadata dictionary: keys are file names; values are ordered lists of strings. First element: sample names (str), second element: descriptive text of experimental condition (str), third element: abbreviated version of experimental condition (str), fourth element: comma-delimited string of binary declarations across samples for computing t statistics between two groups (str dytpe), fifth element: integer replicate number for biological/technical replicates |
| `samplesToExclude` | [ ] (list of strs) | List of sample names (strs) to exclude from the analysis: first elements of sampleMetadata dict values |
| `markersToExclude` | [ ] (list of strs) | List of immunomarkers (strs) to exclude from the analysis (this does not include nuclear dye channels) |

## Module configurations
For module-specific configurations, see [Modules]({{ site.baseurl }}/modules)


<br/>

# Markers.csv
`markers.csv` is a standard input file into the MCMICRO image-processing pipeline that is also used by CyLinter as immunomarker metadata. The file takes the following format and should be included in the top-level CyLinter [input directory](#input-directory-structure):

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
* Additional metadata columns can be present in the file, but are not read by CyLinter.
