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
CyLinter supports analysis of **whole slide image (WSI)** and **tissue microarray (TMA)** data processed by the [MCMICRO](https://mcmicro.org) image-processing pipeline. The program comes with a helper command called `prep` that programmatically transfers and organizes standard MCMICRO output as CyLinter input (currently only available for Mac). However, CyLinter can analyze any multiplex imaging data so long as they conform to the correct file formats: TIFF/OME-TIFF images, TIFF/OME-TIFF cell segmentation outlines, and CSV single-cell feature tables with a cell segmentation area column labeled "Area".

<br/>

Mac users, run the following command to organize **WSI DATA** as CyLinter input:  


``` bash
prep <mcmicro_output_dir> <cylinter_input_dir>  # < > indicates a variable.
```

Use the **"-t"** flag to organize **TMA DATA** as CyLinter input:

``` bash
prep -t <mcmicro_output_dir> <cylinter_input_dir>
```

 * PC and linux users should organize CyLinter input data according to the standard [input directory structure](#input-directory-structure).

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

`config.yml` is the YAML configuration file passed to the `cylinter` command when the program is executed. It specifies general program configuration settings and module-specific parameters for a given analysis and should be included in the top-level CyLinter [input directory](#input-directory-structure). The `config.yml` file downloaded with the program is pre-configured for use with [Example Data]({{ site.baseurl }}/exemplar) used to demonstrate CyLinter.

## General configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `inDir` | /Users/user/Desktop/cylinter_demo | CyLinter input directory. Contains multi-channel image files (OME-TIFF or TIFF), segmentation outlines (OME-TIFF or TIFF), corresponding single-cell feature tables (CSV), `config.yml`, and `markers.csv` organized according to [input directory structure.](#input-directory-structure) |
| `outDir` | /Users/user/Desktop/cylinter_demo/output |CyLinter output directory. Path is created if it does not already exist. |
| `sampleMetadata` | unmicst-\<sample-name>\_cellMask: <br />  ["15", "Glioblastoma", "GBM", "CANCER-TRUE", 1] | Sample metadata dictionary: keys = file names; values = list of strings. First elements: sample names (str). Second elements: descriptive text of experimental condition (str). Third elements: abbreviation of experimental condition (str). Fourth elements: comma-delimited string of arbitrary binary declarations for computing t-statistics between two groups of samples (str dytpe). Fifth elements: replicate number specifying biological or technical replicates (int). |
| `samplesToExclude` | [ ] | (list of strs) Sample names to exclude from analysis specified according to the first elements of sampleMetadata configuration. |
| `markersToExclude` | [ ] | (list of strs) Immunomarkers to exclude from analysis (does not include nuclear dyes). |

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
