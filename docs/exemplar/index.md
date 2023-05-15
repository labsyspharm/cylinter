---
layout: default-cylinter
title: Example data
nav_order: 7
---

# Example data

<img align="right" src="{{ site.baseurl }}/assets/images/cores.png" width="540" style="padding-left: 30px; padding-bottom: 20px;"> Four human tissue samples are provided as CyLinter demonstration data: normal kidney cortex, mesothelioma, glioblastoma, and normal tonsil. These imaging data were collected by [CyCIF](https://www.cycif.org/) and are derived from a tissue microarray collected at the [Laboratory of Systems Pharmacology](https://hits.harvard.edu/the-program/laboratory-of-systems-pharmacology/about/) referred to as Exemplar Microscopy Images of Tissues (EMIT TMA-22, Synapse ID: [syn22345750](https://www.synapse.org/#!Synapse:syn22345750)). The tissues were imaged at 20X magnification using a 0.75 NA objective and 2x2-pixel binning.

Access to the demonstration dataset requires free registration at the Sage Synapse data repository ([https://www.synapse.org/](https://www.synapse.org/)). Once registered, the example dataset can be downloaded with the following commands:

## Step 1: Download
``` bash
# Activate the CyLinter virtual environment.
mamba activate cylinter

# Install Synapse client.
mamba install -c bioconda synapseclient

# Mac users can download the demo dataset using the 'prep' command.
prep cylinter_demo ~/Desktop/cylinter_demo  # Enter Synapse ID and password when prompted.

# PC users, run the following command to download the demo dataset:  
synapse get -r syn25685780 --downloadLocation C:\Users\<username>\Desktop\cylinter_demo --multiThreaded
```
**PC users:** note that the downloaded CSV file names will be in all lowercase. These must be re-formatted to match the way they are specified in the template configuration file (`config.yml`) by capitalizing the 'm' in `cellmask` (e.g., `unmicst-1_cellMask`) for all four files.

## Step 2: Configure
After downloading the exemplar dataset, open the [YAML configuration file]({{ site.baseurl }}/workflow/input#yaml-configuration-file) in `~/Desktop/cylinter_demo/config.yml` and update the `inDir` and `outDir` parameters with user-specific directory paths. All other settings are pre-configured for use with the exemplar dataset.

```yaml
inDir: /Users/<username>/Desktop/cylinter_demo
outDir: /Users/<username>/Desktop/cylinter_demo/output
.
.
.
```

## Step 3: Run
Pass the configuration file to the `cylinter` command to run Cylinter on the example data:

``` bash
# for Mac:
cylinter --module <module-name>(optional) ~/Desktop/cylinter_demo/config.yml  

# for PC:
cylinter --module <module-name>(optional) C:\Users\<username>\Desktop\cylinter_demo\config.yml
```
