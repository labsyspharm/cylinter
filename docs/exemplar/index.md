---
layout: default-cylinter
title: Example data
nav_order: 7
---

# Example data

<img align="right" src="{{ site.baseurl }}/assets/images/cores.png" width="540" style="padding-left: 30px; padding-bottom: 20px;"> A set of four human tissue samples is provided as demonstration data: 2 lung squamous cell carcinomas, 1 mesothelioma, and 1 diverticulitis sample. These imaging data were collected using the [CyCIF method](https://www.cycif.org/) and are derived from a larger TMA dataset collected at the [Laboratory of Systems Pharmacology](https://hits.harvard.edu/the-program/laboratory-of-systems-pharmacology/about/) and referred to as Exemplar Microscopy Images of Tissues (EMIT) (Synapse ID: [syn22345750](https://www.synapse.org/#!Synapse:syn22345750)). The tissues were imaged at 20X magnification using a 0.75 NA objective and 2x2-pixel binning.

Access to the demonstration dataset requires free registration at the Sage Synapse data repository ([https://www.synapse.org/](https://www.synapse.org/)). Once registration is complete, the example dataset can be downloaded with the following commands:

``` bash
# Activate the CyLinter virtual environment.
source ~/cylinter/bin/activate

# Mac users can download the demo dataset using the 'prep' command.
prep cylinter_demo ~/Desktop/cylinter_demo  # Enter Synapse ID and password when prompted.

# PC users, run the following command to download the demo dataset:  
synapse get -r syn25685780 --downloadLocation C:\Users\<username>\Desktop\cylinter_demo --multiThreaded
```

Once downloaded, open the [YAML configuration file](input#yaml-configuration-file) in `~/Desktop/cylinter_demo/config.yml` and update the `in_dir` and `out_dir` parameters with user-specific directory paths; all other parameters are pre-configured for use with the exemplar dataset.

```yaml
in_dir: /Users/<user>/Desktop/cylinter_demo
out_dir: /Users/<user>/Desktop/cylinter_demo/output
.
.
.
```

Pass the configuration file to the `cylinter` command to run Cylinter on the example data:

``` bash
# for Mac:
cylinter --module <module-name> (optional) ~/Desktop/cylinter_demo/config.yml  

# for PC:
cylinter --module <module-name> (optional) C:\Users\<username>\Desktop\cylinter_demo\config.yml
```
