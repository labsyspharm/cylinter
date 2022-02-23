---
layout: default-cylinter
title: Running CyLinter
nav_order: 4
has_children: false
---

# Running CyLinter

## Step 1:
Ensure that the desired configurations for a given analysis have been set in the [YAML configuration file]({{ site.baseurl }}/workflow/input#yaml-configuration-file).

## Step 2:
Activate the dedicated CyLinter virtual environment:  

``` bash
# on Mac:
source ~/cylinter/bin/activate

# on PC
cd C:\Users\<username>\cylinter\scripts
activate
```
## Step 3:
Execute the program from the beginning of the pipeline by passing the YAML configuration file to the `cylinter` command:  

``` bash
cylinter <input_dir>/config.yml
```
## Step 3.1:
CyLinter bookmarks progress by automatically caching partially-redacted feature tables as `.parquet` files in the `checkpoints/` directory of the CyLinter [output directory]({{ site.baseurl }}/workflow/output/index). To re-start the program at a desired module or re-run any of the [Modules]({{ site.baseurl }}/modules/index), pass the `--module` flag followed by a module name:

``` bash
cylinter --module <module-name> <input_dir>/config.yml
```
