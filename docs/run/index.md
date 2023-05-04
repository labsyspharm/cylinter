---
layout: default-cylinter
title: Running CyLinter
nav_order: 4
has_children: false
---

# Running CyLinter

## Step 1:
Ensure that the desired configurations for a given analysis have been set in CyLinter's [YAML configuration file]({{ site.baseurl }}/workflow/input#yaml-configuration-file). A copy of this file can be found in the Mamba environment into which CyLinter was installed (`.../mambaforge/envs/cylinter/config.yml`)

## Step 2:
Activate the Mamba environment:

``` bash
mamba activate cylinter
```

## Step 3:
Execute the program from the beginning of the pipeline by passing the YAML configuration file to the `cylinter` command:  

``` bash
cylinter <input_dir>/config.yml
```

CyLinter bookmarks progress by automatically caching partially-redacted feature tables in the `checkpoints/` directory of the CyLinter [output directory]({{ site.baseurl }}/workflow/output/index). To re-start the program at a desired module or re-run any of the [Modules]({{ site.baseurl }}/modules/index), pass the `--module` flag followed by a module name:

``` bash
cylinter --module <module-name> <input_dir>/config.yml
```
