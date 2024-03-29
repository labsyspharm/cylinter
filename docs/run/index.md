---
layout: default-cylinter
title: Running CyLinter
nav_order: 3
has_children: false
---

# Running CyLinter

## Step 1:
Ensure that the desired configurations for a given analysis have been set in CyLinter's [YAML configuration file]({{ site.baseurl }}/structure/index#yaml-configuration-file). A copy of this file can be found in the virtual environment into which CyLinter was installed (`.../miniconda3/envs/cylinter/lib/pythonXX/site-packages/cylinter/cylinter_config.yml`)

## Step 2:
Activate the CyLinter virtual environment:

``` bash
conda activate cylinter
```

## Step 3:
Execute the program from the beginning of the pipeline by passing the YAML configuration file (which should be stored at the top level of the CyLinter [input directory]({{ site.baseurl }}/structure/index)) to the `cylinter` command:  

``` bash
cylinter <input_dir>/cylinter_config.yml
```

CyLinter bookmarks progress by automatically caching partially-redacted spatial feature tables in the `checkpoints/` directory of the CyLinter [output directory]({{ site.baseurl }}/workflow/index). To re-run any of the [Modules]({{ site.baseurl }}/modules/index), pass the `--module` flag followed by the name of a specific module:

``` bash
cylinter --module <module-name> <input_dir>/cylinter_config.yml
```
