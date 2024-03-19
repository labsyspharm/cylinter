---
layout: default-cylinter
title: frequencyStats
nav_order: 14
parent: Modules
---

14\. `frequencyStats`: this module is fully automated. It computes pairwise statistics for binary declarations specified in the [sampleMetadata]({{ site.baseurl}}/structure/#general-configurations) parameter of `cylinter_config.yml`. Test results are saved to a directory called `frequency_stats` in the clustering subdirectory of the main CyLinter output directory. This path is `clustering/2d/frequency_stats` in the case of 2D clusterings and `clustering/3d/frequency_stats` in the case of 3D clusterings.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `controlGroups` | ["CANCER-FALSE"] | (list of strs) Corresponds to control groups for each binary declaration specified as the fourth elements of [sampleMetadata]({{ site.baseurl }}/workflow/input#general-configurations) values. |
|`denominatorCluster` | null | (null or int) Cluster to be used as the denominator when computing cluster frequency ratios. Set to null first, then change to cluster number (int) to normalize cluster frequencies to a particular identified cluster if desired. |
| `FDRCorrection` | False | (bool) Whether to compute p-values and false discovery rate (FDR)-corrected q-values (True) or compute uncorrected p-values only (False). |