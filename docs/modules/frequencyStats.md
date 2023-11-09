---
layout: default-cylinter
title: frequencyStats
nav_order: 14
parent: Modules
---

14\. `frequencyStats`: this module is fully automated and computes pairwise statistics for binary declarations specified in the [sampleMetadata]({{ site.baseurl }}/workflow/input#general-configurations) parameter of `config.yml`. Test results are saved in `<output_dir>/clustering/2d/frequency_stats` in the case of 2D clustering or `<output_dir>/clustering/3d/frequency_stats` in the case of 3D clustering.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `controlGroups` | ["CANCER-FALSE"] | (list of strs) Corresponds to control groups for each binary declaration specified as the fourth elements of [sampleMetadata]({{ site.baseurl }}/workflow/input#general-configurations) values. |
|`denominatorCluster` | null | (null or int) Cluster to be used as the denominator when computing cluster frequency ratios. Set to null first, then change to cluster number (int) to normalize cluster frequencies to a particular identified cluster if desired. |
| `FDRCorrection` | False | (bool) Whether to compute p-values and false discovery rate (FDR)-corrected q-values (True) or compute uncorrected p-values only (False). |