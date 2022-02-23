---
layout: default-cylinter
title: frequencyStats
nav_order: 13
parent: Modules
---

13\. `frequencyStats`: This module computes pairwise statistics for binary declarations specified in the [sampleMetadata]({{ site.baseurl }}/workflow/input#general-configurations) parameter of `config.yml`. Test results are saved in `<cylinter_output_path>/clustering/final/frequency_stats`. This module is automated and configurable (see YAML configures below).

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `controlGroups` | ["CANCER-FALSE"] | (list of strs) Corresponds to control groups for each binary declaration specified as the fourth elements of [sampleMetadata]({{ site.baseurl }}/workflow/input#general-configurations) values. |
|`denominatorCluster` | null | (null or int) Cluster to be used as the denominator when computing cluster frequency ratios. Set to null first, then change to cluster number (int) to normalize cluster frequencies to a particular identified cluster if desired. |
| `FDRCorrection` | False | (bool) Whether to compute p-vals and false discovery rate (FDR)-corrected q-vals (True) or compute uncorrected p-vals only (False). |
