---
layout: default
title: frequencyStats
nav_order: 14
parent: Modules
---

14\. `frequencyStats`: This module computes pairwise statistics for binary declarations specified in the `sampleMetadata` dictionary in [`config.yml`]({{ site.baseurl }}{% link input/index.md %}) (fourth dictionary values). Test results are saved in `<cylinter_output_path>/clustering/final/frequency_stats`. This module is automated and configurable; see `config.yml` for details.

### YAML configurations (`config.yml`)

| Parameter | Default | Description |
| --- | --- | --- |
| `controlGroups` | ["CANCER-FALSE"] (list or strs) | List of strings corresponding to the control groups for each of the binary declarations in the third elements of sample_metadata dict values |
| `denominatorCluster` | null (none type) | Cluster to be used as the denominator when computing cell type frequency ratios. Set to null first, change to cluster integer number to normalize cluster frequencies to a particular cluster. |
| `FDRCorrection` | False (bool) | If True, compute p-vals and false discovery rate (FDR) q-vals. If False, compute uncorrected p-vals only. |
