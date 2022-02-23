---
layout: default-cylinter
title: areaFilter
nav_order: 4
parent: Modules
---

4\. `areaFilter`: Cell segmentation errors introduce significant error into image-derived, single-cell data. In this module, users assign lower and upper cutoffs on cell segmentation area of cells comprising each tissue as calculated by the quantification module in [MCMICRO](https://mcmicro.org/) to remove under- and over-segmented cells. After gate placement, users can visualize selected cells in their corresponding tissue by clicking the `Plot Points` button beneath the histogram. Closing the Napari window causes the program to apply the current cutoffs and proceed to the next tissue for segmentation area cutoff assignment.

### No YAML configurations
