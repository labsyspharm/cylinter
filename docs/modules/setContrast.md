---
layout: default
title: setContrast
nav_order: 12
parent: Modules
---

12\. `setContrast`: The Napari image viewer is used to adjust image contrast settings on a reference tissue specified by the `viewSample` parameter in `config.yml`. Contrast settings are applied to all other tissues in the analysis before images are curated by the `curatethumbnails` module. Users toggle channels on and off by clicking their respectively labeled buttons at the left of the Napari window and moving the lower and upper sliders on the "contrast limits" control bar in the upper left of the Napari window. The lower slider thresholds background signal intensities, while the upper slider increases channel gain. Closing the Napari viewer saves the current settings as key:value pairs in `<output_dir>/contrast/contrast_limits.yml` and causes the program to proceed to the next module. Remove the `contrast_limits.yml` file (or edit specific key:value pairs) to re-define image contrast settings.

### YAML configurations

| Parameter | Default | Description |
| --- | --- | --- |
| `viewSample` | 39 (str) | Sample used to assign image contrast settings that will be applied to all samples in the analysis. Corresponds to first elements of `sampleMetadata` dictionary values |
