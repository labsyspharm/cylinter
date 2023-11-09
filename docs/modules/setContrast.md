---
layout: default-cylinter
title: setContrast
nav_order: 10
parent: Modules
---

10\. `setContrast`: in this module, image channel contrast is adjusted using the `contrast limits` slider bar in the `layer controls` dock at the top left of the Napari viewer. For each channel, contrast limits are set on a reference tissue whose median channel value is closest to the 85% quantile for the batch of tissues, such that a bright sample (but not the brightest, which tends to be driven by artifacts) is used to set the baseline contrast settings. The lower slider of the `contrast limits` slider bar is used to reduce background signal intensities, while the upper slider is used to increase channel gain. Once lower and upper sliders have been adjusted on the reference sample, the fit can be checked on other samples by entering their name in the `Sample Name` field the `Arbitrary Sample Selection` widget at the right of the Napari viewer and clicking the `Enter` button. Clicking the `Apply Limits and Move to Next Channel` button causes the module to move to the next channel for contrast adjustment. Contrast limits are applied to samples in the `gating` and `curteThumbnails` modules later in the pipeline. Channel contrast limits may be re-adjusted at any time by re-running the `setContrast` module.


### No YAML configurations
