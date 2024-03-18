---
layout: default-cylinter
title: setContrast
nav_order: 10
parent: Modules
---

10\. `setContrast`: in this module, image channel contrast is adjusted using the `contrast limits` slider bar in the `layer controls` dock at the top left of the Napari viewer. For each channel, contrast limits are set on a reference image whose median channel value is nearest to the 85th quantile of tissues in the batch which are applied to that image channel for all tissues in a batch. The 85th quantile (not 100th) is chosen to avoid picking tissue whose channel intensity is drive by bright artifacts outliers sample. The lower slider of the `contrast limits` slider bar is used to reduce background signal intensities by sliding to the right, while the upper slider is used to increase channel gain by sliding to the left. Once lower and upper sliders have been adjusted on the reference sample, the fit can be checked against other tissues in the batch by entering their name in the `Sample Name` field the `Arbitrary Sample Selection` widget at the right of the Napari viewer and clicking the `Enter` button. Clicking the `Apply Limits and Move to Next Channel` button causes the module to move to the next channel for contrast adjustment. To re-define contrast settings, simply re-run the `setContrast` module with `cylinter --module setContrast cylinter_config.yml`.

<!-- Once contrast limits have been defined they will automatically be applied to any module in which image channels are shown (e.g., `selectROIs`, `gating` and `curteThumbnails`, etc.)  -->


### No YAML configurations
