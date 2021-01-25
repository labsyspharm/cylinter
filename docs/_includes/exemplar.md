# Exemplar data

An example is provided for demonstration purposes:

* `exemplar-002` is standard output from the [mcmicro](https://github.com/labsyspharm/mcmicro) multiplex image preprocessing pipeline stored at Sage Synapse under synapseID: syn22241433. It consists of tabular and imaging data from a two-by-two cut-out of a TMA. The four cores are two meningioma tumors, one GI stroma tumor, and one normal colon specimen. A file containing immunomarker metadata is also provided.

* Registration at [Sage Synapse](https://www.synapse.org/) is needed to download the example data.

After registration, download `exemplar-002` and format its contents as CyLinter input.

``` bash
prep exemplar-002 </input/dir>

```
You will then be asked to enter your Synapse userID and password before your download begins.
`</input/dir>` points to the directory where the formatted CyLinter input should be stored.

After download is complete, open the YAML configuration file with a text editor. Users experienced with Vim may use the following command to open the file:

```bash
vim </input/dir>/exemplar-002/config.yml>
```

Replace the template contents with the following:
```yaml
in_dir: </input/dir>/exemplar-002
out_dir: <output/dir>/exemplar-002
random_sample_size: 1.0  # floating point percentage value between 0.0-1.0
mask_object: cellMask  # cellMask, nucleiMask, cellRingMask ...
sample_metadata:
  unmicst-1: ["tissue_a", "a", 1] # <sampleString>: [<fullConditionString>, <abbrConditionString>, <replicateInteger>]
  unmicst-2: ["tissue_b", "b", 1]
  unmicst-3: ["tissue_c", "c", 1]
  unmicst-4: ["tissue_d", "d", 1]
samples_to_exclude: []  # [<sampleString1>, <sampleString2>, ...]
markers_to_exclude: []  # [<markerString1>, <markerString2>, ...]
```

Run CyLinter on `exemplar-002` with following command:

``` bash
cylinter </input/dir>/exemplar-002/config.yml
```
