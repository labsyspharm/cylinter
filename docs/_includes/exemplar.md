# Exemplar data

`exemplar-002` is standard output from the [mcmicro](https://github.com/labsyspharm/mcmicro) multiplex image preprocessing pipeline stored at Sage Synapse under synapseID: syn24193163. It consists of tabular and imaging data from a two-by-two cut-out of a TMA. The four cores are two meningioma tumors, one GI stroma tumor, and one normal colon specimen. The `markers.csv` file containing immunomarker metadata is also provided. Registration at [Sage Synapse](https://www.synapse.org/) is needed to download the example data.

After establishing a Synapse account, download `exemplar-002` and reorganize its data as CyLinter input with the following command:

``` bash
prep exemplar-002 <destination_path>
```

When prompted, enter your Synapse user ID and password, then press return.
`<destination_path>` is equivalent to the CyLinter input directory.

After your download is complete, open the YAML configuration file in the newly-reorganized CyLinter input directory with a text editor. Vim users may use the following command to open the file:

```bash
vim <destination_path>/exemplar-002/config.yml>
```

Replace the template contents with the following:
```yaml
in_dir: <input_path>/exemplar-002
out_dir: <output_path>/exemplar-002
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
source $HOME/cylinter/bin/activate  # activates the CyLinter virtual environment
cylinter <input_dir>/exemplar-002/config.yml
```
