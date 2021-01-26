# Exemplar data

[mcmicro](https://github.com/labsyspharm/mcmicro) tabular and imaging output data from a two-by-two cut-out of a TMA (`exemplar-002`) are provided for demonstration purposes. A `markers.csv` file containing immunomarker metadata is also provided. The four cores are two meningioma tumors, one GI stroma tumor, and one normal colon specimen. The data are stored at Sage Synapse under synapseID: syn24193163. Registration at Sage Synapse is required to download the exemplar data; establish a free account at [Sage Synapse](https://www.synapse.org/).

``` bash
prep exemplar-002 <cylinter_input_dir>

# Enter your Synapse user ID and password when prompted, then press return.
```

Once download is complete, open and edit the YAML configuration file in the CyLinter input directory with a text editor. Vim users may do this from the terminal with the following command:

```bash
vim <cylinter_input_dir>/exemplar-002/config.yml>
```

Replace the configuration template contents with the following:

```yaml
in_dir: <cylinter_input_dir>/exemplar-002
out_dir: <cylinter_output_dir>/exemplar-002
random_sample_size: 1.0  # floating point (0.0-1.0); 1.0 is full dataset
mask_object: cellMask  # cellMask, nucleiMask, cellRingMask ...
sample_metadata:
  unmicst-1: ["tissue_a", "a", 1]  # <sampleString>: [<fullConditionString>, <abbrConditionString>, <replicateInteger>]
  unmicst-2: ["tissue_b", "b", 1]
  unmicst-3: ["tissue_c", "c", 1]
  unmicst-4: ["tissue_d", "d", 1]
samples_to_exclude: []  # [<sampleString1>, <sampleString2>, ...]
markers_to_exclude: ["AF488", "AF555", "AF647"]  # [<markerString1>, <markerString2>, ...]
```

Run CyLinter on `exemplar-002` with the following command:

``` bash
# Activate the CyLinter virtual environment
source $HOME/cylinter/bin/activate

# Pass exemplar-002 configuration file to CyLinter for analysis
cylinter --module (optional) <cylinter_input_dir>/exemplar-002/config.yml  
```
