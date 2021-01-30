# Exemplar data

[mcmicro](https://github.com/labsyspharm/mcmicro) tabular and imaging output data from a two-by-two cut-out of a TMA (`exemplar-002`) are provided for demonstration purposes along with a `markers.csv` file containing immunomarker metadata. The four cores are two meningioma tumors, one GI stroma tumor, and one normal colon specimen. The data are stored at Sage Synapse under synapseID: syn24193163. Registration at Sage Synapse is required to download the exemplar data; you can establish a free account at [https://www.synapse.org/](https://www.synapse.org/).

``` bash
# Activate the CyLinter virtual environment (if not already)
source $HOME/cylinter/bin/activate

prep exemplar-002 <cylinter_input_dir>  # Enter Synapse ID and password when prompted
```

Once the download is complete, edit the template YAML configuration file (`config.yml`) added to the CyLinter input directory by the `prep` command. Vim users may use the following command:

```bash
vim <cylinter_input_dir>/exemplar-002/config.yml>
```

Replace the contents of the template configuration file with the following:

```yaml
in_dir: "<cylinter_input_dir>/exemplar-002"
out_dir: "<cylinter_output_dir>/exemplar-002"
random_sample_size: 1.0  # floating point (0.0-1.0); 1.0 is full dataset
mask_object: cellMask  # cellMask, nucleiMask, cellRingMask ...
sample_metadata:
  unmicst-1: ["tissue_a", "a", 1]  # <sampleString>: [<fullConditionString>, <abbrConditionString>, <replicateInteger>]
  unmicst-2: ["tissue_b", "b", 1]
  unmicst-3: ["tissue_c", "c", 1]
  unmicst-4: ["tissue_d", "d", 1]
samples_to_exclude: []  # [<sampleString1>, <sampleString2>, ...]
markers_to_exclude: ["AF488", "AF555", "AF647", A488_background", "A555_background", "A647_background"]  # [<markerString1>, <markerString2>, ...]
```

Make sure the `in_dir` and `out_dir` fields are updated with the correct `<cylinter_input_dir>` and `cylinter_output_dir` directory paths.

Run CyLinter on `exemplar-002`:

``` bash
# Pass the YAML configuration file for exemplar-002 to CyLinter for analysis
cylinter --module (optional) <cylinter_input_dir>/exemplar-002/config.yml  
```
