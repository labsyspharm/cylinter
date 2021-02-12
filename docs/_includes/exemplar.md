# Exemplar data

[mcmicro](https://github.com/labsyspharm/mcmicro) output data from a two-by-two cut-out of a TMA (`exemplar-002`) are provided for demonstration purposes along with a `markers.csv` file containing immunomarker metadata. The four cores are two meningioma tumors, one GI stroma tumor, and one normal colon specimen. The data are stored at Sage Synapse under synapseID: syn24193163.

Registration at Sage Synapse is required to download the exemplar data; registration is free at [https://www.synapse.org/](https://www.synapse.org/).

``` bash
# Activate the CyLinter virtual environment (if not already).
source ~/cylinter/bin/activate

prep exemplar-002 ~/Desktop  # Enter Synapse ID and password when prompted.
```

Once the download is complete, edit the `in_dir` and `out_dir` parameters of the template YAML configuration file (`~/Desktop/exemplar-002/config.yml`) previously added to the CyLinter input directory by the `prep` command. Note: all other metadata parameters have been pre-configured for `exemplar-002`.

```yaml
in_dir: /Users/<user>/Desktop/exemplar-002
out_dir: /Users/<user>/Desktop/output
.
.
.
```

Run CyLinter on `exemplar-002`:

``` bash
# Pass the YAML configuration file for exemplar-002 to CyLinter for analysis
cylinter --module (optional) ~/Desktop/exemplar-002/config.yml  
```
