# Exemplar data

Multiplex immunofluorescence data of four human tissue cores from a reference TMA dataset: Exemplar Microscopy Images of Tissues (EMIT), are provided for the purpose of demonstrating CyLinter usage (Synapse ID: [syn22345750](https://www.synapse.org/#!Synapse:syn22345750). The data consist of two lung squamous cell carcinoma cores, one colon adenocarcinoma core, and one diverticulitis core.

Registration at the Sage Synapse data repository ([https://www.synapse.org/](https://www.synapse.org/)) is free and required before downloading the exemplar data using the following command:

``` bash
# Activate the CyLinter virtual environment.
source ~/cylinter/bin/activate

prep emit22_demo ~/Desktop/emit22_demo  # Enter Synapse ID and password when prompted.
```

Be sure to populate the `in_dir` and `out_dir` parameters in the YAML configuration file (`~/Desktop/emit22_demo/config.yml`) with user-specific variables before running the program; all other metadata is pre-configured for use with the EMIT-22 demonstration data.

```yaml
in_dir: /Users/<user>/Desktop/emit22_demo
out_dir: /Users/<user>/Desktop/emit22_demo/output
.
.
.
```

Run CyLinter on `emit22_demo`:

``` bash
# Pass the YAML configuration file for emit22_demo to CyLinter for analysis
cylinter --module (optional) ~/Desktop/emit22_demo/config.yml  
```
