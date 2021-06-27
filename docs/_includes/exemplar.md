# Exemplar data

Multiplex immunofluorescence images of four different human tissue types (2 lung squamous cell carcinomas, 1 colon adenocarcinoma, and 1 diverticulitis) are used as a reference dataset for CyLinter testing and demonstration. The images were extracted from a reference dataset called Exemplar Microscopy Images of Tissues (EMIT; Synapse ID: [syn22345750](https://www.synapse.org/#!Synapse:syn22345750).

Registration at the Sage Synapse data repository ([https://www.synapse.org/](https://www.synapse.org/)) is free and required for downloading the CyLinter exemplar dataset.

Once registration is complete, download the exemplar data with the following commands:

``` bash
# Activate the CyLinter virtual environment.
source ~/cylinter/bin/activate

prep emit22_demo ~/Desktop/emit22_demo  # Enter Synapse ID and password when prompted.
```

Be sure to populate the `in_dir` and `out_dir` parameters in `~/Desktop/emit22_demo/config.yml` with user-specific variables; all other metadata is in `config.yml` is pre-configured for use with the exemplar data.

```yaml
in_dir: /Users/<user>/Desktop/emit22_demo
out_dir: /Users/<user>/Desktop/emit22_demo/output
.
.
.
```

Run CyLinter on the exemplar data:

``` bash
# Pass the YAML configuration file for emit22_demo to CyLinter for analysis
cylinter --module (optional) ~/Desktop/emit22_demo/config.yml  
```
