# Exemplar data

Multiplex immunofluorescence images of four different human tissue types (2 lung squamous cell carcinomas, 1 colon adenocarcinoma, and 1 diverticulitis) are used as a reference dataset for CyLinter testing and demonstration. The images were extracted from a reference dataset called Exemplar Microscopy Images of Tissues (EMIT TMA22; Synapse ID: [syn22345750](https://www.synapse.org/#!Synapse:syn22345750).

Registration at the Sage Synapse data repository ([https://www.synapse.org/](https://www.synapse.org/)) is free and required for downloading the CyLinter exemplar dataset.

Once registration is complete, download the exemplar data with the following commands:

``` bash
# Activate the CyLinter virtual environment.
source ~/cylinter/bin/activate

# Download the demo dataset using the prep command.
prep cylinter_demo ~/Desktop/cylinter_demo  # Enter Synapse ID and password when prompted.

# The prep command is currently only available for Mac. PC users, please use the following command to download the demo dataset:  
synapse get -r syn25685780 --downloadLocation C:\Users\<username>\Desktop\cylinter_demo --multiThreaded
```

Be sure to populate the `in_dir` and `out_dir` parameters in `~/Desktop/cylinter_demo/config.yml` with user-specific variables; all other metadata in `config.yml` is pre-configured for use with the exemplar dataset.

```yaml
in_dir: /Users/<user>/Desktop/cylinter_demo
out_dir: /Users/<user>/Desktop/cylinter_demo/output
.
.
.
```

Run CyLinter on the exemplar dataset:

``` bash
# Pass the YAML configuration file for the demonstration to CyLinter for analysis
cylinter --module (optional) ~/Desktop/cylinter_demo/config.yml  
```
