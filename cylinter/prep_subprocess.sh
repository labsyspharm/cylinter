#!/bin/bash

# Transfer mcmicro output files from cluster and format local directory structure for compatibility with CyLinter pipeline (run locally).

# USAGE: sh src/qc_prep.sh -t(required for TMA data) <source_dir> <destination_dir>

# EXAMPLE(TMA): sh /Users/greg/projects/cycif-qc/src/qc_prep.sh -t gjb15@transfer.rc.hms.harvard.edu:/n/scratch3/users/g/gjb15/exemplar-002 /Users/greg/projects/cycif-qc/input

shift $(( OPTIND - 1 ))

if [ "$2" == "exemplar-002" ]; then
  if [ ! -d "$3"/"$2" ]; then
    synapse get -r syn24193163 --downloadLocation "$3" --multiThreaded  # download exemplar-002 from Sage Synapse

    # Rename quantification and dearray subdirectories to "csv" and "tif", respectively.
    ROOT_NAME=$(basename "$2")
    if [ -d "$3"/"$ROOT_NAME"/quantification ]; then
        mv "$3"/"$ROOT_NAME"/quantification "$3"/"$ROOT_NAME"/csv
        mv "$3"/"$ROOT_NAME"/dearray "$3"/"$ROOT_NAME"/tif
    fi

    # copy configuration template to input dir
    cp "$4" "$3"/"$ROOT_NAME"/
  fi
else
  if $1; then
    echo "Transferring TMA data."
    rsync -avP -m "$2" "$3" --include quantification/*.csv --include dearray/*.tif --include markers.csv --exclude work --exclude '*.*'

    # Rename quantification and dearray subdirectories to "csv" and "tif", respectively.
    ROOT_NAME=$(basename "$2")
    if [ -d "$3"/"$ROOT_NAME"/quantification ]; then
        mv "$3"/"$ROOT_NAME"/quantification "$3"/"$ROOT_NAME"/csv
        mv "$3"/"$ROOT_NAME"/dearray "$3"/"$ROOT_NAME"/tif
    fi

    # copy configuration template to input dir
    cp "$4" "$3"/"$ROOT_NAME"/

  else
    echo "Transferring whole tissue data."
    rsync -avP -m "$2" "$3" --include quantification/*.csv --include registration/*.tif --include markers.csv --exclude work --exclude '*.*'

    # Make directories for tabular and imaging data.
    ROOT_NAME=$(basename "$2")
    mkdir -p "$3"/"$ROOT_NAME"/csv
    mkdir -p "$3"/"$ROOT_NAME"/tif

    # Combine sample csv tables and tifs into respective "csv" and "tif" subdirectories.
    for SAMPLE_PATH in "$3"/"$ROOT_NAME"/* ; do
      SAMPLE_NAME=$(basename $SAMPLE_PATH)
      if [ $SAMPLE_NAME != "csv" ] && [ $SAMPLE_NAME != "tif" ]; then
          echo $SAMPLE_NAME
          if [ -d $SAMPLE_PATH/quantification ]; then
              mv $SAMPLE_PATH/quantification/*.csv "$3"/"$ROOT_NAME"/csv/
              mv $SAMPLE_PATH/registration/*.tif "$3"/"$ROOT_NAME"/tif/
              mv $SAMPLE_PATH/markers.csv "$3"/"$ROOT_NAME"/  # overwrites markers files with every sample in loop
          fi
          rm -r $SAMPLE_PATH
      fi
    done

    # copy configuration template to input dir
    cp "$4" "$3"/"$ROOT_NAME"/

  fi
fi
