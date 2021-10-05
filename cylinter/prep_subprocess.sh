#!/bin/bash

# Transfer mcmicro output files from cluster and format local directory structure for compatibility with CyLinter pipeline (run locally).

# USAGE: sh src/qc_prep.sh -t(required for TMA data) <source_dir> <destination_dir>

# EXAMPLE(TMA): sh /Users/greg/projects/cycif-qc/src/qc_prep.sh -t gjb15@transfer.rc.hms.harvard.edu:/n/scratch3/users/g/gjb15/exemplar-002 /Users/greg/projects/cycif-qc/input

shift $(( OPTIND - 1 ))

if [ "$2" == "cylinter_demo" ]; then

  # Transfer emit22_demo from Sage Synapse
  synapse get -r syn25685780 --downloadLocation "$3" --multiThreaded

  # Rename quantification and dearray subdirectories to "csv" and "tif", respectively.
  if [ -d "$3"/quantification ]; then
      mv "$3"/quantification "$3"/csv
      mv "$3"/dearray "$3"/tif
  fi

  # Copy config.yml template to CyLinter input directory
  cp "$4" "$3"

else

  if $1; then

    echo "Transferring TMA data."

    # Transfer mcmicro output files to CyLinter input directory.
    rsync -avP -m "$2"/ "$3" --include quantification/*.csv --include dearray/*.tif --include markers.csv --exclude work --exclude '*.*'

    mkdir -p "$3"/seg

    # Rename quantification and dearray subdirectories to "csv" and "tif", respectively.
    if [ -d "$3"/seg ]; then

      mv "$3"/*/quantification "$3"/csv
      mv "$3"/*/dearray "$3"/tif

      files=("$3"/csv/*)
      SEG=$(echo "$(basename "${files[0]}")" | cut -d'_' -f2 | cut -d'.' -f1)

      rsync -avP -m "$2"/ "$3" --include qc/s3seg/*/"$SEG"Outlines.ome.tif --exclude work --exclude '*.*'

      for RESOLVED_PATH in "$3"/*/qc/s3seg/* ; do
        SAMPLE_NAME=$(basename "$RESOLVED_PATH")
        arrIN=(${SAMPLE_NAME//-/ })
        NAME=${arrIN[1]}
        mv "$RESOLVED_PATH"/"$SEG"Outlines.ome.tif "$RESOLVED_PATH"/"$NAME".ome.tif
        mv "$RESOLVED_PATH"/"$NAME".ome.tif "$3"/seg/
      done

      for SAMPLE_PATH in "$3"/* ; do
        SAMPLE_NAME=$(basename "$SAMPLE_PATH")
        if [ $SAMPLE_NAME != "csv" ] && [ $SAMPLE_NAME != "tif" ] && [ $SAMPLE_NAME != "seg" ] && [ $SAMPLE_NAME != "markers.csv" ]; then
          mv "$SAMPLE_PATH"/markers.csv "$3"/
          rm -r "$SAMPLE_PATH"
        fi
      done
    fi

    # copy configuration template to input dir
    cp "$4" "$3"/

  else
    echo "Transferring whole tissue data."

    # Transfer mcmicro output files to CyLinter input directory.
    rsync -avP -m "$2"/ "$3" --include quantification/*.csv --include registration/*.ome.tif --include markers.csv --exclude work --exclude '*.*'

    # Make directories for images, data tables, and segmentation outlines
    mkdir -p "$3"/csv
    mkdir -p "$3"/tif
    mkdir -p "$3"/seg

    # combine sample tifs, csv files, and their segmentation outlines into respectively-labeled subdirectories.
    mv "$3"/*/quantification/*.csv "$3"/csv/
    mv "$3"/*/registration/*.tif "$3"/tif/

    SEG=$(echo "$(basename "$3"/csv/*)" | cut -d'_' -f2 | cut -d'.' -f1)
    rsync -avP -m "$2"/ "$3" --include qc/s3seg/*/"$SEG"Outlines.ome.tif --exclude work --exclude '*.*'

    for SAMPLE_PATH in "$3"/*/qc/s3seg/* ; do
      SAMPLE_NAME=$(basename "$SAMPLE_PATH")

      # crop off "unmicst-" prefix from outlines directory name
      arrIN=(${SAMPLE_NAME//unmicst-/ })

      mv "$SAMPLE_PATH/$SEG"Outlines.ome.tif "$SAMPLE_PATH"/"${arrIN[0]}".ome.tif
      mv "$SAMPLE_PATH"/"${arrIN[0]}".ome.tif "$3"/seg/
    done

    for SAMPLE_PATH in "$3"/* ; do
      SAMPLE_NAME=$(basename "$SAMPLE_PATH")
      if [ $SAMPLE_NAME != "csv" ] && [ $SAMPLE_NAME != "tif" ] && [ $SAMPLE_NAME != "seg" ] && [ $SAMPLE_NAME != "markers.csv" ]; then
        mv "$SAMPLE_PATH"/markers.csv "$3"/
        rm -r "$SAMPLE_PATH"
      fi
    done

    # Copy configuration template to input dir
    cp "$4" "$3"/

  fi
fi
