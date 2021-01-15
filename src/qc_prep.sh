#!/bin/bash

# Transfer mcmicro output files from cluster and format local directory structure for compatibility with qc pipeline (run locally).

# USAGE: sh src/qc_prep.sh -t(required for TMA data) <source_dir> <destination_dir>

# EXAMPLE(TMA): sh /Users/greg/projects/cycif-qc/src/qc_prep.sh -t gjb15@transfer.rc.hms.harvard.edu:/n/scratch3/users/g/gjb15/exemplar-002 /Users/greg/projects/cycif-qc/input

# ARCHIVED MCMICRO OUTPUTS:
# topacio: /n/files/ImStor/sorger/data/RareCyte/Tuulia/mcmicro/Topacio/21092020_s/pro
# emit tma22: /n/files/ImStor/sorger/analysis/EMIT/mcmicro-reference/TMA22
# sardana: /n/files/sorger/data/RareCyte
# gl261_vs_ct2a: /n/scratch3/users/g/gjb15/gl261_vs_ct2a
# exemplar-002: /n/scratch3/users/g/gjb15/exemplar-002

# Parse TMA flag.
has_t_option=false
while getopts :ht opt; do
    case $opt in
        h) show_some_help; exit ;;
        t) has_t_option=true ;;
        :) echo "Missing argument for option -$OPTARG"; exit 1;;
       \?) echo "Unknown option -$OPTARG"; exit 1;;
    esac
done
shift $(( OPTIND - 1 ))

if $has_t_option; then
  echo "Transferring TMA data."
  rsync -avP -m "$1" "$2" --include quantification/*.csv --include dearray/*.tif --include markers.csv --exclude work --exclude '*.*'

  # Make directories for tabular and imaging data.
  ROOT_NAME=$(basename "$1")

  # Combine sample csv tables and tifs into respective "csv" and "tif" subdirectories.
  for TMA_PATH in "$2"/"$ROOT_NAME"/* ; do
    TMA_NAME=$(basename $TMA_PATH)
    echo $TMA_NAME
    if [ -d $TMA_PATH/quantification ]; then
        mv $TMA_PATH/quantification $TMA_PATH/csv
        mv $TMA_PATH/dearray $TMA_PATH/tif
    fi

    cp /Users/greg/projects/cycif-qc/src/qc_config.yml "$2"/"$ROOT_NAME"/"$TMA_NAME"

  done

else
  echo "Transferring tissue section data."
  rsync -avP -m "$1" "$2" --include quantification/*.csv --include registration/*.tif --include markers.csv --exclude work --exclude '*.*'

  # Make directories for tabular and imaging data.
  ROOT_NAME=$(basename "$1")
  mkdir -p "$2"/"$ROOT_NAME"/csv
  mkdir -p "$2"/"$ROOT_NAME"/tif

  # Combine sample csv tables and tifs into respective "csv" and "tif" subdirectories.
  for SAMPLE_PATH in "$2"/"$ROOT_NAME"/* ; do
    SAMPLE_NAME=$(basename $SAMPLE_PATH)
    if [ $SAMPLE_NAME != "csv" ] && [ $SAMPLE_NAME != "tif" ]; then
        echo $SAMPLE_NAME
        if [ -d $SAMPLE_PATH/quantification ]; then
            mv $SAMPLE_PATH/quantification/*.csv "$2"/"$ROOT_NAME"/csv/
            mv $SAMPLE_PATH/registration/*.tif "$2"/"$ROOT_NAME"/tif/
            mv $SAMPLE_PATH/markers.csv "$2"/"$ROOT_NAME"/  # overwrites markers files with every sample in loop
        fi
        rm -r $SAMPLE_PATH
    fi
  done

  cp /Users/greg/projects/cycif-qc/src/qc_config.yml "$2"/"$ROOT_NAME"/

fi
