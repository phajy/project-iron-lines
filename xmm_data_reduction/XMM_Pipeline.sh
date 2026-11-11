#!/bin/bash

function usage() {
    cat <<USAGE

    Usage: $0 [-o OBS_ID] [-i INSTRUMENT] [-s SIGMA]

    Options:
        -o, --obs-id:       Observation ID
        -i, --instrument:   Instrument
        -s:                 GTI standard deviation
        --source-max:       Maximum radius of source region
        --m1, --margin1:    Margin from primary flare of PSF
        --m2, --margin2:    Margin from edge of window

USAGE
    exit 1
}

# Checking for empty input
if [ $# -eq 0 ]; then
    usage
    exit 1
fi

# Finding the Directory of this script
DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default parameter arguments
OBS_ID=""
INSTRUMENT="pn"
SIGMA=1.5
SOURCE_MAX=1000
MARGIN1=5
MARGIN2=1

# Reading inputs
while [ "$1" != "" ]; do
    case $1 in
    -o | --obs-id)
        shift
        OBS_ID=$1
        ;;
    -i | --instrument)
        shift
        INSTRUMENT=$1
        ;;
    -s)
        shift
        SIGMA=$1
        ;;
    --source-max)
        shift
        SOURCE_MAX=$1
        ;;
    --m1 | --margin1)
        shift
        MARGIN1=$1
        ;;
    --m2 | --margin2)
        shift
        MARGIN2=$1
        ;;
    -h | --help)
        usage # run usage function
        ;;
    *)
        echo "Unknown option $1"
        usage
        exit 1
        ;;
    esac
    shift
done

# Checking for empty OBS_ID
if [ $OBS_ID == "" ]; then
    echo "You must pass an observation ID"
    usage
    exit 1
fi

export FOLDER=$(echo -n $(pwd) | tail -c 3)

if [ FOLDER != "ODF" ]; then
    echo "You must be in an ODF directory"
    usage
    exit 1
fi

# Running pipelines
bash $DIR/Pipelines.sh $SIGMA

cd ../PROC

# Running image and spectrum extraction
bash $DIR/Extraction.sh -o $OBS_ID -i $INSTRUMENT --m1 $MARGIN1 --m2 $MARGIN2 --source-max $SOURCE_MAX
