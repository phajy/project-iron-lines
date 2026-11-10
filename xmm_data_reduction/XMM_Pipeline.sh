#!/bin/bash

# $1 is OBS_ID $2 is INSTRUMENT

# Finding the Directory of this script
DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

bash $DIR/Pipelines.sh

cd ../PROC

bash $DIR/Extraction.sh $1 $2
