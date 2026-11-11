#!/bin/bash

function usage() {
    cat <<USAGE

    Usage: $0 [-o OBS_ID] [-i INSTRUMENT] [--source-max SOURCE_MAX] [--m1 MARGIN1] [--m2 MARGIN2]

    Options:
        -o, --obs-id:       Observation ID
        -i, --instrument:   Instrument
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

# Finding the directory of this script
DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default parameter arguments
OBS_ID=""
INSTRUMENT="pn"
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

# Setting Environment Variables
cd ../ODF
export SAS_ODFPATH="$(pwd)"
export SAS_CCFPATH="/opt/local/XMM/ccf/"
export SAS_CCF="$SAS_ODFPATH/ccf.cif"
export SAS_ODF="$SAS_ODFPATH/$(ls *SUM.SAS)"
cd ../PROC

# ==========================================================================================
# GTI Filtering
# ==========================================================================================

mapfile -t GTIFilters < GTIFilters.txt

if [ $INSTRUMENT=="pn" ]; then
    export GTIExpression="${GTIFilters[0]}"
else
    echo "MOS REDUCTION NOT YET IMPLEMENTED"
    exit 1
fi

# # Selecting the correct filtering expression
# if [ $INSTRUMENT=="pn" ]; then
#     export GTIExpression="${GTIFilters[0]}"
# elif [ $INSTRUMENT=="mos1" ]; then
#     export GTIExpression="${GTIFilters[1]}"
# else
#     export GTIExpression="${GTIFilters[2]}"
# fi
    
# Creating the GTI (Good Time Interval) file for time filtering
tabgtigen table="$INSTRUMENT"_ltcrv.fits gtiset=gtiset.fits expression=$GTIExpression

# Filtering the events file using the GTI
evselect table="$INSTRUMENT"_filt.fits filtertype=expression filteredset="$INSTRUMENT"_filt_time.fits expression='GTI(gtiset.fits,TIME)' 

# ==========================================================================================
# Image Generation
# ==========================================================================================

# Creating an attitude file
atthkgen atthkset=attitude.fits

if [ $INSTRUMENT=="pn" ]; then
    export BINSIZE=82
else
    export BINSIZE=22
fi

echo "Bin size of "$BINSIZE" selected for "$INSTRUMENT""

# Creating an image for soft X-ray detections
evselect table="$INSTRUMENT"_filt_time.fits withimageset=yes imageset="$INSTRUMENT"-s.fits imagebinning=binSize xcolumn=X ximagebinsize=$BINSIZE ycolumn=Y yimagebinsize=$BINSIZE filtertype=expression expression='(FLAG == 0)&&(PI in [300:2000])'

# Creating an image for hard X-ray detections
evselect table="$INSTRUMENT"_filt_time.fits withimageset=yes imageset="$INSTRUMENT"-h.fits imagebinning=binSize xcolumn=X ximagebinsize=$BINSIZE ycolumn=Y yimagebinsize=$BINSIZE filtertype=expression expression='(FLAG == 0)&&(PI in [2000:10000])'

# Creating an image for all X-ray detections
evselect table="$INSTRUMENT"_filt_time.fits withimageset=yes imageset="$INSTRUMENT"-all.fits imagebinning=binSize xcolumn=X ximagebinsize=$BINSIZE ycolumn=Y yimagebinsize=$BINSIZE filtertype=expression expression='(FLAG == 0)&&(PI in [300:10000])'

# Creating an exposure map
eexpmap imageset="$INSTRUMENT"-h.fits attitudeset=attitude.fits eventset="$INSTRUMENT"_filt_time.fits expimageset="$INSTRUMENT"_expmap_2-10keV.fits pimin=2000 pimax=10000

# ==========================================================================================
# Source Finding
# ==========================================================================================

# Finding the source and background regions
python $DIR/FindSource.py $INSTRUMENT $MARGIN1 $MARGIN2 $SOURCE_MAX
mapfile -t regions < $INSTRUMENT-regions.txt

# Set region expressions as variables, replacing X, Y and R
export SOURCE_REGION="${regions[0]}"
export BACKGROUND_REGION_1="${regions[1]}"
export BACKGROUND_REGION_2="${regions[2]}"

# ==========================================================================================
# Spectrum Extraction
# ==========================================================================================

if [ $INSTRUMENT=="pn" ]; then
    # Extract the Source Spectrum
    evselect table=pn_filt_time.fits energycolumn="PI" filteredset="pn_filtered.fits" filtertype="expression" expression="((X,Y) in CIRCLE("$SOURCE_REGION"))&&(FLAG==0)&&(PATTERN<=4)" spectrumset="pn_pi.fits" spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479

    # Extract the Background Spectrum
    evselect table=pn_filt_time.fits energycolumn="PI" filteredset="bkg_filtered.fits" filtertype="expression" spectrumset="bkg_pi.fits" spectralbinsize=5 expression="((X,Y) in CIRCLE("$BACKGROUND_REGION_1"))||((X,Y) in CIRCLE("$BACKGROUND_REGION_2"))&&(FLAG==0)&&(PATTERN<=4)" withspecranges=yes specchannelmin=0 specchannelmax=20479
else
    # Extract the Source Spectrum
    evselect table="$INSTRUMENT"_filt_time.fits energycolumn="PI" filteredset=""$INSTRUMENT"_filtered.fits" filtertype="expression" expression="((X,Y) in CIRCLE("$SOURCE_REGION"))" spectrumset=""$INSTRUMENT"_pi.fits" spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999 > "$INSTRUMENT"_spec_log.txt

    # Extract the Background Spectrum
    evselect table="$INSTRUMENT"_filt_time.fits energycolumn="PI" filteredset="bkg_filtered.fits" filtertype="expression" spectrumset="bkg_pi.fits" spectralbinsize=5 expression="((X,Y) in CIRCLE("$BACKGROUND_REGION_1"))&&!((X,Y) in CIRCLE("$BACKGROUND_REGION_2"))" withspecranges=yes specchannelmin=0 specchannelmax=11999 >> "$INSTRUMENT"_spec_log.txt
fi

# ==========================================================================================
# Pileup
# ==========================================================================================

# Check for Pileup
epatplot set="$INSTRUMENT"_filtered.fits plotfile="$INSTRUMENT"_epat.ps useplotfile=yes withbackgroundset=yes backgroundset=bkg_filtered.fits > "$INSTRUMENT"_epatplot_log.txt

# ==========================================================================================
# Response
# ==========================================================================================

# Find the source and background region areas for the ancillary file
backscale spectrumset="$INSTRUMENT"_pi.fits badpixlocation="$INSTRUMENT"_filt_time.fits > "$INSTRUMENT"_backscale_log.txt
backscale spectrumset=bkg_pi.fits badpixlocation="$INSTRUMENT"_filt_time.fits >> "$INSTRUMENT"_backscale_log.txt

# Generate the Photon Redistribution Matrix
rmfgen rmfset="$INSTRUMENT"_rmf.fits spectrumset="$INSTRUMENT"_pi.fits > "$INSTRUMENT"_rmflog.txt

# Generate the Ancillary File
arfgen arfset="$INSTRUMENT"_arf.fits spectrumset="$INSTRUMENT"_pi.fits withrmfset=yes rmfset="$INSTRUMENT"_rmf.fits withbadpixcorr=yes badpixlocation="$INSTRUMENT"_filt_time.fits > "$INSTRUMENT"_arflog.txt

# ==========================================================================================
# Cleanup
# ==========================================================================================

# Making a folder for the products
mkdir ../PRODUCTS

# Moving and renaming the product files
cp "$INSTRUMENT"_pi.fits ../PRODUCTS/"$OBS_ID"_"$INSTRUMENT"_pi.fits
cp bkg_pi.fits ../PRODUCTS/"$OBS_ID"_"$INSTRUMENT"_bkg_pi.fits
cp "$INSTRUMENT"_rmf.fits ../PRODUCTS/"$OBS_ID"_"$INSTRUMENT"_rmf.fits
cp "$INSTRUMENT"_arf.fits ../PRODUCTS/"$OBS_ID"_"$INSTRUMENT"_arf.fits

# Enter the products directory
cd ../PRODUCTS

# Group the files
specgroup spectrumset="$OBS_ID"_"$INSTRUMENT"_pi.fits addfilenames=yes arfset="$OBS_ID"_"$INSTRUMENT"_arf.fits rmfset="$OBS_ID"_"$INSTRUMENT"_rmf.fits backgndset="$OBS_ID"_"$INSTRUMENT"_bkg_pi.fits
