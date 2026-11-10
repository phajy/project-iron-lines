#!/bin/bash

# $1 is OBS_ID $2 is INSTRUMENT

# Finding the directory of this script
DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

# Selecting the correct filtering expression
if [ $2=="pn" ]; then
    export GTIExpression="${GTIFilters[0]}"
elif [ $2=="mos1" ]; then
    export GTIExpression="${GTIFilters[1]}"
else
    export GTIExpression="${GTIFilters[2]}"
fi
    
# Creating the GTI (Good Time Interval) file for time filtering
tabgtigen table="$2"_ltcrv.fits gtiset=gtiset.fits expression='$GTIExpression'

# Filtering the events file using the GTI
evselect table="$2"_filt.fits filtertype=expression filteredset="$2"_filt_time.fits expression='GTI(gtiset.fits,TIME)' 

# ==========================================================================================
# Image Generation
# ==========================================================================================

# Creating an attitude file
atthkgen atthkset=attitude.fits

if [ $2=="pn" ]; then
    export binsize=82
else
    export binsize=22
fi

echo "Bin size of "$binsize" selected for "$2""

# Creating an image for soft X-ray detections
evselect table="$2"_filt_time.fits withimageset=yes imageset="$2"-s.fits imagebinning=binSize xcolumn=X ximagebinsize=$binsize ycolumn=Y yimagebinsize=$binsize filtertype=expression expression='(FLAG == 0)&&(PI in [300:2000])'

# Creating an image for hard X-ray detections
evselect table="$2"_filt_time.fits withimageset=yes imageset="$2"-h.fits imagebinning=binSize xcolumn=X ximagebinsize=$binsize ycolumn=Y yimagebinsize=$binsize filtertype=expression expression='(FLAG == 0)&&(PI in [2000:10000])'

# Creating an image for all X-ray detections
evselect table="$2"_filt_time.fits withimageset=yes imageset="$2"-all.fits imagebinning=binSize xcolumn=X ximagebinsize=$binsize ycolumn=Y yimagebinsize=$binsize filtertype=expression expression='(FLAG == 0)&&(PI in [300:10000])'

# Creating an exposure map
eexpmap imageset="$2"-h.fits attitudeset=attitude.fits eventset="$2"_filt_time.fits expimageset="$2"_expmap_2-10keV.fits pimin=2000 pimax=10000

# ==========================================================================================
# Source Finding
# ==========================================================================================

python $DIR/FindSource.py $2

mapfile -t regions < $2-srcregion.txt

# Set region expressions as variables, replacing X, Y and R
export SOURCE_CENTRE="${regions[0]}"
export SOURCE_REGION=""$SOURCE_CENTRE",${regions[1]}"
export BACKGROUND_REGION_INNER=""$SOURCE_CENTRE",${regions[2]}"
export BACKGROUND_REGION_OUTER=""$SOURCE_CENTRE",${regions[3]}"

# ==========================================================================================
# Spectrum Extraction
# ==========================================================================================

if [ $2=="pn" ]; then
    # Extract the Source Spectrum
    evselect table="$2"_filt_time.fits energycolumn="PI" filteredset=""$2"_filtered.fits" filtertype="expression" expression="((X,Y) in CIRCLE("$SOURCE_REGION"))&&(FLAG==0)&&(PATTERN<=4)" spectrumset=""$2"_pi.fits" spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479

    # Extract the Background Spectrum
    evselect table="$2"_filt_time.fits energycolumn="PI" filteredset="bkg_filtered.fits" filtertype="expression" spectrumset="bkg_pi.fits" spectralbinsize=5 expression="((X,Y) in CIRCLE("$BACKGROUND_REGION_OUTER"))&&!((X,Y) in CIRCLE("$BACKGROUND_REGION_INNER"))&&(FLAG==0)&&(PATTERN<=4)" withspecranges=yes specchannelmin=0 specchannelmax=20479
else
    # Extract the Source Spectrum
    evselect table="$2"_filt_time.fits energycolumn="PI" filteredset=""$2"_filtered.fits" filtertype="expression" expression="((X,Y) in CIRCLE("$SOURCE_REGION"))" spectrumset=""$2"_pi.fits" spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999

    # Extract the Background Spectrum
    evselect table="$2"_filt_time.fits energycolumn="PI" filteredset="bkg_filtered.fits" filtertype="expression" spectrumset="bkg_pi.fits" spectralbinsize=5 expression="((X,Y) in CIRCLE("$BACKGROUND_REGION_OUTER"))&&!((X,Y) in CIRCLE("$BACKGROUND_REGION_INNER"))" withspecranges=yes specchannelmin=0 specchannelmax=11999
fi

# ==========================================================================================
# Pileup
# ==========================================================================================

# Check for Pileup
epatplot set="$2"_filtered.fits plotfile="$2"_epat.ps useplotfile=yes withbackgroundset=yes backgroundset=bkg_filtered.fits

# ==========================================================================================
# Response
# ==========================================================================================

# Find the source and background region areas for the ancillary file
backscale spectrumset="$2"_pi.fits badpixlocation="$2"_filt_time.fits
backscale spectrumset=bkg_pi.fits badpixlocation="$2"_filt_time.fits

# Generate the Photon Redistribution Matrix
rmfgen rmfset="$2"_rmf.fits spectrumset="$2"_pi.fits > rmflog.txt

# Generate the Ancillary File
arfgen arfset="$2"_arf.fits spectrumset="$2"_pi.fits withrmfset=yes rmfset="$2"_rmf.fits withbadpixcorr=yes badpixlocation="$2"_filt_time.fits > arflog.txt

# ==========================================================================================
# Cleanup
# ==========================================================================================

# Making a folder for the products
mkdir ../PRODUCTS

# Moving and renaming the product files
cp "$2"_pi.fits ../PRODUCTS/"$1"_"$2"_pi.fits
cp bkg_pi.fits ../PRODUCTS/"$1"_"$2"_bkg_pi.fits
cp "$2"_rmf.fits ../PRODUCTS/"$1"_"$2"_rmf.fits
cp "$2"_arf.fits ../PRODUCTS/"$1"_"$2"_arf.fits

# Enter the products directory
cd ../PRODUCTS

# Group the files
specgroup spectrumset="$1"_"$2"_pi.fits addfilenames=yes arfset="$1"_"$2"_arf.fits rmfset="$1"_"$2"_rmf.fits backgndset="$1"_"$2"_bkg_pi.fits
