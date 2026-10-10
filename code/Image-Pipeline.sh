#!/bin/bash

# $1 = Instrument

# Filtering the events file using the GTI
evselect table="$1"_filt.fits filtertype=expression filteredset="$1"_filt_time.fits expression='GTI(gtiset.fits,TIME)' 

# Creating an attitude file
atthkgen atthkset=attitude.fits
W
if [ $1=="pn" ]; then
    export binsize=82
    echo "Bin size of "$binsize" selected for "$1""
else
    export binsize=22
    echo "Bin size of "$binsize" selected for "$1""
fi

# Creating an image for soft X-ray detections
evselect table="$1"_filt_time.fits withimageset=yes imageset="$1"-s.fits imagebinning=binSize xcolumn=X ximagebinsize=$binsize ycolumn=Y yimagebinsize=$binsize filtertype=expression expression='(FLAG == 0)&&(PI in [300:2000])'

# Creating an image for hard X-ray detections
evselect table="$1"_filt_time.fits withimageset=yes imageset="$1"-h.fits imagebinning=binSize xcolumn=X ximagebinsize=$binsize ycolumn=Y yimagebinsize=$binsize filtertype=expression expression='(FLAG == 0)&&(PI in [2000:10000])'

# Creating an image for all X-ray detections
evselect table="$1"_filt_time.fits withimageset=yes imageset="$1"-all.fits imagebinning=binSize xcolumn=X ximagebinsize=$binsize ycolumn=Y yimagebinsize=$binsize filtertype=expression expression='(FLAG == 0)&&(PI in [300:10000])'

# Creating an exposure map
eexpmap imageset="$1"-h.fits attitudeset=attitude.fits eventset="$1"_filt_time.fits expimageset="$1"_expmap_2-10keV.fits pimin=2000 pimax=10000


