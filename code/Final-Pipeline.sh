#!/bin/bash

# $1 is OBS_ID $2 is Instrument

# Setting environment variables
export SAS_ODFPATH="../ODF"
export SAS_CCFPATH="/opt/local/XMM/ccf/"
export SAS_CCF="$SAS_ODFPATH/ccf.cif"
export SAS_ODF="$SAS_ODFPATH/$(ls *SUM.SAS)"

# Find the source and background region areas for the ancillary file
backscale spectrumset="$2"_pi.fits badpixlocation="$2"_filt_time.fits
backscale spectrumset=bkg_pi.fits badpixlocation="$2"_filt_time.fits

# Generate the Photon Redistribution Matrix
rmfgen rmfset="$2"_rmf.fits spectrumset="$2"_pi.fits > rmflog.txt

# Generate the Ancillary File
arfgen arfset="$2"_arf.fits spectrumset="$2"_pi.fits withrmfset=yes rmfset="$2"_rmf.fits withbadpixcorr=yes badpixlocation="$2"_filt_time.fits > arflog.txt

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
