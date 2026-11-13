# How to use the XMM-Newton Pipeline

## Preface

In the `xmm_data_reduction` directory there are 5 files for reducing data from the XMM-Newton observatory. The code is currently only targeted at data taken in small window mode with a single source. The function of these is laid out in the following structure

```
XMM_Pipeline.sh
    ├── Pipelines.sh
    |   └── GTIFilter.py
    └── Extraction.sh
        └── FindSource.py
```

where 

```
a
└── b
```

means that file `a` calls file `b`. Each file completes the following task(s):

1. `XMM_Pipeline.sh`
    - The main script
    - Calls all of the other files
    - Must be run within the `ODF` directory
2. `Pipelines.sh`
    - Runs all of the XMM-SAS pipeline scripts
    - Calls `GTIFilter.py`
3. `GTIFilter.py`
    - Determines the GTI filtering expression (output as a `.txt` file)
4. `Extraction.sh` 
    - Handles the generation of images, spectra, and response files
    - Calls `FindSource.py`
5. `FindSource.py`
    - Determines the source/background regions
    - The algorithm for the region selection for the pn instruments is described [here](https://bradadday.github.io/ProjectIronLine_Labbook/Week8.html).

## Usage

All of the scripts produce assume access to a HEASoft installation for the data reduction, as well as a python installation containing the necessary dependencies. The python dependencies can be installed with pip.

```
pip install astropy numpy matplotlib photutils
```

The main pipeline script can be run as follows. Note that it must be ran in the ODF directory of the XMM data.

```
Usage: /PATH/TO/SCRIPT/XMM_Pipeline.sh [-o OBS_ID] [-i INSTRUMENT] [-s SIGMA] [--source-max SOURCE_MAX] [--m1 MARGIN1] [--m2 MARGIN2]

    Options:
        -o, --obs-id:       Observation ID
        -i, --instrument:   Instrument
        -s:                 GTI standard deviation
        --source-max:       Maximum radius of source region
        --m1, --margin1:    Margin from primary flare of PSF
        --m2, --margin2:    Margin from edge of window
```

here

- `OBS_ID` is the observation ID of the data being reduced
- `INSTRUMENT` is the instrument to use for data processing
    - Either `pn`, `mos1` or `mos2`
    - Default is `pn`
- `SIGMA` is the number of standard deviations from the mean to accept when filtering for the Good Time Interval
    Default is `1.5`
- `SOURCE_MAX` is the maximum radius in physical units (detector pixels) to use
    - Default is `1000`
- `MARGIN1` is the distance the background regions must be separated from the primary flare of the source point spread function
    - The default is `5`
- `MARGIN2` is the distance all regions must be from the edge of the small window
    - The default is `1`

This runs all of the data reduction for a single instrument. To process another instrument, you can enter the `PROC` directory produced by the main pipeline (it will be in the same directory as the `ODF` directory) and run the `Extraction.sh` script as follows.

```
Usage: /home/jq23394/Data/Scripts/new_pipeline/Extraction.sh [-o OBS_ID] [-i INSTRUMENT] [--source-max SOURCE_MAX] [--m1 MARGIN1] [--m2 MARGIN2]

    Options:
        -o, --obs-id:       Observation ID
        -i, --instrument:   Instrument
        --source-max:       Maximum radius of source region
        --m1, --margin1:    Margin from primary flare of PSF
        --m2, --margin2:    Margin from edge of window
```

Here the arguments are the same as before.