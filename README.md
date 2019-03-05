# README #

### What is this repository for? ###

This is a set of programs to simulate large-footprint full-waveform lidar from airborne lidar and to process it and perform various other tasks. The parts are:


**gediRat**: simulates GEDI waveforms from ALS .las files and outputs ASCII or HDF5 waveforms.

gediMetric: processes full-waveform data (LVIS or simulated GEDI) and outputs metrics.

mapLidar: produces geotiffs from .las files of different properties. Also produces the bounding boxes of .las files for use in overlapLasFiles.csh.

lasPoints: outputs .pts files from .las files for selected areas.

lvisBullseye: produces the correlation bullseye plots of ALS to full-waveform data following Blair and Hofton (1999).

addNoiseHDF: Reads waveform data from HDF5 files and adds noise of a chosen level.


To find options, type the above command with "-help". These are explained in full detail later. The other .c files are either small test programs in the development of the above or are libraries called by the above. Important libraries are:


gediIO.c: contains the GEDI simulator.


There are some C and bash shells to control the above:

gediRatList.csh: batch processes gediRat.

filtForR.csh: converts the output .txt files into .csv files for reading in to R.

overlapLasFiles.csh: determines which .las files are needed for a given simulation.

orbitTracks.bash: produces lists of footprints from GEDI orbital simulations.

listALS.csh: produces the lists of .las files needed to read multiple files.



### How do I get set up? ###

Clone the repository and point to its location with the environment variable:

GEDIRAT_ROOT


All programs depend on these libraries:

Gnu Scientific Library

Geotiff

HDF5

GDAL

Minpack (Levenberg-Maruqardt):  https://www.physics.wisc.edu/~craigm/idl/cmpfit.html


Point to these with the environment variable:

GSL_ROOT

HDF5_LIB

CMPFIT_ROOT


Once they're installed, it also requires two libraries without package managers:

C-tools: https://bitbucket.org/StevenHancock/tools

libClidar: https://bitbucket.org/StevenHancock/libclidar

Clone these from bitbucket and point to their locations with the environment variables:

HANCOCKTOOLS_ROOT

LIBCLIDAR_ROOT

To compile, type:

make THIS=gediRat

Make sure that ~/bin/$ARCH exists and is in the path, where $ARCH is the result of `uname -m`.

make THIS=gediRat install

Replace "gediRat" with each of the commands above to compile and install.


Make sure that all .csh and .bash files are also in your path.

# Function operation #

## gediRat ##

Program to create GEDI waveforms from ALS las or pts files. laz not yet supported

### Input output filenames and format
-input name;     lasfile input filename

-inList list;    input file list (ASCII file) for multiple files

-output name;    output filename

-ground;         record separate ground and canopy waveforms

-hdf;            write output as HDF5. Best with gridded or list of coords

-ascii;          write output as ASCII (default). Good for quick tests

-waveID id;      supply a waveID to pass to the output (only for single footprints)

### Single footprint, list of footprints, or grid of footprints
-coord lon lat;  footprint coordinate in same system as lasfile

-listCoord name; list of coordinates

-gridBound minX maxX minY maxY;     make a grid of waveforms in this box

-gridStep res;   grid step size

### Lidar characteristics. Defaults are expected GEDI values.-pSigma sig;     set Gaussian pulse width as 1 sigma
-pFWHM fhwm;     set Gaussian pulse width as FWHM in ns

-readPulse file; read pulse shape and width from a file insteda of making Gaussian

-fSigma sig;     set footprint width

-wavefront file; read wavefront shape from file instead of setting Gaussian. Note that footprint width is still set by fSigma

-res res;        range resolution of waveform digitisation to output, in units of ALS data

-LVIS;           use LVIS pulse length, sigma=6.25m

-topHat;         use a top hat wavefront

-sideLobe;       use side lobes

-lobeAng ang;    lobe axis azimuth


### Input data quality filters
-checkCover;     check that the footprint is covered by ALS data. Do not output if not

-maxScanAng ang; maximum scan angle, degrees

-decimate x;     probability of accepting an ALS beam


### Computational speed options
-pBuff s;        point reading buffer size in Gbytes

-maxBins;        Optional: for HDF5, limit number of bins to save trimming.

-countOnly;      only use count method

-pulseAfter;     apply the pulse smoothing after binning for computational speed, at the risk of aliasing (default)

-pulseBefore;    apply the pulse smoothing before binning to avoid the risk of aliasing, at the expense of computational speed

-noNorm;         don't normalise for ALS density

### Octree
-noOctree;       do not use an octree

-octLevels n;    number of octree levels to use

-nOctPix n;      number of octree pixels along a side for the top level


### Using full-waveform input data (not tested)
-decon;          deconvolve

-indDecon;       deconvolve individual beams

-readWave;       read full-waveform where available


### Miscellaneous
-listFiles;      list files. Do not read them

-keepOld;        do not overwrite old files, if they exist

-useShadow;      account for shadowing in discrete return data through voxelisation

-polyGround;     find mean ground elevation and slope through fitting a polynomial

-nnGround;       find mean ground elevation and slope through nearest neighbour



### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

svenhancock@gmail.com
