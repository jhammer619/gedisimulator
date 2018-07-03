# README #

### What is this repository for? ###

This is a set of programs to simulate large-footprint full-waveform lidar from airborne lidar and to process it and perform various other tasks. The parts are:


gediRat: simulates GEDI waveforms from ALS .las files and outputs ASCII or HDF5 waveforms.

gediMetric: processes full-waveform data (LVIS or simulated GEDI) and outputs metrics.

mapLidar: produces geotiffs from .las files of different properties. Also produces the bounding boxes of .las files for use in overlapLasFiles.csh.

lasPoints: outputs .pts files from .las files for selected areas.

lvisBullseye: produces the correlation bullseye plots of ALS to full-waveform data following Blair and Hofton (1999).

addNoiseHDF: Reads waveform data from HDF5 files and adds noise of a chosen level.


To find options, type the above command with "-help". These are exp,ained in full detail later. The other .c files are either small test programs in the development of the above or are libraries called by the above. Important libraries are:


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



### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

svenhancock@gmail.com
