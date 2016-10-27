# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

This is a set of programs to simulate large-footprint full-waveform lidar from airborne lidar and to process it and perform various other tasks.




### How do I get set up? ###

There is an inelegant setup script to install the simulator and dependencies by brute force. Do not use if you are happy installing on unix yourself. Script:

https://bitbucket.org/umdgedi/gedisimulator/src/2d34db078333b56f42d8c28ef95bacdcb2f23455/setupGediRat


It depends on these libraries:

Gnu Scientific Library

Minpack (Levenberg-Maruqardt)

Geotiff

Once they're installed, it also requires two less common, "special", libraries:

C-tools: https://bitbucket.org/StevenHancock/tools

libClidar: https://bitbucket.org/StevenHancock/libclidar


Once these have been downloaded, modify the makefile to point to these special libraries. Make sure that the common libraries are in the LD_LIBRARY_PATH and type "make". Add "THIS=gediMetric" to compile the other .c files in that directory.

There are some Cshells to control jobs. These point to awk scripts and so the "set bin=..." lines should be modified to point to the directory in which this simulator is downloaded.


### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

svenhancock@gmail.com