

/*##############################*/
/*# Structures for storing     #*/
/*# simulated GEDI waveforms   #*/
/*# 2017 svenhancock@gmail.com #*/
/*##############################*/

/*#######################################*/
/*# Copyright 2015-2017, Steven Hancock #*/
/*# The program is distributed under    #*/
/*# the terms of the GNU General Public #*/
/*# License.    svenhancock@gmail.com   #*/
/*#######################################*/


/*########################################################################*/
/*# This file is part of the NASA GEDI simulator, gediRat.               #*/
/*#                                                                      #*/
/*# gediRat is free software: you can redistribute it and/or modify      #*/
/*# it under the terms of the GNU General Public License as published by #*/
/*# the Free Software Foundation, either version 3 of the License, or    #*/
/*#  (at your option) any later version.                                 #*/
/*#                                                                      #*/
/*# gediRat is distributed in the hope that it will be useful,           #*/
/*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
/*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
/*#   GNU General Public License for more details.                       #*/
/*#                                                                      #*/
/*#    You should have received a copy of the GNU General Public License #*/
/*#    along with gediRat.  If not, see <http://www.gnu.org/licenses/>.  #*/
/*########################################################################*/


/*###########################################################*/
/*data structure*/

typedef struct{
  int nBins;        /*number of wavefor bins*/
  float *wave;      /*original waveform*/
  float *ground;    /*ground wave if we are to read it*/
  float *noised;    /*noised waveform, if needed*/
  double *z;        /*wave bin elevation*/
  double gElev;     /*mean ground elevation*/
  double tElev;     /*top elevation*/
  double lon;       /*footprint centre longitude*/
  double lat;       /*footprint centre latitude*/
  float totE;       /*total waveform energy (not integral)*/
  float gStdev;     /*measure of ground width*/
  float slope;      /*ground effective slope*/
  float cov;        /*ALS canopy cover*/
  float gLap;       /*ground overlap with canopy*/
  float gMinimum;   /*amplitude of minimum between ground and canopy*/
  float gInfl;      /*amplitude of inflection point between ground and canopy*/
  float pSigma;     /*pulse length*/
  float fSigma;     /*footprint width*/
  char useID;       /*use waveID or not*/
  char waveID[200]; /*wave ID for labelling*/
  char usable;      /*is the data usable*/
  float beamDense;  /*beam density*/
  float pointDense; /*point denisty*/
  float res;        /*range resolution*/
  float zen;        /*beam zenith angle (degrees)*/
  char demGround;   /*use the defined DEM ground values*/
  /*for LVIS HDF5 and level2*/
   uint32_t lfid;   /*LVIS file identifier*/
   uint32_t shotN;  /*LVIS shotnumber*/
}dataStruct;


/*###########################################################*/
/*GEDI IO options*/

typedef struct{
  /*input files*/
  int nFiles;   /*number of waveforms*/
  char **inList;

  /*switches*/
  char ground;      /*read separateground wave or not*/
  char useInt;      /*use discrete intensity instead of count*/
  char useFrac;     /*use fraction of hits per beam for weighting*/
  char dontTrustGround; /*don't trust ground included with waveforms*/

  /*denoising parameters*/
  denPar *den;   /*for denoising*/
  denPar *gFit;  /*for Gaussian fitting*/

  /*pulse parameters*/
  float pSigma;    /*pulse length*/
  float fSigma;    /*footprint width*/
  float res;      /*range resolution*/

  /*others*/
  int nMessages;  /*number of progress messages*/
}gediIOstruct;


/*###########################################################*/
/*functions*/

dataStruct **tidyAsciiStruct(dataStruct **,int);
dataStruct *readASCIIdata(char *,gediIOstruct *);
gediHDF *arrangeGEDIhdf(dataStruct **,gediIOstruct *);

/*# the end*/
/*###########################################################*/

