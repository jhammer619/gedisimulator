

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


/*####################################*/
/*pulse structure*/

typedef struct{
  int nBins;
  int centBin;  /*peak bin*/
  float *y;
  float *x;
}pulseStruct;


/*####################################*/
/*wavefront structure*/

typedef struct{
  int nX;         /*number of x pixels*/
  int nY;         /*number of y pixels*/
  int x0;         /*centre coordinate*/
  int y0;         /*centre coordinate*/
  float res;      /*resolution of wavefront grid*/
  float **front;  /*wavefront array*/
  char frontFile[200];/*file containing wavefront*/
}wFrontStruct;
               

/*####################################*/
/*lobe structure*/

typedef struct{
  double coord[2];  /*central coordinate*/
  float E;          /*fraction of energy here*/
  float fSigma;     /*footprint sigma*/
  double maxSepSq;  /*maximum distance from footprint needed*/
}lobeStruct;


/*###########################################################*/
/*data structure*/

typedef struct{
  int nBins;        /*number of wavefor bins*/
  int nWaveTypes;   /*number of types of waves, int, count, frac*/
  int useType;      /*the index to use for all analysis*/
  float **wave;     /*original waveform*/
  float **ground;   /*ground wave if we are to read it*/
  float *noised;    /*noised waveform, if needed*/
  double *z;        /*wave bin elevation*/
  double gElev;     /*mean ground elevation*/
  double tElev;     /*top elevation*/
  double lon;       /*footprint centre longitude*/
  double lat;       /*footprint centre latitude*/
  float *totE;      /*total waveform energy (not integral) per type*/
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
  char useInt;      /*use discrete intensity for weighting*/
  char useCount;    /*use no weighting*/
  char useFrac;     /*use fraction of hits per beam for weighting*/
  char dontTrustGround; /*don't trust ground included with waveforms*/
  char readPsigma;   /*read psigma from files or not*/

  /*denoising parameters*/
  denPar *den;   /*for denoising*/
  denPar *gFit;  /*for Gaussian fitting*/

  /*lidar parameters*/
  float pFWHM;     /*pulse width in ns*/
  float pSigma;    /*pulse length*/
  float fSigma;    /*footprint width*/
  float res;      /*range resolution*/

  /*system pulse*/
  char readPulse;      /*read pulse to simulate with*/
  char pulseFile[200];
  pulseStruct *pulse;
  float pRes;

  /*others*/
  int nMessages;  /*number of progress messages*/
}gediIOstruct;


/*###########################################################*/
/*GEDI simulator options*/

typedef struct{
  /*switches*/
  char readALSonce;    /*read all ALS data once*/
  char readWave;       /*read waveform switch*/
  char useShadow;      /*account for shadowing through voxelisation*/
  char topHat;       /*use a top hat wavefront rather than Gaussian*/
  char normCover;      /*normalise for variable ALS coverage*/
  char checkCover;     /*check that the whole footprit is covered by data*/
  char useFootprint;   /*use footprint or not flag*/
  char cleanOut;       /*clean subterranean outliers*/

  /*coordinates*/
  double coord[2];

  /*read a batch of coords*/
  char coordList[200]; /*list of coordinates*/
  char **waveIDlist;   /*list of waveform IDs*/
  double **coords;     /*list of coordinates*/

  /*GEDI footprint parameters, per footprint*/
  char sideLobe;     /*side lobe switch*/
  float lobeAng;     /*lobe major axis, degrees*/
  int nLobes;        /*number of side lobes*/
  lobeStruct *lobe;  /*lobe structure*/
  double minX;       /*minimum latitude of interest*/
  double maxX;       /*maximum latitude of interest*/
  double minY;       /*minimum longitude of interest*/
  double maxY;       /*maximum longitude of interest*/
  /*non-Gaussian wavefront if needed*/
  char defWfront;            /*define wavefront switch*/
  wFrontStruct *wavefront;   /*wavefront structure*/

  /*global area of interest*/
  double globMinX;
  double globMaxX;
  double globMinY;
  double globMaxY;

  /*grid to output multiple waveforms per run*/
  char doGrid;         /*gridded switch*/
  double gRes;         /*grid resolution*/
  double gMinX;        /*minimum x of grid*/
  double gMaxX;        /*maximum x of grid*/
  double gMinY;        /*minimum y of grid*/
  double gMaxY;        /*maximum y of grid*/
  int gNx;             /*number of x steps*/
  int gNy;             /*number of y steps*/

  /*grid to normalise sampling density*/
  int *nGrid;    /*beam per grid cell*/
  int gX;
  int gY;
  float gridRes;
  double g0[2];  /*grid origin*/

  /*point and beam density*/
  float denseRadSq; /*radius to calculate density within*/
  float pointDense; /*point density within 2 sigma*/
  float beamDense;  /*beam density within 2 sigma*/

  /*waveform processing*/
  char doDecon;    /*deconolution switch*/
  char indDecon;   /*deconolve individual ALS waves*/
  float meanN;
  denPar *decon;   /*denoising parameters*/

  /*others*/
  double maxSep;   /*maximum acceptable separation*/
  float maxScanAng;    /*maximum scan angle*/
  float iThresh;   /*intensity threshold*/

  /*for voxel shadows*/
  float beamRad;   /*beam radius at ground*/
  float vRes[3];   /*resolution along each axis*/
}gediRatStruct;


/*###########################################################*/
/*GEDI HDF5 structure*/

typedef struct{
  /*header*/
  int nWaves;      /*number of waveforms*/
  int nBins;       /*number of waveform bins*/
  int idLength;    /*length of wavefor ID strings*/
  float pSigma;    /*pulse length*/
  float fSigma;    /*footprint width*/
  float *pulse;    /*pulse*/
  float pRes;      /*pulse resolution*/
  int nPbins;      /*number ofpulse bins*/
  /*beams*/
  int nTypeWaves;  /*number of waveform types (frac, count and int)*/
  float **wave;    /*waveform*/
  float **ground;  /*ground waveforms*/
  float *z0;       /*wave top elevations*/
  float *zN;       /*wave bottom elevations*/
  double *lon;     /*longitudes*/
  double *lat;     /*latitudes*/
  float *slope;    /*ground slope*/
  float *gElev;    /*ground elevation, CofG*/
  float *demElev;  /*ground elevation, DEM*/
  char *waveID;    /*waveform ID*/
  float *beamDense;/*beam density*/
  float *pointDense;/*point density*/
  float *zen;      /*scan angles, or mean angles*/
}gediHDF;


/*####################################*/
/*waveform structure*/

typedef struct{
  float **wave;  /*waveforms*/
  float **canopy;/*canopy waveform*/
  float **ground;/*ground waveform*/
  double gElev;   /*ground elevation if calculated*/
  float gSlope;    /*ground sope if calculated*/
  double gElevSimp; /*simple ground elevation if calculated*/
  float gSlopeSimp;  /*simple ground sope if calculated*/
  float meanScanAng;/*mean ALS scan angle*/
  double minZ;    /*elevation bounds*/
  double maxZ;   /*elevation bounds*/
  int nBins;     /*number of wave bins*/
  int nWaves;    /*number of different ways*/
  double groundBreakElev;  /*break in ground*/
}waveStruct;


/*###########################################################*/
/*functions*/

dataStruct **tidyAsciiStruct(dataStruct **,int);
dataStruct *readASCIIdata(char *,gediIOstruct *);
dataStruct *unpackHDFlvis(char *,lvisHDF **,gediIOstruct *,int);
dataStruct *readBinaryLVIS(char *,lvisLGWstruct *,int,gediIOstruct *);
gediHDF *arrangeGEDIhdf(dataStruct **,gediIOstruct *);
gediHDF *readGediHDF(char *,gediIOstruct *);
gediHDF *tidyGediHDF(gediHDF *);
pCloudStruct *readALSdata(lasFile *las,gediRatStruct *gediRat);
waveStruct *makeGediWaves(gediRatStruct *,gediIOstruct *,pCloudStruct **);
void setGediGrid(gediIOstruct *,gediRatStruct *);
void setGediPulse(gediIOstruct *,gediRatStruct *);
void writeGEDIhdf(gediHDF *,char *);
void setGediFootprint(gediRatStruct *,gediIOstruct *);
void updateGediCoord(gediRatStruct *,int,int);
wFrontStruct *copyFrontFilename(char *);

/*# the end*/
/*###########################################################*/
