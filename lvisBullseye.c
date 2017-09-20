#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "hdf5.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "gediIO.h"


/*#############################*/
/*# Uses simulated waveforms #*/
/*# to correlte to LVIS to   #*/
/*# coregister     2017      #*/
/*# svenhancock@gmail.com    #*/
/*############################*/

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
/*control structure*/

typedef struct{
  /*input/output*/
  gediIOstruct simIO; /*input/output structure*/
  gediIOstruct lvisIO; /*input/output structure*/
  char outNamen[200];

  /*options*/
  char useLvisHDF;
  char useLvisLGW;
  float offset;        /*vertical datum offset*/

  /*bullseye settings*/
  float maxShift;      /*maximum distance to shift*/
  float shiftStep;     /*distance to shift steps*/

  /*bounds for subsets*/
  char useBounds;
  double minX;
  double maxX;
  double minY;
  double maxY;
}control;


/*########################################################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct *lvis=NULL;
  dataStruct *readMultiLVIS(control *);


  /*read command Line*/
  dimage=readCommands(argc,argv);

  /*read LVIS data*/
  lvis=readMultiLVIS(dimage);


  /*read ALS data*/


  /*loop over, simulating*/

  /*tidy up*/
  if(dimage){
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################################*/
/*read multiple LVIS files and save relevant data*/

dataStruct *readMultiLVIS(control *dimage)
{
  int i=0;
  dataStruct *lvis=NULL;
  lvisHDF *hdf=NULL;

  /*allocate structures*/
  if(!(lvis=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  lvis->nWaveTypes=0;
  lvis->wave=NULL;
  lvis->ground=NULL;
  lvis->noised=NULL;
  lvis->z=NULL;


  for(i=0;i<dimage->lvisIO.nFiles;i++){
    if(dimage->useLvisHDF){
      /*read*/
      hdf=readLVIShdf(dimage->lvisIO.inList[i]);
      /*unpack*/

      /*tidy up*/
      hdf=tidyLVISstruct(hdf);
    }else if(dimage->useLvisLGW){

    }


    /*unpack relevant data*/
  }/*file loop*/



  return(lvis);
}/*readMultiLVIS*/


/*####################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  void setDenoiseDefault(denPar *);
  void readPulse(denPar *);

  /*allocate structures*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  if(!(dimage->simIO.den=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  dimage->simIO.gFit=NULL;
  if(!(dimage->lvisIO.den=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  dimage->lvisIO.gFit=NULL;


  /*defaults*/
  dimage->useLvisHDF=1;
  dimage->useLvisLGW=0;



  /*read the command line*/
  for(i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-lvis",5)){
        checkArguments(1,i,argc,"-lvis");
        TTIDY((void **)dimage->lvisIO.inList,dimage->lvisIO.nFiles);
        dimage->lvisIO.inList=NULL;
        dimage->lvisIO.nFiles=1;
        dimage->lvisIO.inList=chChalloc(dimage->lvisIO.nFiles,"input name list",0);
        dimage->lvisIO.inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->lvisIO.inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-lvisList",9)){
        checkArguments(1,i,argc,"-lvisList");
        TTIDY((void **)dimage->lvisIO.inList,dimage->lvisIO.nFiles);
        dimage->lvisIO.inList=readInList(&dimage->lvisIO.nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-als",4)){
        checkArguments(1,i,argc,"-als");
        TTIDY((void **)dimage->simIO.inList,dimage->simIO.nFiles);
        dimage->simIO.inList=NULL;
        dimage->simIO.nFiles=1;
        dimage->simIO.inList=chChalloc(dimage->simIO.nFiles,"input name list",0);
        dimage->simIO.inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->simIO.inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-alsList",8)){
        checkArguments(1,i,argc,"-alsList");
        TTIDY((void **)dimage->simIO.inList,dimage->simIO.nFiles);
        dimage->simIO.inList=readInList(&dimage->simIO.nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(4,i,argc,"-bounds");
        dimage->useBounds=1;
        dimage->minX=atof(argv[++i]);
        dimage->minY=atof(argv[++i]);
        dimage->maxX=atof(argv[++i]);
        dimage->maxY=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-lgw",4)){
        dimage->useLvisHDF=0;
        dimage->useLvisLGW=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to calculate GEDI waveform metrics\n#####\n\n-input name;     waveform  input filename\n-output name;   output filename\n-inList list;    input file list for multiple files\n-writeFit;       write fitted waveform\n-writeGauss;     write Gaussian parameters\n-ground;         read true ground from file\n-useInt;         use discrete intensity instead of count\n-useFrac;        use fractional hits rather than counts\n-readBinLVIS;    input is an LVIS binary file\n-readHDFlvis;    read LVIS HDF5 input\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry gediRat -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }

  return(dimage);
}/*readCommands*/

/*the end*/
/*########################################################################*/

