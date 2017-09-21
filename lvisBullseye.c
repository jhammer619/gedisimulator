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
#include "ogr_srs_api.h"


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

  /**/
  int nLvis;           /*number of LVIS*/

  /*bounds for subsets*/
  char useBounds;
  double minX;
  double maxX;
  double minY;
  double maxY;
  int aEPSG;           /*ALS epsg*/
  int lEPSG;           /*LVIS epsg*/
}control;


/*########################################################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct **lvis=NULL;
  dataStruct **readMultiLVIS(control *);


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

dataStruct **readMultiLVIS(control *dimage)
{
  int i=0;
  dataStruct **lvis=NULL;
  dataStruct **copyLVIShdf(lvisHDF *,dataStruct **,control *,double *);
  dataStruct **copyLVISlgw(char *,dataStruct **,control *,double *);
  double bounds[4];
  lvisHDF *hdf=NULL;
  void reprojectBounds(control *,double *);


  dimage->nLvis=0;

  /*reproject bounds*/
  reprojectBounds(dimage,bounds);

  /*loop over lvus files*/
  for(i=0;i<dimage->lvisIO.nFiles;i++){
    if(dimage->useLvisHDF){
      /*read*/
      hdf=readLVIShdf(dimage->lvisIO.inList[i]);
      /*unpack*/
      lvis=copyLVIShdf(hdf,lvis,dimage,bounds);
      /*tidy up*/
      hdf=tidyLVISstruct(hdf);
    }else if(dimage->useLvisLGW){
      /*unpack*/
      lvis=copyLVISlgw(dimage->lvisIO.inList[i],lvis,dimage,bounds);
    }
  }/*file loop*/

fprintf(stdout,"Found %d\n",dimage->nLvis);

  return(lvis);
}/*readMultiLVIS*/


/*####################################################*/
/*reproject the bounds*/

void reprojectBounds(control *dimage,double *bounds)
{
  double *x=NULL,*y=NULL,*z=NULL;
  OGRCoordinateTransformationH hTransform;
  OGRSpatialReferenceH hSourceSRS,hTargetSRS;
  OGRErr err;

  if(dimage->aEPSG!=dimage->lEPSG){
    x=dalloc(2,"x trans",9);
    y=dalloc(2,"y trans",9);
    z=dalloc(2,"z trans",9);

    x[0]=dimage->minX;
    x[1]=dimage->maxX;
    y[0]=dimage->minY;
    y[1]=dimage->maxY;
    z[0]=z[1]=0.0;

    hSourceSRS=OSRNewSpatialReference(NULL);
    hTargetSRS=OSRNewSpatialReference(NULL);
    err=OSRImportFromEPSG(hTargetSRS,dimage->aEPSG);
    err=OSRImportFromEPSG(hSourceSRS,dimage->lEPSG);
    hTransform=OCTNewCoordinateTransformation(hSourceSRS,hTargetSRS);
    OCTTransform(hTransform,2,x,y,z);
    OCTDestroyCoordinateTransformation(hTransform);
    OSRDestroySpatialReference(hSourceSRS);
    OSRDestroySpatialReference(hTargetSRS);

    bounds[0]=x[0];
    bounds[1]=y[0];
    bounds[2]=x[1];
    bounds[3]=y[1];
  }else{
    bounds[0]=dimage->minX;
    bounds[1]=dimage->minY;
    bounds[2]=dimage->maxX;
    bounds[3]=dimage->maxY;
  }

  TIDY(x);
  TIDY(y);
  TIDY(z);
  return;
}/*reprojectBounds*/


/*####################################################*/
/*unpack relevant data from LGW*/

dataStruct **copyLVISlgw(char *namen,dataStruct **lvis,control *dimage,double *bounds)
{
  int i=0,nNew=0;
  double x=0,y=0;
  lvisLGWstruct lgw;
  dataStruct *tempData=NULL;

  /*count number of new files*/
  nNew=0;
  tempData=readBinaryLVIS(namen,&lgw,0,&dimage->lvisIO);
  for(i=0;i<lgw.nWaves;i++){
    x=(lgw.data[i].lon0+lgw.data[i].lon431)/2.0;
    y=(lgw.data[i].lat0+lgw.data[i].lat431)/2.0;

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3]))nNew++;
  }

  /*allocate space*/
  if(lvis==NULL){
    if(!(lvis=(dataStruct **)calloc(nNew,sizeof(dataStruct *)))){
      fprintf(stderr,"error in lvis data allocation.\n");
      exit(1);
    }
  }else{
    if(!(lvis=(dataStruct **)realloc(lvis,(dimage->nLvis+nNew)*sizeof(dataStruct *)))){
      fprintf(stderr,"Balls\n");
      exit(1);
    }
  }

  /*copy data*/
  nNew=0;
  for(i=0;i<lgw.nWaves;i++){
    x=(lgw.data[i].lon0+lgw.data[i].lon431)/2.0;
    y=(lgw.data[i].lat0+lgw.data[i].lat431)/2.0;

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3])){
      lvis[nNew+dimage->nLvis]=readBinaryLVIS(namen,&lgw,i,&dimage->lvisIO);
      nNew++;
    }
  }

  /*keep count*/
  dimage->nLvis+=nNew;
  TIDY(lgw.data);

  return(lvis);
}/*copyLVISlgw*/


/*####################################################*/
/*unpack relevant data from HDF*/

dataStruct **copyLVIShdf(lvisHDF *hdf,dataStruct **lvis,control *dimage,double *bounds)
{
  int i=0,nNew=0;
  double x=0,y=0;

  /*count number within*/
  nNew=0;
  for(i=0;i<hdf->nWaves;i++){
    x=(hdf->lon0[i]+hdf->lon1023[i])/2.0;
    y=(hdf->lat0[i]+hdf->lat1023[i])/2.0;

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3]))nNew++;
  }/*wave loop*/
  if(!nNew)return(lvis);

  /*allocate space*/
  if(lvis==NULL){
    if(!(lvis=(dataStruct **)calloc(nNew,sizeof(dataStruct *)))){
      fprintf(stderr,"error in lvis data allocation.\n");
      exit(1);
    }
  }else{
    if(!(lvis=(dataStruct **)realloc(lvis,(dimage->nLvis+nNew)*sizeof(dataStruct *)))){
      fprintf(stderr,"Balls\n");
      exit(1);
    }
  }

  /*copy data*/
  nNew=0;
  for(i=0;i<hdf->nWaves;i++){
    x=(hdf->lon0[i]+hdf->lon1023[i])/2.0;
    y=(hdf->lat0[i]+hdf->lat1023[i])/2.0;

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3])){
      lvis[nNew+dimage->nLvis]=unpackHDFlvis(NULL,hdf,&dimage->lvisIO,i);
      nNew++;
    }
  }

  /*keep count*/
  dimage->nLvis+=nNew;

  return(lvis);
}/*copyLVIShdf*/


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
  dimage->aEPSG=32632;
  dimage->lEPSG=4326;



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
      }else if(!strncasecmp(argv[i],"-aEPSG",6)){
        checkArguments(1,i,argc,"-aEPSG");
        dimage->aEPSG=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-lEPSG",6)){
        checkArguments(1,i,argc,"-lEPSG");
        dimage->lEPSG=atoi(argv[++i]);
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

