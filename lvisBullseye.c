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
#include "libOctree.h"
#include "gediIO.h"
#include "ogr_srs_api.h"
#include "gsl/gsl_multimin.h"



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
  gediIOstruct simIO;    /*input/output structure*/
  gediIOstruct lvisIO;   /*input/output structure*/
  gediRatStruct gediRat; /*simulator options*/
  char outNamen[200];    /*correlation output filename*/
  char waveNamen[200];   /*waveform output filename*/
  int nLvis;             /*number of LVIS*/

  /*options*/
  char useLvisHDF;     /*input data format switch*/
  char useLvisLGW;     /*input data format switch*/
  char useGediHDF;     /*input data format switch*/
  char writeFinWave;   /*write out final waveforms switch*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/
  char filtOutli;      /*filter outliers to avoid falling trees*/
  float maxZen;        /*maximum LVIS zenith angle to use*/
  float minDense;      /*minimum ALS beam density*/
  float minSense;      /*minimum waveform beam sensitivity to use*/
  char leaveEmpty;     /*exit if no usable footprints at any time*/
  FILE *opoo;          /*output file*/

  /*bullseye settings*/
  float maxShift;      /*maximum horizontal distance to shift*/
  float shiftStep;     /*distance to shift steps horizontally*/
  float maxVshift;     /*maximum vertical distance to shift*/
  float vShiftStep;    /*distance to shift steps vertically*/
  double origin[3];    /*origin to shift around*/

  /*bullseye optimisation*/
  char fullBull;       /*do the full bullseye. If not, do simplex*/
  char findFsig;       /*solve for fSigma*/
  int maxIter;         /*maximum number of interations*/
  float optTol;        /*tolerance*/
  int nUsed;           /*number of footprints used in optimum*/
  char writeSimProg;   /*write simplex progress switch*/

  /*for large geolocation errors*/
  char largeErr;       /*switch for large error method*/

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
/*structure containg data for optimisation*/

typedef struct{
  control *dimage;
  float **denoised;
  int nTypeWaves;
  double **coords;
  dataStruct **lvis;
  pCloudStruct **als;
}optStruct;


/*####################################*/
/*global functions*/

float **getCorrelStats(control *,dataStruct **,pCloudStruct **,int *,double,double,double,double **,float **,int,char);
float *waveCorrel(waveStruct *,float *,dataStruct *,gediIOstruct *,double);
double findMeanCorr(const gsl_vector *, void *);    /*error fuction for optimisation*/
double **shiftPrints(double **,double,double,int);
void fullBullseyePlot(control *,float **,int,dataStruct **,pCloudStruct **,float *);
void simplexBullseye(control *,float **,int,dataStruct **,pCloudStruct **);


/*########################################################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct **lvis=NULL;
  dataStruct **readMultiLVIS(control *,float *);
  pCloudStruct **als=NULL;
  pCloudStruct **readMultiALS(control *,dataStruct **);
  void bullseyeCorrel(dataStruct **,pCloudStruct **,control *);
  void setALSbounds(control *);

  /*read command Line*/
  dimage=readCommands(argc,argv);

  /*read ALS bounds if none defined*/
  if(dimage->useBounds==0)setALSbounds(dimage);

  /*read LVIS data*/
  lvis=readMultiLVIS(dimage,&dimage->simIO.res);

  /*read ALS data*/
  als=readMultiALS(dimage,lvis);

  /*loop over, simulating*/
  bullseyeCorrel(lvis,als,dimage);

  /*tidy up allocated memory*/
  if(als){
    for(i=0;i<dimage->simIO.nFiles;i++)als[i]=tidyPointCloud(als[i]);
    TIDY(als);
  }
  lvis=tidyAsciiStruct(lvis,dimage->nLvis);
  if(dimage){
    TTIDY((void **)dimage->lvisIO.inList,dimage->lvisIO.nFiles);
    TTIDY((void **)dimage->simIO.inList,dimage->simIO.nFiles);
    if(dimage->lvisIO.den){
      TTIDY((void **)dimage->lvisIO.den->pulse,2);
      TIDY(dimage->lvisIO.den->matchPulse);
      TIDY(dimage->lvisIO.den->hardPulse);
      TIDY(dimage->lvisIO.den);
    }
    if(dimage->lvisIO.gFit){
      TTIDY((void **)dimage->lvisIO.gFit->pulse,2);
      TIDY(dimage->lvisIO.gFit);
    }
    if(dimage->simIO.den){
      TTIDY((void **)dimage->simIO.den->pulse,2);
      TIDY(dimage->simIO.den->matchPulse);
      TIDY(dimage->simIO.den->hardPulse);
      TIDY(dimage->simIO.den);
    }
    if(dimage->simIO.gFit){
      TTIDY((void **)dimage->simIO.gFit->pulse,2);
      TIDY(dimage->simIO.gFit);
    }
    TTIDY((void **)dimage->gediRat.coords,dimage->gediRat.gNx);
    TTIDY((void **)dimage->gediRat.waveIDlist,dimage->gediRat.gNx);
    TIDY(dimage->gediRat.nGrid);
    dimage->gediRat.octree=tidyOctree(dimage->gediRat.octree);
    TIDY(dimage);
  }/*tidying up allocated memory*/
  return(0);
}/*main*/


/*####################################################*/
/*simulate and calculate correlation*/

void bullseyeCorrel(dataStruct **lvis,pCloudStruct **als,control *dimage)
{
  int nTypeWaves=0;
  float **denoiseAllLvis(dataStruct **,control *);
  void rapidGeolocation(control *,float **,int,dataStruct **,pCloudStruct **);
  float **denoised=NULL;

  /*how mamy types of simuation methods*/
  nTypeWaves=(int)(dimage->simIO.useCount+dimage->simIO.useFrac+dimage->simIO.useInt);

  /*set up pulse*/
  setGediPulse(&dimage->simIO,&dimage->gediRat);

  /*denoise LVIS*/
  denoised=denoiseAllLvis(lvis,dimage);

  /*are we doing the full bullseye plot?*/
  if(dimage->fullBull){
    fullBullseyePlot(dimage,denoised,nTypeWaves,lvis,als,NULL);
  }else if(dimage->largeErr){ /*do a rough bullseye followed by a simplex*/
    rapidGeolocation(dimage,denoised,nTypeWaves,lvis,als);
  }else{ /*do simplex*/
    simplexBullseye(dimage,denoised,nTypeWaves,lvis,als);
  }/*are we doing the full bullseye plot?*/


  /*tidy up*/
  TTIDY((void **)denoised,dimage->nLvis);
  if(dimage->simIO.pulse){
    TIDY(dimage->simIO.pulse->y);
    TIDY(dimage->simIO.pulse->x);
    TIDY(dimage->simIO.pulse);
  }
  return;
}/*bullseyeCorrel*/


/*####################################################*/
/*do a rough bullseye followed by a simplex*/

void rapidGeolocation(control *dimage,float **denoised,int nTypeWaves,dataStruct **lvis,pCloudStruct **als)
{
  int nX=0,nZ=0;
  float x=0,y=0,z=0;
  float *roughCorrel=NULL;   /*mean correlation for rough bullseye*/
  void bestRoughGeo(float *,float *,float *,float *,int,int,double *,float,float,control *);

  /*perform rough bullseye run*/
  /*allocate space*/
  nX=(int)(2.0*dimage->maxShift/dimage->shiftStep+1);
  if(dimage->maxVshift>0.0)nZ=(int)(2.0*dimage->maxVshift/dimage->vShiftStep+1);
  else                     nZ=1;
  roughCorrel=falloc((uint64_t)nX*(uint64_t)nX*(uint64_t)nZ,"mean correlation",0);
  fullBullseyePlot(dimage,denoised,nTypeWaves,lvis,als,roughCorrel);

  /*open output if needed*/
  if((dimage->opoo==NULL)&&(dimage->writeSimProg)){
    if((dimage->opoo=fopen(dimage->outNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
      exit(1);
    }
  }

  /*find optimium location*/
  bestRoughGeo(&x,&y,&z,roughCorrel,nX,nZ,dimage->origin,dimage->shiftStep,dimage->vShiftStep,dimage);
  TIDY(roughCorrel);

  /*start simplex from optimum*/
  dimage->origin[0]=(double)x;
  dimage->origin[1]=(double)y;
  dimage->origin[2]=(double)z;
  simplexBullseye(dimage,denoised,nTypeWaves,lvis,als);

  return;
}/*rapidGeolocation*/


/*####################################################*/
/*find optimium correlation shift*/

void bestRoughGeo(float *x,float *y,float *z,float *roughCorrel,int nX,int nZ,double *origin,float xStep,float zStep,control *dimage)
{
  int i=0,j=0,k=0,place=0;
  float maxCorrel=0;

  /*set nonesense value*/
  maxCorrel=-1000.0;

  /*loop over all correls*/
  for(i=0;i<nX;i++){
    for(j=0;j<nX;j++){
      for(k=0;k<nZ;k++){
        place=i*nX*nZ+j*nZ+k;
        if(dimage->writeSimProg){
          fprintf(dimage->opoo,"%f %f %f %f roughGrid\n",(float)(i-nX/2)*xStep+(float)origin[0],(float)(j-nX/2)*xStep+(float)origin[1],(float)(k-nZ/2)*zStep+(float)origin[2],roughCorrel[place]);
        }
        if(roughCorrel[place]>maxCorrel){
          maxCorrel=roughCorrel[place];
          *x=(float)(i-nX/2)*xStep+(float)origin[0];
          *y=(float)(j-nX/2)*xStep+(float)origin[1];
          *z=(float)(k-nZ/2)*zStep+(float)origin[2];
        }/*max value check*/
      }/*z loop*/
    }/*y loop*/
  }/*x loop*/

  fprintf(stdout,"Starting from %f %f %f\n",*x,*y,*z);
  return;
}/*bestRoughGeo*/


/*####################################################*/
/*use simplex to optimise geolocation*/

void simplexBullseye(control *dimage,float **denoised,int nTypeWaves,dataStruct **lvis,pCloudStruct **als)
{
  int i=0,nPar=0;
  int status;
  float fSig=0;
  double size;
  size_t iter=0;
  const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2;
  gsl_vector *start=NULL,*ss=NULL;    /*initial conditions*/
  gsl_multimin_fminimizer *s=NULL;
  gsl_multimin_function minex_func;  /*optimiser options*/
  optStruct optBits;         /*to pass needed pieces to optimiser*/
  void writeFinalWaves(control *,dataStruct **,pCloudStruct **,double **,float,float,float,float);


  /*load needed bits into structure*/
  optBits.dimage=dimage;
  optBits.denoised=denoised;
  optBits.coords=dimage->gediRat.coords;
  dimage->gediRat.coords=NULL;
  optBits.nTypeWaves=nTypeWaves;
  optBits.lvis=lvis;
  optBits.als=als;
  dimage->gediRat.coords=NULL;

  /*allocate space*/
  if(dimage->findFsig)nPar=4;
  else                nPar=3;
  start=gsl_vector_alloc(nPar);
  /*initial estimate*/
  gsl_vector_set(start,0,dimage->origin[0]);
  gsl_vector_set(start,1,dimage->origin[1]);
  gsl_vector_set(start,2,dimage->origin[2]);
  if(dimage->findFsig)gsl_vector_set(start,3,dimage->simIO.fSigma);
  minex_func.n=nPar;
  minex_func.f=findMeanCorr;   /*the opimisation function*/
  minex_func.params =(void *)(&optBits);

  /* Set initial step sizes to expected step over four*/
  ss=gsl_vector_alloc(nPar);
  gsl_vector_set_all(ss,dimage->shiftStep/4.0);
  if(dimage->findFsig)gsl_vector_set(ss,3,0.2);

  /*first run*/
  s=gsl_multimin_fminimizer_alloc(T,nPar);
  gsl_multimin_fminimizer_set(s,&minex_func,start,ss);
  /*iterate*/
  do{
    /*count iterations*/
    iter++;

    status=gsl_multimin_fminimizer_iterate(s);
    if(status)break;
    size=gsl_multimin_fminimizer_size (s);
    status=gsl_multimin_test_size(size,dimage->optTol);

    /*write progress if needed*/
    if(dimage->writeSimProg){
      for(i=0;i<nPar;i++)fprintf(dimage->opoo,"%f ",(float)gsl_vector_get (s->x,i));
      fprintf(dimage->opoo,"%f %d simplex\n",1.0-(float)s->fval,dimage->nUsed);
    }
  }while((status==GSL_CONTINUE)&&(iter<dimage->maxIter));

  /*write results*/
  /*open output*/
  if(dimage->nUsed>0){
    if(dimage->opoo==NULL){
      if((dimage->opoo=fopen(dimage->outNamen,"w"))==NULL){
        fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
        exit(1);
      }
    }
    for(i=0;i<nPar;i++)fprintf(dimage->opoo,"%f ",(float)gsl_vector_get (s->x,i));
    fprintf(dimage->opoo,"%f %d\n",1.0-(float)s->fval,dimage->nUsed);
    if(dimage->opoo){
      fclose(dimage->opoo);
      dimage->opoo=NULL;
    }
    fprintf(stdout,"Written to %s\n",dimage->outNamen);
  }

  /*output waveforms if needed*/
  if(dimage->writeFinWave){
    if(nPar==4)fSig=(float)gsl_vector_get (s->x,3);
    else       fSig=dimage->simIO.fSigma;
    writeFinalWaves(dimage,lvis,als,optBits.coords,(float)gsl_vector_get(s->x,0),(float)gsl_vector_get(s->x,1),(float)gsl_vector_get(s->x,2),fSig);
  }

  /*tidy up*/
  gsl_vector_free(start);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
  dimage->gediRat.coords=optBits.coords;
  optBits.coords=NULL;
  optBits.dimage=NULL;
  optBits.denoised=NULL;
  optBits.lvis=NULL;
  optBits.als=NULL;
  return;
}/*simplexBullseye*/


/*####################################################*/
/*Error function for optimisation*/

double findMeanCorr(const gsl_vector *v, void *params)
{
  int i=0,nTypeWaves=0,contN=0;
  float **denoised=NULL;
  float **correl=NULL,meanCorrel=0;
  double xOff=0,yOff=0,zOff=0;
  double **coords=NULL;
  control *dimage=NULL;
  dataStruct **lvis=NULL;
  pCloudStruct **als=NULL;
  optStruct *optBits=NULL;

  /*unpack structure*/
  optBits=(optStruct *)params;
  dimage=optBits->dimage;
  denoised=optBits->denoised;
  coords=optBits->coords;
  nTypeWaves=optBits->nTypeWaves;
  lvis=optBits->lvis;
  als=optBits->als;

  /*estimate*/
  xOff=gsl_vector_get(v,0);
  yOff=gsl_vector_get(v,1);
  zOff=gsl_vector_get(v,2);
  if(dimage->findFsig){
    dimage->simIO.fSigma=(float)gsl_vector_get(v,3);
    if(dimage->simIO.fSigma<1.0)dimage->simIO.fSigma=1.0;  /*prevent non-physical*/
  }

  /*get correlations*/
  correl=getCorrelStats(dimage,lvis,als,&contN,xOff,yOff,zOff,coords,denoised,nTypeWaves,dimage->leaveEmpty);

  /*get mean correlation*/
  meanCorrel=0.0;
  for(i=0;i<contN;i++)meanCorrel+=correl[i][0];
  if(contN>0)meanCorrel/=(float)contN;

  /*record number used*/
  dimage->nUsed=contN;

  /*tidy up*/
  TTIDY((void **)dimage->gediRat.coords,dimage->gediRat.gNx);
  TTIDY((void **)correl,dimage->gediRat.gNx);
  dimage=NULL;
  coords=NULL;
  denoised=NULL;
  lvis=NULL;
  als=NULL;
  optBits=NULL; 
  return(1.0-(double)meanCorrel);
}/*findMeanCorr*/


/*####################################################*/
/*write out finale waveform file*/

void writeFinalWaves(control *dimage,dataStruct **lvis,pCloudStruct **als,double **coords,float xOff,float yOff,float zOff,float fSig)
{
  int k=0;
  int hdfCount=0;
  waveStruct *waves=NULL;
  gediHDF *hdfData=NULL;
  char waveID[500];

  /*how mamy types of simuation methods*/
  dimage->simIO.useCount=1;
  dimage->simIO.useFrac=dimage->simIO.useInt=0;
  dimage->simIO.nTypeWaves=(int)(dimage->simIO.useCount+dimage->simIO.useFrac+dimage->simIO.useInt);
  dimage->simIO.ground=1;
  dimage->gediRat.doGrid=1;
  dimage->gediRat.waveIDlist=NULL;

  /*allocate space*/
  sprintf(waveID,"%s.100000000",lvis[0]->waveID);
  hdfData=setUpHDF(&dimage->simIO,&dimage->gediRat,1,waveID,&hdfCount,1024);

  /*set to optimum sigma*/
  dimage->simIO.fSigma=fSig;

  /*simulate waveforms*/
  for(k=0;k<dimage->gediRat.gNx;k++){
    if((lvis[k]->zen>dimage->maxZen)||(lvis[k]->beamSense<dimage->minSense))continue;
    /*set one coordinate*/
    dimage->gediRat.coord[0]=coords[k][0]+xOff;
    dimage->gediRat.coord[1]=coords[k][1]+yOff;

    /*simulate waves*/
    setGediFootprint(&dimage->gediRat,&dimage->simIO);
    waves=makeGediWaves(&dimage->gediRat,&dimage->simIO,als);

    /*save to HDF5 file*/
    if(dimage->gediRat.useFootprint){
      /*apply z offset*/
      waves->minZ-=zOff;
      waves->maxZ-=zOff;
      waves->gElev-=zOff;
      waves->gElevSimp-=zOff;
      /*load to HDF structure*/
      packGEDIhdf(waves,hdfData,k,&dimage->simIO,&dimage->gediRat,&hdfCount,1,lvis[k]->waveID);
    }

    /*tidy up*/
    if(waves){
      TTIDY((void **)waves->wave,waves->nWaves);
      TTIDY((void **)waves->canopy,dimage->simIO.nTypeWaves);
      TTIDY((void **)waves->ground,dimage->simIO.nTypeWaves);
      TIDY(waves);
    }
    TIDY(dimage->gediRat.lobe);
    TIDY(dimage->gediRat.nGrid);
  }/*waveform loop*/

  /*write HDF5 file*/
  hdfData->nWaves=hdfCount;  /*account for unusable footprints*/
  writeGEDIhdf(hdfData,dimage->waveNamen,&(dimage->simIO));

  return;
}/*writeFinalWaves*/

/*####################################################*/
/*loop over full bullseye plot and write*/

void fullBullseyePlot(control *dimage,float **denoised,int nTypeWaves,dataStruct **lvis,pCloudStruct **als,float *meanCorrel)
{
  int indOff=0,place=0;
  int i=0,j=0,m=0,k=0;
  int nX=0,nZ=0,contN=0;
  float **correl=NULL;
  double **coords=NULL;
  double xOff=0,yOff=0,zOff=0;
  void writeCorrelStats(float **,int,int,FILE *,double,double,double,control *);
  FILE *opoo=NULL;

  /*open output file and write header*/
  if(meanCorrel==NULL){
    if((opoo=fopen(dimage->outNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
      exit(1);
    }
    fprintf(opoo,"# 1 xOff, 2 yOff");
    if(dimage->maxVshift>0.0){
      fprintf(opoo,", 3 zOff");
      indOff=1;
    }else indOff=0;
    if(dimage->simIO.useInt)fprintf(opoo,", %d correlInt, %d stdev, %d deltaCofGint, %d numb",3+indOff,4+indOff,5+indOff,6+indOff);
    if(dimage->simIO.useCount)fprintf(opoo,", %d correlCount, %d stdev, %d deltaCofGcount, %d numb",3+4*dimage->simIO.useInt+indOff,3+4*dimage->simIO.useInt+1+indOff,3+4*dimage->simIO.useInt+2+indOff,3+4*dimage->simIO.useInt+3+indOff);
    if(dimage->simIO.useFrac)fprintf(opoo,", %d correlFrac, %d stdev, %d deltaCofGfrac, %d numb",3+4*(dimage->simIO.useInt+dimage->simIO.useCount)+indOff,3+4*(dimage->simIO.useInt+dimage->simIO.useCount)+1+indOff\
                                                                                                ,3+4*(dimage->simIO.useInt+dimage->simIO.useCount)+2+indOff,3+4*(dimage->simIO.useInt+dimage->simIO.useCount)+3+indOff);
    fprintf(opoo,"\n");
  }/*open output file and write header*/

  /*calculate number of steps*/
  nX=(int)(2.0*dimage->maxShift/dimage->shiftStep+1);
  if(dimage->maxVshift>0.0)nZ=(int)(2.0*dimage->maxVshift/dimage->vShiftStep+1);
  else                     nZ=1;

  /*save original coordinates*/
  coords=dimage->gediRat.coords;
  dimage->gediRat.coords=NULL;

  /*loop over shifts and get correltion*/
  for(i=0;i<nX;i++){  /*x loop*/
    xOff=(double)(i-nX/2)*(double)dimage->shiftStep+dimage->origin[0];
    for(j=0;j<nX;j++){  /*y loop*/
      yOff=(double)(j-nX/2)*(double)dimage->shiftStep+dimage->origin[1];
      for(m=0;m<nZ;m++){  /*z loop*/
        if(nZ>1)zOff=(double)(m-nZ/2)*(double)dimage->vShiftStep;   /*datum offset has already been applied*/
        else    zOff=dimage->origin[2];
        /*get correlation stats*/
        correl=getCorrelStats(dimage,lvis,als,&contN,xOff,yOff,zOff,coords,denoised,nTypeWaves,dimage->leaveEmpty);

        /*which mode are we using*/
        if(meanCorrel==NULL){  /*write out all points*/
          /*output results*/
          writeCorrelStats(correl,contN,nTypeWaves,opoo,xOff,yOff,zOff,dimage);
        }else{                /*rough cut; keep track of mean correlation*/
          place=i*nX*nZ+j*nZ+m;
          meanCorrel[place]=0.0;
          for(k=0;k<contN;k++)meanCorrel[place]+=correl[k][0];
          if(contN>0)meanCorrel[place]/=(float)contN;
        }/*which mode are we using*/

        /*tidy up*/
        TTIDY((void **)dimage->gediRat.coords,dimage->gediRat.gNx);
        TTIDY((void **)correl,dimage->gediRat.gNx);
      }/*z loop*/
    }/*y loop*/
  }/*x loop*/

  /*return original coords*/
  dimage->gediRat.coords=coords;
  coords=NULL;
  /*close output file*/
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  /*let them know output file, if made*/
  if(meanCorrel==NULL)fprintf(stdout,"Written to %s\n",dimage->outNamen);
  return;
}/*fullBullseyePlot*/


/*####################################################*/
/*get correlation across all waves*/

float **getCorrelStats(control *dimage,dataStruct **lvis,pCloudStruct **als,int *contN,double xOff,double yOff,double zOff,double **coords,float **denoised,int nTypeWaves,char leaveEmpty)
{
  int k=0;
  float **correl=NULL;
  waveStruct *waves=NULL;

  /*progress report*/
  fprintf(stdout,"Testing x %.2f y %.2f z %.2f fSig %.2f\n",xOff,yOff,zOff,dimage->simIO.fSigma);

  /*shift prints*/
  dimage->gediRat.coords=shiftPrints(coords,xOff,yOff,dimage->gediRat.gNx);
  correl=fFalloc(dimage->gediRat.gNx,"correl",0);

  /*loop over footprints*/
  *contN=0;
  for(k=0;k<dimage->gediRat.gNx;k++){
    if((lvis[k]->zen>dimage->maxZen)||(lvis[k]->beamSense<dimage->minSense))continue;
    /*set one coordinate*/
    updateGediCoord(&dimage->gediRat,k,0);

    /*simulate waves*/
    setGediFootprint(&dimage->gediRat,&dimage->simIO);
    waves=makeGediWaves(&dimage->gediRat,&dimage->simIO,als);

    /*calculate correlation*/
    if(dimage->gediRat.useFootprint&&(dimage->gediRat.beamDense>=dimage->minDense)){
      correl[*contN]=waveCorrel(waves,denoised[k],lvis[k],&dimage->simIO,zOff);
      (*contN)++;
    }//else correl[*contN]=NULL;

    /*tidy up*/
    if(waves){
      TTIDY((void **)waves->wave,waves->nWaves);
      TTIDY((void **)waves->canopy,nTypeWaves);
      TTIDY((void **)waves->ground,nTypeWaves);
      TIDY(waves);
    }
    TIDY(dimage->gediRat.lobe);
    TIDY(dimage->gediRat.nGrid);
  }/*footprint loop*/

  /*check that we have something*/
  if((*contN)==0){
    fprintf(stderr,"No usable footprints contained\n");
    TTIDY((void **)correl,dimage->gediRat.gNx);
    correl=NULL;
    if(leaveEmpty)exit(1);
  }
  return(correl);
}/*getCorrelStats*/


/*####################################################*/
/*correlation stats and write*/

void writeCorrelStats(float **correl,int numb,int nTypes,FILE *opoo,double xOff,double yOff,double zOff,control *dimage)
{
  int i=0,k=0,nUsed=0;
  int usedNew=0;
  float mean=0,stdev=0;
  float newMean=0,meanCofG=0;
  float thresh=0;
  char writtenCoords=0;

  /*loop over types*/
  for(k=0;k<nTypes;k++){
    mean=stdev=0.0;
    nUsed=0;
    for(i=0;i<numb;i++){
      if(correl[i]){
        mean+=correl[i][2*k];
        nUsed++;
      }
    }
    mean/=(float)nUsed;
    meanCofG/=(float)nUsed;
    for(i=0;i<numb;i++){
      if(correl[i]){
        stdev+=pow(correl[i][2*k]-mean,2.0);
      }
    }
    stdev=sqrt(stdev/(float)nUsed);

    /*check for outliers*/
    if(dimage->filtOutli)thresh=2.5*stdev;
    else                 thresh=1000000000.0;

    usedNew=0;
    newMean=meanCofG=0.0;
    for(i=0;i<numb;i++){
      if(correl[i]){
        if((mean-correl[i][2*k])<thresh){
          newMean+=correl[i][2*k];
          meanCofG+=correl[i][2*k+1];
          usedNew++;
        }
      }
    }

    /*check that there is data*/
    if(nUsed==0)return;

    if(writtenCoords==0){
      fprintf(opoo,"%f %f",xOff,yOff);
      if(dimage->maxVshift>0.0)fprintf(opoo," %f",zOff+dimage->origin[2]);  /*for backwards compatability*/
      writtenCoords=1;
    }

    newMean/=(float)usedNew;
    meanCofG/=(float)usedNew;
    stdev=0.0;
    for(i=0;i<numb;i++){
      if(correl[i]){
        if((mean-=correl[i][2*k])<thresh){
          stdev+=pow(correl[i][2*k]-newMean,2.0);
        }
      }
    }
    stdev=sqrt(stdev/(float)usedNew);
    fprintf(opoo," %f %f %f %d",newMean,stdev,meanCofG,usedNew);
  }
  fprintf(opoo," %d\n",nUsed);

  return;
}/*writeCorrelStats*/


/*####################################################*/
/*calculate correlation*/

float *waveCorrel(waveStruct *sim,float *truth,dataStruct *lvis,gediIOstruct *simIO,double zOff)
{
  int i=0,j=0,k=0,bin=0;
  int numb=0;
  float *correl=NULL;
  float totS=0,totL=0;
  float CofGl=0,CofGs=0;
  float z=0,thresh=0;
  float sLx=0,eLx=0;
  float sSx=0,eSx=0;
  float startX=0,endX=0;
  float sumProd=0,minSepSq=0;
  float sepSq=0;
  float *smooSim=NULL;
  float meanL=0,meanS=0;
  float stdevL=0,stdevS=0;
  float *zShift=NULL;

  /*allocate space fgor correlation and CofG shift*/
  correl=falloc(2*(uint64_t)sim->nWaves,"correlation",0);

  /*apply datum offset*/
  zShift=falloc((uint64_t)lvis->nBins,"shifted LVIS elevation",0);
  for(i=0;i<lvis->nBins;i++)zShift[i]=(float)lvis->z[i]+zOff;

  /*total energies*/
  totL=CofGl=0.0;
  for(i=0;i<lvis->nBins;i++){
    totL+=truth[i];
    CofGl+=truth[i]*zShift[i];
  }
  CofGl/=totL;

  /*find lvis bounds*/
  thresh=0.0001*totL;
  for(i=0;i<lvis->nBins;i++){
    if(truth[i]>thresh){
      eLx=zShift[i];
      break;
    }
  }
  for(i=lvis->nBins-1;i>=0;i--){
    if(truth[i]>thresh){
      sLx=zShift[i];
      break;
    }
  }

  /*loop over three wave types*/
  for(k=0;k<sim->nWaves;k++){
    smooSim=processFloWave(sim->wave[k],sim->nBins,simIO->den,1.0);


    totS=CofGs=0.0;
    for(i=0;i<sim->nBins;i++){
      totS+=smooSim[i];
      z=(float)sim->maxZ-(float)i*simIO->res;
      CofGs+=smooSim[i]*z;
    }

    thresh=0.0001*totS;
    CofGs/=totS;
    correl[2*k+1]=CofGl-CofGs;

    /*find sim bounds*/
    for(i=0;i<sim->nBins;i++){
      if(smooSim[i]>thresh){
        z=(float)sim->maxZ-(float)i*simIO->res;
        eSx=z;
        break;
      } 
    }
    for(i=sim->nBins-1;i>=0;i--){
      if(smooSim[i]>thresh){
        z=(float)sim->maxZ-(float)i*simIO->res;
        sSx=z;
        break;
      }
    }

    startX=(sSx<sLx)?sSx:sLx;
    endX=(eSx>eLx)?eSx:eLx;
    numb=(int)(fabs(endX-startX)/simIO->res);

    /*means*/
    meanS=meanL=0.0;
    for(i=0;i<lvis->nBins;i++)if((zShift[i]>=startX)&&(zShift[i]<=endX))meanL+=truth[i];
    for(i=0;i<sim->nBins;i++){
      z=(float)sim->maxZ-(float)i*simIO->res;
      if((z>=startX)&&(z<=endX))meanS+=smooSim[i];
    }
    meanL/=(float)numb;
    meanS/=(float)numb;
    /*stdev*/
    stdevS=stdevL=0.0;
    for(i=0;i<lvis->nBins;i++)if((zShift[i]>=startX)&&(zShift[i]<=endX))stdevL+=(truth[i]-meanL)*(truth[i]-meanL);
    for(i=0;i<sim->nBins;i++){
      z=(float)sim->maxZ-(float)i*simIO->res;
      if((z>=startX)&&(z<=endX))stdevS+=(smooSim[i]-meanS)*(smooSim[i]-meanS);
    }
    stdevL=sqrt(stdevL/(float)numb);
    stdevS=sqrt(stdevS/(float)numb);

    /*shared variance*/
    sumProd=0.0;
    for(i=lvis->nBins-1;i>=0;i--){
      if((zShift[i]<startX)||(zShift[i]>endX))continue;
      minSepSq=100000.0;
      for(j=0;j<sim->nBins;j++){
        z=(float)sim->maxZ-(float)j*simIO->res;
        if((z<startX)||(z>endX))continue;
        sepSq=pow((z-zShift[i]),2.0);
        if(sepSq<minSepSq){
          minSepSq=sepSq;
          bin=j;
        }
      }
      sumProd+=(truth[i]-meanL)*(smooSim[bin]-meanS);
    }
    correl[2*k]=(sumProd/(float)numb)/(stdevL*stdevS);
    TIDY(smooSim);
  }/*wave type loop*/

  /*tidy up*/
  TIDY(zShift);
  return(correl);
}/*waveCorrel*/


/*####################################################*/
/*denoiseall LVIS data*/

float **denoiseAllLvis(dataStruct **lvis,control *dimage)
{
  int i=0,j=0;
  float **denoised=NULL,tot=0;

  /*allocate space*/
  denoised=fFalloc(dimage->nLvis,"denoised waveforms",0);

  /*make sure fSigma is not 0*/
  dimage->simIO.linkFsig=dimage->simIO.fSigma;

  /*loop over LVIS footprints*/
  for(i=0;i<dimage->nLvis;i++){
    denoised[i]=processFloWave(lvis[i]->wave[0],lvis[i]->nBins,dimage->lvisIO.den,1.0);
    /*rescale*/
    tot=0.0;
    for(j=0;j<lvis[i]->nBins;j++)tot+=denoised[i][j];
    for(j=0;j<lvis[i]->nBins;j++)denoised[i][j]/=tot;
    /*determine beam sensitivity*/
    if(dimage->minSense>0.0)lvis[i]->beamSense=findBlairSense(lvis[i],&dimage->simIO);
    else                    lvis[i]->beamSense=0.0;
  }/*LVIS footprint loop*/

  return(denoised);
}/*denoiseAllLvis*/


/*####################################################*/
/*shift coordinate centres*/

double **shiftPrints(double **coords,double xOff,double yOff,int numb)
{
  int i=0;
  double **newCoords=NULL;

  /*allocate space*/
  newCoords=dDalloc(numb,"newCoords",0);

  /*loop overall points*/
  for(i=0;i<numb;i++){
    newCoords[i]=dalloc(2,"newCoords",i+1);
    newCoords[i][0]=coords[i][0]+xOff;
    newCoords[i][1]=coords[i][1]+yOff;
  }

  return(newCoords);
}/*shiftPrints*/


/*####################################################*/
/*read multiple ALS files and save relevant data*/

pCloudStruct **readMultiALS(control *dimage,dataStruct **lvis)
{
  int i=0;
  uint32_t totPoints=0;
  pCloudStruct **als=NULL;
  lasFile *las=NULL;
  void copyLvisCoords(gediRatStruct *,dataStruct **,int,int,int);


  /*determine bounds*/
  copyLvisCoords(&dimage->gediRat,lvis,dimage->nLvis,dimage->aEPSG,dimage->lEPSG);
  setGediGrid(&dimage->simIO,&dimage->gediRat);
 /*account for jittering*/
  dimage->gediRat.globMinX+=dimage->origin[0]-(double)dimage->maxShift;
  dimage->gediRat.globMaxX+=dimage->origin[0]+(double)dimage->maxShift;
  dimage->gediRat.globMinY+=dimage->origin[1]-(double)dimage->maxShift;
  dimage->gediRat.globMaxY+=dimage->origin[1]+(double)dimage->maxShift;
  dimage->gediRat.maxSep+=+(double)dimage->maxShift;

  /*allocate space*/
  if(!(als=(pCloudStruct **)calloc(dimage->simIO.nFiles,sizeof(pCloudStruct *)))){
    fprintf(stderr,"error in lvis data allocation.\n");
    exit(1);
  }

  /*loop over ALS files*/
  for(i=0;i<dimage->simIO.nFiles;i++){
    /*read header*/
    las=readLasHead(dimage->simIO.inList[i],dimage->pBuffSize);

    /*read data*/
    als[i]=readALSdata(las,&dimage->gediRat,i);
    if(als[i]->nPoints>0)fprintf(stdout,"Found %d ALS points in file %d of %d\n",als[i]->nPoints,i+1,dimage->simIO.nFiles);

    /*tidy up*/
    las=tidyLasFile(las);

    totPoints+=als[i]->nPoints;
  }/*ALS file loop*/

  /*is there data?*/
  if(totPoints==0){
    fprintf(stderr,"No ALS data contained\n");
    exit(1);
  }

  return(als);
}/*readMultiALS*/


/*####################################################*/
/*copy LVIS coords into ALS structure*/

void copyLvisCoords(gediRatStruct *gediRat,dataStruct **lvis,int nLvis,int aEPSG,int lEPSG)
{
  int i=0;
  double *x=NULL,*y=NULL,*z=NULL;
  OGRCoordinateTransformationH hTransform;
  OGRSpatialReferenceH hSourceSRS,hTargetSRS;
  OGRErr err;

  gediRat->gNx=nLvis;
  gediRat->gNy=1;
  gediRat->coords=dDalloc(gediRat->gNx,"coord list",0);
  gediRat->waveIDlist=NULL;


  if(aEPSG!=lEPSG){
    x=dalloc(nLvis,"x",0);
    y=dalloc(nLvis,"y",0);
    z=dalloc(nLvis,"z",0);
    for(i=0;i<nLvis;i++){
      x[i]=lvis[i]->lon;
      y[i]=lvis[i]->lat;
      z[i]=0.0;
    }

    hSourceSRS=OSRNewSpatialReference(NULL);
    hTargetSRS=OSRNewSpatialReference(NULL);
    err=OSRImportFromEPSG(hTargetSRS,aEPSG);
    err=OSRImportFromEPSG(hSourceSRS,lEPSG);
    hTransform=OCTNewCoordinateTransformation(hSourceSRS,hTargetSRS);
    OCTTransform(hTransform,nLvis,x,y,z);
    OCTDestroyCoordinateTransformation(hTransform);
    OSRDestroySpatialReference(hSourceSRS);
    OSRDestroySpatialReference(hTargetSRS);

    for(i=0;i<nLvis;i++){
      gediRat->coords[i]=dalloc(2,"coords",i+1);
      gediRat->coords[i][0]=x[i];
      gediRat->coords[i][1]=y[i];
    }
  }else{
    for(i=0;i<nLvis;i++){
      gediRat->coords[i]=dalloc(2,"coords",i+1);
      gediRat->coords[i][0]=lvis[i]->lon;
      gediRat->coords[i][1]=lvis[i]->lat;
    }
  }

  TIDY(x);
  TIDY(y);
  TIDY(z);
  return;
}/*copyLvisCoords*/


/*####################################################*/
/*read multiple LVIS files and save relevant data*/

dataStruct **readMultiLVIS(control *dimage,float *res)
{
  int i=0;
  dataStruct **lvis=NULL;
  dataStruct **copyLVIShdf(lvisHDF *,dataStruct **,control *,double *);
  dataStruct **copyLVISlgw(char *,dataStruct **,control *,double *);
  dataStruct **copyGEDIhdf(gediHDF *,dataStruct **,control *,double *);
  double *bounds=NULL,tempBound[4];
  lvisHDF *hdf=NULL;        /*LVIS HDF5 structure*/
  gediHDF *simHDF=NULL;     /*GEDI HDF5 structure*/

  dimage->nLvis=0;

  /*reproject bounds*/
  tempBound[0]=dimage->minX;
  tempBound[1]=dimage->minY;
  tempBound[2]=dimage->maxX;
  tempBound[3]=dimage->maxY;
  bounds=reprojectWaveBounds(&(tempBound[0]),dimage->aEPSG,dimage->lEPSG);

  /*loop over lvis files*/
  for(i=0;i<dimage->lvisIO.nFiles;i++){
    if(dimage->useLvisHDF){
      /*read HDF5*/
      hdf=readLVIShdf(dimage->lvisIO.inList[i]);
      /*unpack*/
      lvis=copyLVIShdf(hdf,lvis,dimage,bounds);
      /*tidy up*/
      hdf=tidyLVISstruct(hdf);
    }else if(dimage->useLvisLGW){
      /*unpack*/
      lvis=copyLVISlgw(dimage->lvisIO.inList[i],lvis,dimage,bounds);
    }else if(dimage->useGediHDF){
      /*read HDF*/
      simHDF=readGediHDF(dimage->lvisIO.inList[i],&dimage->lvisIO);
      /*unpack simulated HDF data*/
      lvis=copyGEDIhdf(simHDF,lvis,dimage,bounds);
      /*tidy up*/
      simHDF=tidyGediHDF(simHDF);
    }
  }/*file loop*/

  /*res for simulations*/
  *res=0.0;
  for(i=0;i<dimage->nLvis;i++)*res+=lvis[i]->res;
  *res/=(float)dimage->nLvis;


  fprintf(stdout,"Found %d LVIS\n",dimage->nLvis);
  if(dimage->nLvis==0){
    fprintf(stderr,"No large-footprints found\n");
    exit(1);
  }
  return(lvis);
}/*readMultiLVIS*/


/*####################################################*/
/*unpack relevant data from LGW*/

dataStruct **copyLVISlgw(char *namen,dataStruct **lvis,control *dimage,double *bounds)
{
  int i=0,nNew=0;
  int oldN=0;
  double x=0,y=0;
  lvisLGWstruct lgw;
  dataStruct *tempData=NULL;

  /*as this is overwritten later*/
  oldN=dimage->lvisIO.nFiles;

  /*count number of new files*/
  nNew=0;
  lgw.data=NULL;
  lgw.ipoo=NULL;
  tempData=readBinaryLVIS(namen,&lgw,0,&dimage->lvisIO);
  for(i=0;i<lgw.nWaves;i++){
    x=(lgw.data[i].lon0+lgw.data[i].lon431)/2.0;
    y=(lgw.data[i].lat0+lgw.data[i].lat431)/2.0;
    if((dimage->lEPSG==4326)&&(x>180.0))x-=360.0;

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3]))nNew++;
  }
  if(tempData){
    TTIDY((void **)tempData->wave,tempData->nWaveTypes);
    TTIDY((void **)tempData->ground,tempData->nWaveTypes);
    TIDY(tempData->noised);
    TIDY(tempData->totE);
    TIDY(tempData->z);
    TIDY(tempData);
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
    if((dimage->lEPSG==4326)&&(x>180.0))x-=360.0;

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3])){
      lvis[nNew+dimage->nLvis]=readBinaryLVIS(namen,&lgw,i,&dimage->lvisIO);
      nNew++;
    }
  }

  /*keep count*/
  dimage->nLvis+=nNew;
  dimage->lvisIO.nFiles=oldN;
  TIDY(lgw.data);
  if(lgw.ipoo){
    fclose(lgw.ipoo);
    lgw.ipoo=NULL;
  }
  return(lvis);
}/*copyLVISlgw*/


/*####################################################*/
/*copy data from simulated HDF5*/

dataStruct **copyGEDIhdf(gediHDF *hdf,dataStruct **lvis,control *dimage,double *bounds)
{
  int i=0,nNew=0;
  double x=0,y=0;


  /*count number within*/
  nNew=0;
  for(i=0;i<hdf->nWaves;i++){
    x=hdf->lon[i];
    y=hdf->lat[i];
    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3]))nNew++;
  }/*point loop*/

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
    x=hdf->lon[i];
    y=hdf->lat[i];

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3])){
      lvis[nNew+dimage->nLvis]=unpackHDFgedi(NULL,&dimage->lvisIO,&hdf,i);
      nNew++;
    }
  }

  /*keep count*/
  dimage->nLvis+=nNew;

  return(lvis);
}/*copyGEDIhdf*/


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
      lvis[nNew+dimage->nLvis]=unpackHDFlvis(NULL,&hdf,&dimage->lvisIO,i);
      nNew++;
    }
  }

  /*keep count*/
  dimage->nLvis+=nNew;

  return(lvis);
}/*copyLVIShdf*/


/*####################################################*/
/*set bounds from ALS data*/

void setALSbounds(control *dimage)
{
  int i=0;
  lasFile *las=NULL;
  double buff=0;

  /*set nonesense values*/
  dimage->minX=dimage->minY=10000000000.0;
  dimage->maxX=dimage->maxY=-10000000000.0;

  /*loop over ALS files and read bounds from header*/
  for(i=0;i<dimage->simIO.nFiles;i++){
    las=readLasHead(dimage->simIO.inList[i],dimage->pBuffSize);

    if(las->minB[0]<dimage->minX)dimage->minX=dimage->lvisIO.bounds[0]=las->minB[0];
    if(las->minB[1]<dimage->minY)dimage->minY=dimage->lvisIO.bounds[1]=las->minB[1];
    if(las->maxB[0]>dimage->maxX)dimage->maxX=dimage->lvisIO.bounds[2]=las->maxB[0];
    if(las->maxB[1]>dimage->maxY)dimage->maxY=dimage->lvisIO.bounds[3]=las->maxB[1];

    las=tidyLasFile(las);
  }/*ALS file loop*/

  /*add on a buffer for uncertainty*/
  buff=(double)dimage->maxShift;
  dimage->minX-=buff;
  dimage->minY-=buff;
  dimage->maxX+=buff;
  dimage->maxY+=buff;
  dimage->lvisIO.bounds[0]-=buff;
  dimage->lvisIO.bounds[1]-=buff;
  dimage->lvisIO.bounds[2]+=buff;
  dimage->lvisIO.bounds[3]+=buff;

  dimage->useBounds=1;
  return;
}/*setALSbounds*/


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
  strcpy(dimage->outNamen,"teast.correl");
  dimage->useLvisHDF=1;
  dimage->useLvisLGW=0;
  dimage->useGediHDF=0;
  dimage->aEPSG=dimage->lvisIO.bEPSG=32632;
  dimage->lEPSG=dimage->lvisIO.wEPSG=4326;
  dimage->pBuffSize=(uint64_t)200000000;
  dimage->filtOutli=1;    /*filter outliers*/
  dimage->maxZen=100000.0;
  dimage->minDense=0.0;
  dimage->minSense=0.0;
  dimage->leaveEmpty=0;
  dimage->opoo=NULL;
  dimage->useBounds=0;     /*read bounds from ALS file*/
  dimage->writeFinWave=0;  /*don't write out the final waveforms*/

  /*octree*/
  dimage->gediRat.useOctree=1;
  dimage->gediRat.octLevels=0;
  dimage->gediRat.nOctTop=30;
  dimage->gediRat.octree=NULL;

  /*LVIS params for sim*/
  dimage->simIO.fSigma=5.5;    /*default value for GEDI. 4.31 if LVIS*/;
  dimage->simIO.pSigma=0.6893;
  dimage->gediRat.iThresh=0.002;
  dimage->simIO.useCount=1;
  dimage->simIO.useInt=0;
  dimage->simIO.useFrac=0;

  /*waveform reading settings*/
  dimage->lvisIO.useInt=dimage->lvisIO.useFrac=0;
  dimage->lvisIO.useCount=1;
  dimage->lvisIO.ground=0;
  dimage->lvisIO.nMessages=200;

  /*LVIS denoising*/
  setDenoiseDefault(dimage->lvisIO.den);
  dimage->lvisIO.den->varNoise=1;
  dimage->lvisIO.den->statsLen=12.0;
  dimage->lvisIO.den->noiseTrack=1;
  dimage->lvisIO.den->threshScale=4.0;

  /*bullseyte params*/
  dimage->maxShift=10.0;     /*maximum distance to shift*/
  dimage->shiftStep=1.0;     /*distance to shift steps*/
  dimage->maxVshift=-10.0;   /*maximum vertical distance to shift. Ignored if negative*/
  dimage->vShiftStep=0.2;    /*distance to shift steps vertically*/
  for(i=0;i<3;i++)dimage->origin[i]=0.0;    /*origin to shift around*/

  /*bulseye optimsation*/
  dimage->fullBull=1;        /*do full bullseye plot*/
  dimage->largeErr=0;        /*don't do a large error search*/
  dimage->findFsig=1;
  dimage->maxIter=100;
  dimage->optTol=0.01;
  dimage->writeSimProg=0;

  /*all data*/
  dimage->minX=-100000000.0;
  dimage->minY=-100000000.0;
  dimage->maxX=100000000.0;
  dimage->maxY=100000000.0;

  /*simulation settings*/
  dimage->gediRat.readALSonce=1;    /*read all ALS data once*/
  dimage->gediRat.readWave=0;       /*do not read waveform switch*/
  dimage->gediRat.useShadow=0;      /*do not account for shadowing through voxelisation*/
  dimage->gediRat.maxScanAng=1000000000.0;    /*maximum scan angle*/
  dimage->gediRat.coords=NULL;     /*list of coordinates*/
  dimage->gediRat.readWave=0;
  dimage->gediRat.pulseAfter=1;  /*smooth after for speed*/
  dimage->gediRat.sideLobe=0;   /*no side lobes*/
  dimage->gediRat.lobeAng=0.0;
  dimage->gediRat.checkCover=1;
  dimage->gediRat.normCover=0;
  dimage->gediRat.cleanOut=0;
  dimage->gediRat.topHat=0;
  dimage->gediRat.useShadow=0;
  dimage->gediRat.maxScanAng=1000000.0;   /*maximum scan angle*/
  dimage->gediRat.coords=NULL;       /*list of coordinates*/
  dimage->gediRat.waveIDlist=NULL;   /*list of waveform IDs*/
  dimage->gediRat.doDecon=0;
  dimage->gediRat.indDecon=0;
  dimage->gediRat.decimate=1.0;  /*accept all ALS beams*/
  dimage->simIO.nMessages=200;
  dimage->simIO.pRes=0.01;
  dimage->simIO.ground=0;
  dimage->simIO.readPulse=0;
  setDenoiseDefault(dimage->simIO.den);
  dimage->simIO.den->meanN=0.0;
  dimage->simIO.den->thresh=0.0;
  dimage->simIO.den->minWidth=0.0;


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
      }else if(!strncasecmp(argv[i],"-listLvis",9)){
        checkArguments(1,i,argc,"-listLvis");
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
      }else if(!strncasecmp(argv[i],"-listAls",8)){
        checkArguments(1,i,argc,"-listAls");
        TTIDY((void **)dimage->simIO.inList,dimage->simIO.nFiles);
        dimage->simIO.inList=readInList(&dimage->simIO.nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(4,i,argc,"-bounds");
        dimage->useBounds=1;
        dimage->minX=dimage->lvisIO.bounds[0]=atof(argv[++i]);
        dimage->minY=dimage->lvisIO.bounds[1]=atof(argv[++i]);
        dimage->maxX=dimage->lvisIO.bounds[2]=atof(argv[++i]);
        dimage->maxY=dimage->lvisIO.bounds[3]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-lgw",4)){
        dimage->useLvisHDF=0;
        dimage->useLvisLGW=1;
        dimage->useGediHDF=0;
      }else if(!strncasecmp(argv[i],"-readHDFgedi",11)){
        dimage->useLvisHDF=0;
        dimage->useLvisLGW=0;
        dimage->useGediHDF=1;
      }else if(!strncasecmp(argv[i],"-aEPSG",6)){
        checkArguments(1,i,argc,"-aEPSG");
        dimage->aEPSG=dimage->lvisIO.bEPSG=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-lEPSG",6)){
        checkArguments(1,i,argc,"-lEPSG");
        dimage->lEPSG=dimage->lvisIO.wEPSG=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pBuff",6)){
        checkArguments(1,i,argc,"-pBuff");
        dimage->pBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0);
      }else if(!strncasecmp(argv[i],"-readPulse",10)){
        checkArguments(1,i,argc,"-readPulse");
        dimage->simIO.readPulse=1;
        strcpy(dimage->simIO.pulseFile,argv[++i]);
      }else if(!strncasecmp(argv[i],"-pSigma",7)){
        checkArguments(1,i,argc,"-pSigma");
        dimage->simIO.pSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fSigma",7)){
        checkArguments(1,i,argc,"-fSigma");
        dimage->simIO.fSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-smooth",7)){
        checkArguments(1,i,argc,"-smooth");
        dimage->lvisIO.den->sWidth=dimage->simIO.den->sWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxShift",9)){
        checkArguments(1,i,argc,"-maxShift");
        dimage->maxShift=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxVshift",10)){
        checkArguments(1,i,argc,"-maxVshift");
        dimage->maxVshift=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-step",5)){
        checkArguments(1,i,argc,"-step");
        dimage->shiftStep=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-vStep",6)){
        checkArguments(1,i,argc,"-vStep");
        dimage->vShiftStep=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-offset",8)){
        checkArguments(1,i,argc,"-offset");
        dimage->origin[2]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-hOffset",9)){
        checkArguments(2,i,argc,"-hOffset");
        dimage->origin[0]=atof(argv[++i]);
        dimage->origin[1]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noNorm",7)){
        dimage->gediRat.normCover=0;
      }else if(!strncasecmp(argv[i],"-norm",5)){
        dimage->gediRat.normCover=1;
      }else if(!strncasecmp(argv[i],"-decimate",9)){
        checkArguments(1,i,argc,"-decimate");
        dimage->gediRat.decimate=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noFilt",7)){
        dimage->filtOutli=0;
      }else if(!strncasecmp(argv[i],"-noOctree",9)){
        dimage->gediRat.useOctree=0;
      }else if(!strncasecmp(argv[i],"-octLevels",9)){
        checkArguments(1,i,argc,"-octLevels");
        dimage->gediRat.octLevels=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-nOctPix",8)){
        checkArguments(1,i,argc,"-nOctPix");
        dimage->gediRat.nOctTop=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxZen",7)){
        checkArguments(1,i,argc,"-maxZen");
        dimage->maxZen=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pulseBefore",12)){
        dimage->gediRat.pulseAfter=0;  /*smooth befpre to prevent aliasing*/
      }else if(!strncasecmp(argv[i],"-allSimMeth",11)){
        dimage->simIO.useCount=dimage->simIO.useInt=dimage->simIO.useFrac=1;
      }else if(!strncasecmp(argv[i],"-simplex",8)){
        dimage->fullBull=0;
      }else if(!strncasecmp(argv[i],"-maxIter",8)){
        checkArguments(1,i,argc,"-maxIter");
        dimage->maxIter=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-optTol",7)){
        checkArguments(1,i,argc,"-optTol");
        dimage->optTol=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minDense",9)){
        checkArguments(1,i,argc,"-minDense");
        dimage->minDense=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minSense",9)){
        checkArguments(1,i,argc,"-minSense");
        dimage->minSense=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-quickGeo",9)){
        dimage->largeErr=1;
        dimage->fullBull=0;
        dimage->maxShift=100.0;    /*expected geolocation error, plus a bit*/
        dimage->shiftStep=20.0;    /*expected width of correlation hole*/
      }else if(!strncasecmp(argv[i],"-geoError",9)){
        checkArguments(2,i,argc,"-geoError");
        dimage->largeErr=1;  
        dimage->fullBull=0;
        dimage->maxShift=atof(argv[++i]);     /*expected geolocation error, plus a bit*/
        dimage->shiftStep=atof(argv[++i]);    /*expected width of correlation hole*/
      }else if(!strncasecmp(argv[i],"-leaveEmpty",11)){
        dimage->leaveEmpty=1;
      }else if(!strncasecmp(argv[i],"-writeSimProg",13)){
        dimage->writeSimProg=1;
      }else if(!strncasecmp(argv[i],"-writeWaves",11)){
        checkArguments(1,i,argc,"-writeWaves");
        strcpy(dimage->waveNamen,argv[++i]);
        dimage->writeFinWave=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to colocate large-footprint and small-footprint lidar data\n#####\n\n-output name;     output filename\n-listAls list;    input file list for multiple als files\n-als file;        input als file\n-lvis file;       single input LVIS file\n-listLvis file;   list of multiple LVIS files\n-lgw;             LVIS is in lgw (default is LVIS hdf5)\n-readHDFgedi;     read GEDI HDF5 input (default is LVIS hdf5)\n-lEPSG epsg;      LVIS projection\n-aEPSG epsg;      ALS projection\n-pSigma x;        pulse length, sigma in metres\n-fSigma x;        footprint width, sigma in metres\n-readPulse file;  pulse shape\n-pulseBefore;     apply pulse shape before binning to prevent aliasing\n-minDense x;      minimum ALS beam density to accept\n-minSense x;      minimum waveform beam sensitivity to accept\n-smooth sig;      smooth both waves before comparing\n-maxShift x;      horizontal distance to search over\n-step x;          horizontal step size\n-maxVshift x;     vertical distance to search over\n-vStep z;         vertical step size\n-hOffset dx dy;   centre of horizontal offsets\n-offset z;        vertical datum offset\n-bounds minX minY maxX maxY;    bounds to use, in ALS projection\n-noNorm;          don't correct sims for ALS densiy variations\n-norm;            correct sims for ALS densiy variations\n-decimate f;      decimate ALS point cloud by a factor, to save RAM\n-noFilt;          don't filter outliers from correlation\n-allSimMeth;      use all simulation methods\n\n# Optimisation\n-simplex;         use simplex optimisation rather than doing the full bullseye plot\n-maxIter n;       maximum number of iterations\n-optTol x;        tolerance for optimisation\n-quickGeo;        perform a rapid geolocation of initial GEDI data, using default error values\n-geoError expError correlDist;       perform a rapid geolocation of initial GEDI data, providing an expected geolocation error and an expected correlation distance\n-writeSimProg;    write progress of simplex to output\n-writeWaves name;  write out final waveforms as HDF5 when using simplex\n\n# Octree\n-noOctree;       do not use an octree\n-octLevels n;    number of octree levels to use\n-nOctPix n;      number of octree pixels along a side for the top level\n-maxZen zen;     maximum zenith angle to use, degrees\n-leaveEmpty;     exit if there are no usable footprints\n\n");
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

