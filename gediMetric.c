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
#include "gediNoise.h"

#define USEPHOTON

#ifdef USEPHOTON
#include "photonCount.h"
#endif



/*##############################*/
/*# Generates metrics from     #*/
/*# simulated GEDI waveforms   #*/
/*# 2015 svenhancock@gmail.com #*/
/*##############################*/

/*#######################################*/
/*# Copyright 2015-2016, Steven Hancock #*/
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



/*tolerances*/
#define TOL 0.00001
#define MINERR 0.0000001

/*element reflectance*/
float rhoG;
float rhoC;

/*###########################################################*/
/*function definition used here only*/

float *findLAIprofile(float *,float,int,float,int *,double,float,double *,float);



/*###########################################################*/
/*LVIS level2 data*/

typedef struct{
  uint64_t numb;      /*number of records*/
  uint32_t *lfid;     /*LVIS file identifier*/
  uint32_t *shotN;    /*LVIS shotnumber*/
  float *zG;          /*ground elevation*/
}lvisL2struct;


/*####################################*/
/*emtpy structure if photon counting not provided*/
#ifndef USEPHOTON
typedef struct{
  void *nothing;
}photonStruct;
#endif


/*####################################*/
/*control structure*/

typedef struct{
  /*input/output*/
  gediIOstruct gediIO; /*input/output structure*/
  char outRoot[200];
  FILE *opooGauss;  /*Gaussian parameter output*/
  FILE *opooMet;    /*waveform metric output*/
  int maxGauss;     /*maximum number of Gaussians for output*/

  /*level2 LVIS for ZG*/
  char l2namen[200]; /*list of level2 filenames*/
  char readL2;      /*switch to read L2 or not*/

  /*switches*/
  char writeFit;    /*write fitted wave switch*/
  float rhRes;      /*rh resolution*/
  char bayesGround; /*Bayseian ground finding*/
  char noRHgauss;   /*do not do Gaussian fitting*/
  char renoiseWave; /*remove noise before adding*/
  char readBinLVIS;  /*read binary LVIS rather than a list of ASCII files*/
  char readHDFlvis;  /*read HDF5 LVIS rather than ASCII*/
  char readHDFgedi;  /*read HDF5 GEDI rather than ASCII*/
  char coord2dp;     /*round up coords to 2dp when writing*/
  char useBounds;    /*when we will process only a subset of bounds*/
  char writeGauss;   /*write Gaussian parameters*/
  float laiRes;      /*LAI profile resolution*/
  float maxLAIh;     /*maximum height bin of LAI profile. Put all above this in top bin*/

  /*noise parameters*/
  noisePar noise;  /*noise adding structure*/
  float bThresh;   /*bounds threshold*/

  /*LVIS or HDF data*/
  lvisLGWstruct lvis;   /*LVIS lgw structure*/
  lvisHDF *hdfLvis;     /*LVIS HDF5 structure*/
  lvisL2struct *lvisL2; /*LVIS level2 data*/
  gediHDF *hdfGedi;     /*GEDI HDF5 structure*/

  /*bounds for subsets*/
  double minX;
  double maxX;
  double minY;
  double maxY;

  /*photon counting*/
  char ice2;         /*ICESat-2 mode. GEDI by default*/
  photonStruct photonCount;  /*photon counting structure*/

  /*others*/
  float rhoRatio; /*ration of canopy to ground reflectance*/
  float gTol;     /*toleranve used to label ALS ground finding*/
  float zen;      /*zenith angle*/
  float fhdHistRes;/*resolution for FHD histogram method*/
}control;


/*###########################################################*/
/*Bayseian groiund structure*/

typedef struct{
  double gHeight;   /*ground elevation*/
  float cov;        /*canopy cover*/
  float slope;      /*slope, degrees*/
}bGround;


/*###########################################################*/
/*metric structure*/

typedef struct{
  float *rh;        /*rh metrics using Gaussian ground*/
  float *rhMax;     /*rh metrics using max ground*/
  float *rhInfl;    /*rh metrics using inflection ground*/
  float *rhReal;    /*rh metric from real ground*/
  int nRH;          /*number of RH metrics*/
  float FHD;        /*foliage height diversity, all waveform, wave*/
  float FHDhist;    /*foliage height diversity, all waveform, hist*/
  float FHDcan;     /*foliage height diversity, canopy, wave*/
  float FHDcanH;    /*foliage height diversity, canopy, hist*/
  float FHDcanGauss;/*foliage height diversity, canopy from Gaussian fitting, wave*/
  float FHDcanGhist;/*foliage height diversity, canopy from Gaussian fitting, hist*/
  int nLm;          /*number of L-moments*/
  //float *LmomGau;   /*L-moments from Gaussian fit*/
  //float *LmomRea;   /*L-moments from ALS ground*/
  //float *LmomInf;   /*L-moments from inflection point*/
  //float *LmomMax;   /*L-moments from maximum*/
  float cov;        /*canopy cover for gaussian fitting*/
  double gHeight;   /*ground height from Gaussians*/
  double maxGround; /*ground height from maximum*/
  double inflGround;/*ground height from inflection*/
  double tElev;     /*top elevation*/
  double bElev;     /*bottom elevation*/
  float leExt;      /*Lefsky's leading edge extent*/
  float teExt;      /*Lefsky's trailing edge extent*/
  float covHalfG;   /*cover from Bryan's half, Gaussian*/
  float covHalfI;   /*cover from Bryan's half, Inflection*/
  float covHalfM;   /*cover from Bryan's half, maximum*/
  float covHalfB;   /*cover from Bryan's half, Bayesian*/
  float totE;       /*total energy after denoising*/
  float blairSense; /*Blair sensitivity metric*/
  float niM2;       /*Ni metric with c=2*/
  float niM21;      /*Ni metric with c=2.1*/
  float *tLAI;      /*true LAI profile*/
  float *gLAI;      /*LAI profile with Gaussian ground removal*/
  float *hgLAI;     /*LAI profile with halp width ground removal, Gaussian elevation*/
  float *hiLAI;     /*LAI profile with halp width ground removal, inflection elevation*/
  float *hmLAI;     /*LAI profile with halp width ground removal, maximum elevation*/
  int laiBins;      /*number of LAI bins*/

  int nBgr;         /*number of ground estimates*/
  bGround *bGr;     /*Bayesian ground structure*/
  double bayGround; /*Bayesian ground elevation*/
}metStruct;


/*###########################################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;
int j=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct *data=NULL;
  metStruct *metric=NULL;
  void setL2ground(dataStruct *,int,control *);
  void findMetrics(metStruct *,float *,int,float *,float *,int,double *,control *,dataStruct *);
  void tidySMoothPulse();
  void alignElevation(double,double,float *,int);
  void writeResults(dataStruct *,control *,metStruct *,int,float *,float *,char *);
  void determineTruth(dataStruct *,control *);
  void modifyTruth(dataStruct *,noisePar *);
  void checkWaveformBounds(dataStruct *,control *);
  void photonCountCloud(float *,dataStruct *,photonStruct *,char *,int);
  float *processed=NULL,*denoised=NULL;;

  /*read command Line*/
  dimage=readCommands(argc,argv);

  /*set link noise if needed*/
  dimage->noise.linkSig=setNoiseSigma(dimage->noise.linkM,dimage->noise.linkCov,dimage->gediIO.linkPsig,dimage->gediIO.linkFsig,rhoC,rhoG);

  /*allocate metric array*/
  if(!(metric=(metStruct *)calloc(1,sizeof(metStruct)))){
    fprintf(stderr,"error metric structure allocation.\n");
    exit(1);
  }

  /*loop over files*/
  for(i=0;i<dimage->gediIO.nFiles;i++){
    if((i%dimage->gediIO.nMessages)==0)fprintf(stdout,"Wave %d of %d\n",i+1,dimage->gediIO.nFiles);

    /*read waveform*/
    if(dimage->readBinLVIS)     data=readBinaryLVIS(dimage->gediIO.inList[0],&dimage->lvis,i,&dimage->gediIO);
    else if(dimage->readHDFlvis)data=unpackHDFlvis(dimage->gediIO.inList[0],&dimage->hdfLvis,&dimage->gediIO,i);
    else if(dimage->readHDFgedi)data=unpackHDFgedi(dimage->gediIO.inList[0],&dimage->gediIO,&dimage->hdfGedi,i);
    else                        data=readASCIIdata(dimage->gediIO.inList[i],&(dimage->gediIO));
    if(dimage->readL2)setL2ground(data,i,dimage);

    /*check bounds if needed*/
    if(dimage->useBounds)checkWaveformBounds(data,dimage);

    /*is the data usable*/
    if(data->usable){
      /*denoise and change pulse if needed*/
      if(dimage->renoiseWave)modifyTruth(data,&dimage->noise);

      /*determine truths before noising*/
      determineTruth(data,dimage);

      /*add noise if needed*/
      addNoise(data,&dimage->noise,dimage->gediIO.fSigma,dimage->gediIO.pSigma,dimage->gediIO.res,rhoC,rhoG);

      /*process waveform*/
      /*denoise*/
      denoised=processFloWave(data->noised,data->nBins,dimage->gediIO.den,1.0);

      /*are we in GEDI mode?*/
      if(!dimage->ice2){
        /*Gaussian fit*/
        if(dimage->noRHgauss==0)processed=processFloWave(denoised,data->nBins,dimage->gediIO.gFit,1.0);

        /*shift Gaussian centres to align to absolute elevation*/
        alignElevation(data->z[0],data->z[data->nBins-1],dimage->gediIO.gFit->gPar,dimage->gediIO.gFit->nGauss);

        /*determine metrics*/
        findMetrics(metric,dimage->gediIO.gFit->gPar,dimage->gediIO.gFit->nGauss,denoised,data->noised,data->nBins,data->z,dimage,data);

        /*write results*/
        if(dimage->readBinLVIS||dimage->readHDFlvis||dimage->readHDFgedi)writeResults(data,dimage,metric,i,denoised,processed,dimage->gediIO.inList[0]);
        else                                                             writeResults(data,dimage,metric,i,denoised,processed,dimage->gediIO.inList[i]);
      }else{  /*ICESat-2 mode*/
        photonCountCloud(denoised,data,&dimage->photonCount,dimage->outRoot,i);
      }/*operation mode switch*/
    }/*is the data usable*/


    /*tidy as we go along*/
    TIDY(processed);
    TIDY(denoised);
    if(data){
      TIDY(data->noised);
      if(dimage->readHDFgedi){  /*pointer to array. do not free*/
        data->wave[0]=NULL;
        if(data->ground)data->ground[0]=NULL;
      }
      TTIDY((void **)data->ground,data->nWaveTypes);
      TTIDY((void **)data->wave,data->nWaveTypes);
      TIDY(data->totE);
      TIDY(data->z);
      TIDY(data);
    }
    TIDY(dimage->gediIO.gFit->gPar);
    TIDY(dimage->gediIO.den->gPar);
    dimage->gediIO.den->nGauss=0;
    dimage->gediIO.gFit->nGauss=0;
    TIDY(metric->rhMax);
    TIDY(metric->rhInfl);
    TIDY(metric->rhReal);
    TIDY(metric->rh);
    TIDY(metric->bGr);
    TIDY(metric->tLAI);
    TIDY(metric->gLAI);
    TIDY(metric->hgLAI);
    TIDY(metric->hiLAI);
    TIDY(metric->hmLAI);
    //TIDY(metric->LmomGau);
    //TIDY(metric->LmomRea);
    //TIDY(metric->LmomInf);
    //TIDY(metric->LmomMax);
  }/*file loop*/

  /*TIDY LVIS data if it was read*/
  if(dimage->readBinLVIS)TIDY(dimage->lvis.data);
  if(dimage->readHDFgedi)dimage->hdfGedi=tidyGediHDF(dimage->hdfGedi);


  if(dimage->writeGauss)fprintf(stdout,"Written to %s.gauss.txt\n",dimage->outRoot);
  if(!dimage->ice2)fprintf(stdout,"Written to %s.metric.txt\n",dimage->outRoot);
  #ifdef USEPHOTON
  else             fprintf(stdout,"Written to %s\n",dimage->photonCount.outNamen);
  #endif


  /*tidy up arrays*/
  tidySMoothPulse();
  TIDY(metric);
  if(dimage){
    if(dimage->lvisL2){
      TIDY(dimage->lvisL2->lfid);
      TIDY(dimage->lvisL2->shotN);
      TIDY(dimage->lvisL2->zG);
      TIDY(dimage->lvisL2);
    }
    if(dimage->readBinLVIS||dimage->readHDFlvis||dimage->readHDFgedi)TTIDY((void **)dimage->gediIO.inList,1);
    else                                        TTIDY((void **)dimage->gediIO.inList,dimage->gediIO.nFiles);
    dimage->gediIO.inList=NULL;
    if(dimage->opooMet){
      fclose(dimage->opooMet);
      dimage->opooMet=NULL;
    }
    if(dimage->opooGauss){
      fclose(dimage->opooGauss);
      dimage->opooGauss=NULL;
    }
    #ifdef USEPHOTON
    if(dimage->photonCount.opoo){
      fclose(dimage->photonCount.opoo);
      dimage->photonCount.opoo=NULL;
    }
    TIDY(dimage->photonCount.prob);
    #endif
    if(dimage->gediIO.den){
      TTIDY((void **)dimage->gediIO.den->pulse,2);
      TIDY(dimage->gediIO.den->matchPulse);
      TIDY(dimage->gediIO.den->hardPulse);
      TIDY(dimage->gediIO.den);
    }
    if(dimage->gediIO.gFit){
      TTIDY((void **)dimage->gediIO.gFit->pulse,2);
      TIDY(dimage->gediIO.gFit);
    }
    dimage->hdfLvis=tidyLVISstruct(dimage->hdfLvis);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################################*/
/*check waveform bounds*/

void checkWaveformBounds(dataStruct *data,control *dimage)
{
  if((data->lon<dimage->minX)||(data->lat<dimage->minY)||(data->lon>dimage->maxX)||(data->lat>dimage->maxY))data->usable=0;

  return;
}/*checkWaveformBounds*/


/*####################################################*/
/*determine true variables*/

void determineTruth(dataStruct *data,control *dimage)
{
  int i=0;
  float totE=0,cumul=0;
  float groundOverlap(float *,float *,int);
  float groundMinAmp(float *,float *,int);
  float groundInflection(float *,float *,int);

  /*determine ground*/
  if(!data->demGround){  /*unless it's already been calculalated from the DEM*/
    if(dimage->gediIO.ground){
      totE=0.0;
      data->gElev=0.0;
      for(i=0;i<data->nBins;i++){
        totE+=data->ground[data->useType][i];
        data->gElev+=(double)data->ground[data->useType][i]*data->z[i];
      }
      if(totE>0.0)data->gElev/=(double)totE;
      else        data->gElev=-1000000.0;

      /*standard deviation as a measure of slope*/
      data->gStdev=0.0;
      for(i=0;i<data->nBins;i++){
        data->gStdev+=(float)((data->z[i]-data->gElev)*(data->z[i]-data->gElev))*data->ground[data->useType][i];
      }
      if(totE>0.0){
        data->gStdev=sqrt(data->gStdev/totE);
        if(data->gStdev>(data->pSigma+dimage->gTol)){
          data->slope=atan2(sqrt(data->gStdev*data->gStdev-(data->pSigma+dimage->gTol)*(data->pSigma+dimage->gTol)),data->fSigma)*180.0/M_PI;
        }else{
          data->slope=0.0;
        }
      }else        data->gStdev=data->slope=-1000000.0;
    }else{/*ground finding*/
      data->gElev=data->gStdev=data->slope=-1000000.0;
    }/*no ground finding*/
  }/*is the ground already defined from the DEM*/

  /*canopy top*/
  totE=0.0;
  for(i=0;i<data->nBins;i++)totE+=data->wave[data->useType][i];
  cumul=0.0;
  for(i=0;i<data->nBins;i++){
    cumul+=data->wave[data->useType][i]/totE;
    if(cumul>=dimage->bThresh){
      data->tElev=data->z[i];
      break;
    }
  }/*top finding*/

  /*canopy cover*/
  data->cov=waveformTrueCover(data,&dimage->gediIO,dimage->rhoRatio);

  /*understorey metrics*/
  if(dimage->gediIO.ground){
    data->gLap=groundOverlap(data->wave[data->useType],data->ground[data->useType],data->nBins);
    data->gMinimum=groundMinAmp(data->wave[data->useType],data->ground[data->useType],data->nBins);
    data->gInfl=groundInflection(data->wave[data->useType],data->ground[data->useType],data->nBins);
  }

  return;
}/*determineTruth*/


/*####################################################*/
/*ground inflection point amplitude*/

float groundInflection(float *wave,float *ground,int nBins)
{
  int i=0,maxInd=0;
  float gInfl=0,maxG=0;
  float *d2x=NULL;

  d2x=falloc(nBins,"d2x",0);

  /*find maximum of ground return*/
  maxG=0.0;
  for(i=nBins-1;i>=0;i--){
    if(ground[i]>maxG){
      maxG=ground[i];
      maxInd=i;
    }
  }

  /*determine derivatives*/
  for(i=1;i<nBins-1;i++){
    d2x[i]=2.0*wave[i]-(wave[i+1]+wave[i-1]);
  }

  /*find first d2x crossing point after max ground*/
  gInfl=-1.0;
  for(i=maxInd;i<nBins-1;i++){
    if(((d2x[i]<0.0)&&(d2x[i-1]>=0.0))||((d2x[i]>0.0)&&(d2x[i-1]<=0.0))){
      gInfl=d2x[i]-d2x[i-1];
      break;
    }
  }

  TIDY(d2x);
  return(gInfl);
}/*groundInflection*/


/*####################################################*/
/*amplitude of minimum between ground and canopy*/

float groundMinAmp(float *wave,float *ground,int nBins)
{
  int i=0,maxInd=0;
  float gMinimum=0;
  float max=0,min=0;

  /*determine max*/
  max=-1000.0;
  for(i=nBins-1;i>=0;i--){
    if(ground[i]>0.0){
      if(wave[i]>max){
        max=wave[i];
        maxInd=i;
      }
    }
  }

  /*find minimum after max*/
  min=max;
  for(i=maxInd;i>=0;i--){
    if(ground[i]>0.0){
      if(wave[i]<min){
        min=wave[i];
      }
    }
  }

  gMinimum=max-min;
  return(gMinimum);
}/*groundMinAmp*/


/*####################################################*/
/*determine ground canopy overlap*/

float groundOverlap(float *wave,float *ground,int nBins)
{
  int i=0;
  float canopy=0;
  float gLap=0,gTot=0;

  gLap=gTot=0.0;
  for(i=0;i<nBins;i++){
    canopy=wave[i]-ground[i];
    gLap+=(canopy>ground[i])?ground[i]:canopy;
    gTot+=ground[i];
  }

  if(gTot>0.0)gLap/=gTot;
  return(gLap);
}/*groundOverlap*/


/*####################################################*/
/*write results*/

void writeResults(dataStruct *data,control *dimage,metStruct *metric,int numb,float *denoised,float *processed,char *inNamen)
{
  int i=0,j=0;
  char waveNamen[500];
  char namen[420];
  float gauss(float,float,float);
  FILE *opoo=NULL;


  /*open file if needed*/
  if((dimage->opooGauss==NULL)&&(dimage->writeGauss)){
    sprintf(namen,"%s.gauss.txt",dimage->outRoot);
    if((dimage->opooGauss=fopen(namen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",namen);
      exit(1);
    }
    fprintf(dimage->opooGauss,"# 1 wave ID, 2 nGauss");
    for(i=0;i<dimage->maxGauss;i++)fprintf(dimage->opooGauss,", %d gauss %d mu, %d A, %d sig",3*i+3,i+1,3*i+4,3*i+5);
    fprintf(dimage->opooGauss,", %d wave name\n",3*i+6);
  }
  if(dimage->opooMet==NULL){
    sprintf(namen,"%s.metric.txt",dimage->outRoot);
    if((dimage->opooMet=fopen(namen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",namen);
      exit(1);
    }
    fprintf(dimage->opooMet,"# 1 wave ID, 2 true ground, 3 true top, 4 ground slope, 5 ALS cover, 6 gHeight, 7 maxGround, 8 inflGround, 9 signal top, 10 signal bottom, 11 cover, 12 leading edge ext, 13 trailing edge extent");
    for(i=0;i<metric->nRH;i++)fprintf(dimage->opooMet,", %d rhGauss %g",14+i,(float)i*dimage->rhRes);
    for(i=0;i<metric->nRH;i++)fprintf(dimage->opooMet,", %d rhMax %g",14+i+metric->nRH,(float)i*dimage->rhRes);
    for(i=0;i<metric->nRH;i++)fprintf(dimage->opooMet,", %d rhInfl %g",14+i+2*metric->nRH,(float)i*dimage->rhRes);
    for(i=0;i<metric->nRH;i++)fprintf(dimage->opooMet,", %d rhReal %g",14+i+3*metric->nRH,(float)i*dimage->rhRes);
    fprintf(dimage->opooMet,", %d filename",14+4*metric->nRH);
    if(dimage->bayesGround)fprintf(dimage->opooMet,", %d bayesGround",14+4*metric->nRH+1);
    fprintf(dimage->opooMet,", %d gaussHalfCov, %d maxHalfCov, %d infHalfCov, %d bayHalfCov",14+4*metric->nRH+1+dimage->bayesGround,14+4*metric->nRH+1+dimage->bayesGround+1,14+4*metric->nRH+1+dimage->bayesGround+2,14+4*metric->nRH+1+dimage->bayesGround+3);
    fprintf(dimage->opooMet,", %d pSigma, %d fSigma",14+4*metric->nRH+1+dimage->bayesGround+4,14+4*metric->nRH+1+dimage->bayesGround+5);
    fprintf(dimage->opooMet,", %d linkM, %d linkCov",14+4*metric->nRH+1+dimage->bayesGround+6,14+4*metric->nRH+1+dimage->bayesGround+7);
    fprintf(dimage->opooMet,", %d lon, %d lat",14+4*metric->nRH+1+dimage->bayesGround+8,14+4*metric->nRH+1+dimage->bayesGround+9);
    fprintf(dimage->opooMet,", %d groundOverlap, %d groundMin, %d groundInfl",14+4*metric->nRH+1+dimage->bayesGround+10,\
                             14+4*metric->nRH+1+dimage->bayesGround+11,14+4*metric->nRH+1+dimage->bayesGround+12);
    fprintf(dimage->opooMet,", %d waveEnergy, %d blairSense",14+4*metric->nRH+1+dimage->bayesGround+13,14+4*metric->nRH+1+dimage->bayesGround+14);
    fprintf(dimage->opooMet,", %d pointDense, %d beamDense",14+4*metric->nRH+1+dimage->bayesGround+15,14+4*metric->nRH+1+dimage->bayesGround+16);
    fprintf(dimage->opooMet,", %d zenith, %d FHD",14+4*metric->nRH+1+dimage->bayesGround+17,14+4*metric->nRH+1+dimage->bayesGround+18);
    fprintf(dimage->opooMet,", %d niM2, %d niM2.1",14+4*metric->nRH+1+dimage->bayesGround+19,14+4*metric->nRH+1+dimage->bayesGround+20);
    fprintf(dimage->opooMet,", %d meanNoise, %d noiseStdev, %d noiseThresh",14+4*metric->nRH+1+dimage->bayesGround+21,14+4*metric->nRH+1+dimage->bayesGround+22,14+4*metric->nRH+1+dimage->bayesGround+23);
    fprintf(dimage->opooMet,", %d FHDhist, %d FHDcan, %d FHDcanHist",14+4*metric->nRH+1+dimage->bayesGround+24,14+4*metric->nRH+1+dimage->bayesGround+25,14+4*metric->nRH+1+dimage->bayesGround+26);
    fprintf(dimage->opooMet,", %d FHDcanGauss, %d FHDcanGhist,",14+4*metric->nRH+1+dimage->bayesGround+27,14+4*metric->nRH+1+dimage->bayesGround+28);
    for(i=0;i<metric->laiBins;i++)fprintf(dimage->opooMet," %d tLAI%dt%d,",14+4*metric->nRH+1+dimage->bayesGround+26+i+3,i*(int)dimage->laiRes,(i+1)*(int)dimage->laiRes);
    for(i=0;i<metric->laiBins;i++)fprintf(dimage->opooMet," %d gLAI%dt%d,",14+4*metric->nRH+1+dimage->bayesGround+26+i+3+metric->laiBins,i*(int)dimage->laiRes,(i+1)*(int)dimage->laiRes);
    for(i=0;i<metric->laiBins;i++)fprintf(dimage->opooMet," %d hgLAI%dt%d,",14+4*metric->nRH+1+dimage->bayesGround+26+i+3+2*metric->laiBins,i*(int)dimage->laiRes,(i+1)*(int)dimage->laiRes);
    for(i=0;i<metric->laiBins;i++)fprintf(dimage->opooMet," %d hiLAI%dt%d,",14+4*metric->nRH+1+dimage->bayesGround+26+i+3+3*metric->laiBins,i*(int)dimage->laiRes,(i+1)*(int)dimage->laiRes);
    for(i=0;i<metric->laiBins;i++)fprintf(dimage->opooMet," %d hmLAI%dt%d,",14+4*metric->nRH+1+dimage->bayesGround+26+i+3+4*metric->laiBins,i*(int)dimage->laiRes,(i+1)*(int)dimage->laiRes);
    //for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet,", %d LmomGauss%d",14+4*metric->nRH+1+dimage->bayesGround+21+i,i+1);
    //for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet,", %d LmomInfl%d",14+4*metric->nRH+1+dimage->bayesGround+21+metric->nLm+i,i+1);
    //for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet,", %d LmomMax%d",14+4*metric->nRH+1+dimage->bayesGround+21+2*metric->nLm+i,i+1);
    //for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet,", %d LmomReal%d",14+4*metric->nRH+1+dimage->bayesGround+21+3*metric->nLm+i,i+1);
    fprintf(dimage->opooMet,"\n");
  }

  if(dimage->gediIO.gFit->nGauss>dimage->maxGauss)fprintf(stderr,"More Gaussians than header entries %d\n",dimage->gediIO.gFit->nGauss);

  /*fitted Gaussians*/
  if(dimage->writeGauss){
    fprintf(dimage->opooGauss,"%d %d",numb,dimage->gediIO.gFit->nGauss);
    for(i=0;i<dimage->gediIO.gFit->nGauss;i++){
      if((dimage->gediIO.gFit->gPar[3*i]>=0.0)&&(dimage->gediIO.gFit->gPar[3*i+1]>=0.0)&&(dimage->gediIO.gFit->gPar[3*i+2]>=0.0)){
        fprintf(dimage->opooGauss," %f %f %f",dimage->gediIO.gFit->gPar[3*i],dimage->gediIO.gFit->gPar[3*i+1],dimage->gediIO.gFit->gPar[3*i+2]);
      }else fprintf(dimage->opooGauss," ? ? ?");
    }
    for(i=dimage->gediIO.gFit->nGauss;i<dimage->maxGauss;i++)fprintf(dimage->opooGauss," ? ? ?");
    fprintf(dimage->opooGauss," %s\n",inNamen);
  }

  /*waveform metrics*/
  if(data->useID==0)fprintf(dimage->opooMet,"%d",numb);
  else              fprintf(dimage->opooMet,"%s",data->waveID);
  fprintf(dimage->opooMet," %.2f %.2f %.4f %.4f %.2f %.2f %.2f %.4f %.2f %.4f %.2f %.2f",data->gElev,data->tElev,data->slope,\
    data->cov,metric->gHeight,metric->maxGround,metric->inflGround,metric->tElev,metric->bElev,metric->cov,metric->leExt,metric->teExt);
  for(i=0;i<metric->nRH;i++)fprintf(dimage->opooMet," %.2f",metric->rh[i]);
  for(i=0;i<metric->nRH;i++)fprintf(dimage->opooMet," %.2f",metric->rhMax[i]);
  for(i=0;i<metric->nRH;i++)fprintf(dimage->opooMet," %.2f",metric->rhInfl[i]);
  for(i=0;i<metric->nRH;i++)fprintf(dimage->opooMet," %.2f",metric->rhReal[i]);
  fprintf(dimage->opooMet," %s",inNamen);
  if(dimage->bayesGround)fprintf(dimage->opooMet," %.2f",metric->bayGround);
  fprintf(dimage->opooMet," %.3f %.3f %.3f %.3f",metric->covHalfG,metric->covHalfM,metric->covHalfI,metric->covHalfB);
  fprintf(dimage->opooMet," %f %f",data->pSigma,data->fSigma);
  if(dimage->noise.linkNoise)fprintf(dimage->opooMet," %f %f",dimage->noise.linkM,dimage->noise.linkCov);
  else                 fprintf(dimage->opooMet," ? ?");
  if(dimage->coord2dp)fprintf(dimage->opooMet," %.2f %.2f",data->lon,data->lat);
  else                fprintf(dimage->opooMet," %.10f %.10f",data->lon,data->lat);
  fprintf(dimage->opooMet," %f %f %f",data->gLap,data->gMinimum,data->gInfl); 
  fprintf(dimage->opooMet," %f %f",metric->totE,metric->blairSense);
  fprintf(dimage->opooMet," %f %f",data->pointDense,data->beamDense);
  fprintf(dimage->opooMet," %f %f",data->zen,metric->FHD);
  fprintf(dimage->opooMet," %f %f",metric->niM2,metric->niM21);
  fprintf(dimage->opooMet," %f %f %f",dimage->gediIO.den->meanN,(dimage->gediIO.den->thresh-dimage->gediIO.den->meanN)/dimage->gediIO.den->threshScale,dimage->gediIO.den->thresh);
  fprintf(dimage->opooMet," %f %f %f",metric->FHDhist,metric->FHDcan,metric->FHDcanH);
  fprintf(dimage->opooMet," %f %f",metric->FHDcanGauss,metric->FHDcanGhist);
  for(i=0;i<metric->laiBins;i++)fprintf(dimage->opooMet," %f",metric->tLAI[i]);
  for(i=0;i<metric->laiBins;i++)fprintf(dimage->opooMet," %f",metric->gLAI[i]);
  for(i=0;i<metric->laiBins;i++)fprintf(dimage->opooMet," %f",metric->hgLAI[i]);
  for(i=0;i<metric->laiBins;i++)fprintf(dimage->opooMet," %f",metric->hiLAI[i]);
  for(i=0;i<metric->laiBins;i++)fprintf(dimage->opooMet," %f",metric->hmLAI[i]);
  /*for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet," %f",metric->LmomGau[i]);
  for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet," %f",metric->LmomInf[i]);
  for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet," %f",metric->LmomMax[i]);
  for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet," %f",metric->LmomRea[i]);*/
  fprintf(dimage->opooMet,"\n");

  /*fitted wave if required*/
  if(dimage->writeFit){
    if(data->useID==0)sprintf(waveNamen,"%s.%d.fit",dimage->outRoot,numb);
    else              sprintf(waveNamen,"%s.%s.fit",dimage->outRoot,data->waveID);
    if((opoo=fopen(waveNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",waveNamen);
      exit(1);
    }

    fprintf(opoo,"# 1 elevation, 2 noised, 3 denoised, 4 processed, 5 original, 6 ground, 7 canopy");
    for(i=0;i<dimage->gediIO.gFit->nGauss;i++)fprintf(opoo,", %d Gauss %d",i+8,i+1);
    fprintf(opoo,"\n");
    fprintf(opoo,"# fSigma %f pSigma %f res %f\n",data->fSigma,data->pSigma,dimage->gediIO.res);
    if(dimage->coord2dp)fprintf(opoo,"# coord %.2f %.2f\n",data->lon,data->lat);
    else                fprintf(opoo,"# coord %.10f %.10f\n",data->lon,data->lat);
    fprintf(opoo,"# cover %f rhoG %f rhoC %f\n",data->cov,rhoG,rhoC);
    fprintf(opoo,"# ground %.2f slope %f\n",data->gElev,data->slope);
    for(i=0;i<data->nBins;i++){
      if(dimage->noRHgauss==0)fprintf(opoo,"%f %f %f %f %f",data->z[i],data->noised[i],denoised[i],processed[i],data->wave[data->useType][i]);
      else                    fprintf(opoo,"%f %f %f ? %f",data->z[i],data->noised[i],denoised[i],data->wave[data->useType][i]);
      if(dimage->gediIO.ground)fprintf(opoo," %f %f",data->ground[data->useType][i],data->wave[data->useType][i]-data->ground[data->useType][i]);
      else              fprintf(opoo," 0 0");
      for(j=0;j<dimage->gediIO.gFit->nGauss;j++)fprintf(opoo," %f",dimage->gediIO.gFit->gPar[j*3+1]*gauss((float)data->z[i],dimage->gediIO.gFit->gPar[3*j+2],dimage->gediIO.gFit->gPar[3*j]));
      fprintf(opoo,"\n");
    }
    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
    fprintf(stdout,"Wave to %s\n",waveNamen);
  }/*fitted wave if required*/

  return;
}/*writeResults*/


/*####################################################*/
/*align Gaussian centres with elevation*/

void alignElevation(double z0,double zEnd,float *gPar,int nGauss)
{
  int i=0;
  double mu=0;

  /*convert Gaussian ranges*/
  for(i=0;i<nGauss;i++){
    if(z0>zEnd)mu=z0-(double)gPar[3*i];
    else       mu=(double)gPar[3*i]-z0;
    gPar[3*i]=mu;
  }
  return;
}/*alignElevation*/


/*####################################################*/
/*calculate metrics*/

void findMetrics(metStruct *metric,float *gPar,int nGauss,float *denoised,float *wave,int nBins,double *z,control *dimage,dataStruct *data)
{
  int i=0,gInd=0;
  float tot=0;
  float *mu=NULL,*sig=NULL;
  float *energy=NULL,*A=NULL;
  float gaussCover(float *,int,float *,float *,int,int);
  double gaussianGround(float *,float *,int *,int,float);
  double maxGround(float *,double *,int);
  double inflGround(float *,double *,int);
  double bayesGround(float *,int,control *,metStruct *,double *,dataStruct *);
  float *blankRH(float,int *);
  float *smoothed=NULL,*canWave=NULL;
  float halfCover(float *,double *,int,double,float);
  void findSignalBounds(float *,double *,int,double *,double *,control *);
  void findWaveExtents(float *,double *,int,double,double,float *,float *);
  void setDenoiseDefault(denPar *);
  float *trueLAIprofile(float *,float *,double *,int,float,float,double,float,int *);
  float *gaussLAIprofile(float *,double *,int,float,float,double,float,int *,float,float,float,float);
  float *halfEnergyLAIprofile(float *,double *,int,float,float,double,float,int *);
  denPar den;


  /*total energy. The sqrt(2.pi) can be normalised out*/
  tot=0.0;
  if(nGauss>0){  /*leave NULL if no Gaussians found to prevent memory leaks*/
    energy=falloc(nGauss,"energy",0);
    mu=falloc(nGauss,"mu",0);
    A=falloc(nGauss,"A",0);
    sig=falloc(nGauss,"sigma",0);
  }
  for(i=0;i<nGauss;i++){
    mu[i]=gPar[3*i];
    A[i]=gPar[3*i+1];
    sig[i]=gPar[3*i+2];
    if(A[i]>=0.0){
      energy[i]=A[i]*sig[i];
      tot+=energy[i];
    }else energy[i]=0.0;
  }
  /*normalise energy*/
  for(i=0;i<nGauss;i++)energy[i]/=tot;

  /*waveform energy*/
  metric->totE=0.0;
  for(i=0;i<nBins;i++)metric->totE+=denoised[i]*dimage->gediIO.res;

  /*Blair sensitivity*/
  metric->blairSense=findBlairSense(data,&dimage->gediIO);

  /*smooth waveform for fonding ground by max and inflection*/
  setDenoiseDefault(&den);
  den.varNoise=0;
  den.meanN=0;
  den.thresh=0;
  den.sWidth=0.76*3.0/4.0;  /*according to Bryan Blair*/
  den.noiseTrack=0;
  smoothed=processFloWave(denoised,nBins,&den,1.0);

  /*ground by Gaussian fit*/
  if(dimage->noRHgauss==0)metric->gHeight=gaussianGround(energy,mu,&gInd,nGauss,tot);
  else                    metric->gHeight=-1.0;

  /*canopy cover*/
  metric->cov=gaussCover(denoised,nBins,mu,energy,nGauss,gInd);

  /*ground by maximum*/
  metric->maxGround=maxGround(smoothed,z,nBins);

  /*ground by inflection*/
  metric->inflGround=inflGround(smoothed,z,nBins);

  /*rh metrics with Gaussian ground*/
  if(dimage->noRHgauss==0)metric->rh=findRH(denoised,z,nBins,metric->gHeight,dimage->rhRes,&metric->nRH);
  else                    metric->rh=blankRH(dimage->rhRes,&metric->nRH);

  /*rh metrics with maximum ground*/
  metric->rhMax=findRH(denoised,z,nBins,metric->maxGround,dimage->rhRes,&metric->nRH);

  /*rh metrics with inflection ground*/
  metric->rhInfl=findRH(denoised,z,nBins,metric->inflGround,dimage->rhRes,&metric->nRH);

  /*rh metrics with real ground, if we have the ground*/
  if(dimage->gediIO.ground||data->demGround){
    if(dimage->noise.linkNoise)metric->rhReal=findRH(data->wave[data->useType],z,nBins,data->gElev,dimage->rhRes,&metric->nRH);  /*original was noiseless*/
    else                 metric->rhReal=findRH(denoised,z,nBins,data->gElev,dimage->rhRes,&metric->nRH);  /*origina was noisy*/
  }else{
    metric->rhReal=falloc(metric->nRH,"rhReal",0);
    for(i=0;i<metric->nRH;i++)metric->rhReal[i]=-1.0;
  }

  /*foliage height diversity*/
  metric->FHD=foliageHeightDiversity(denoised,nBins);
  metric->FHDhist=foliageHeightDiversityHist(denoised,nBins,dimage->fhdHistRes);
  /*from ground removed canopy*/
  if(dimage->gediIO.ground)canWave=subtractGroundFromCan(data->wave[data->useType],data->ground[data->useType],nBins);
  else                     canWave=NULL;  /*no ground estimate. Leave blank*/
  metric->FHDcanH=foliageHeightDiversityHist(canWave,nBins,dimage->fhdHistRes);
  TIDY(canWave);
  /*from Gaussian removed canopy*/
  if(nGauss>0)canWave=subtractGaussFromCan(denoised,nBins,mu[gInd],A[gInd],sig[gInd],z);
  else        canWave=NULL;  /*no Gaussian ground estimate*/
  metric->FHDcanGauss=foliageHeightDiversity(canWave,nBins);
  metric->FHDcanGhist=foliageHeightDiversityHist(canWave,nBins,dimage->fhdHistRes);
  TIDY(canWave);


  /*lai profiles*/
  if(data->ground)metric->tLAI=trueLAIprofile(data->wave[data->useType],data->ground[data->useType],z,nBins,dimage->laiRes,dimage->rhoRatio,data->gElev,dimage->maxLAIh,&metric->laiBins);
  else metric->tLAI=trueLAIprofile(data->wave[data->useType],NULL,z,nBins,dimage->laiRes,dimage->rhoRatio,data->gElev,dimage->maxLAIh,&metric->laiBins);
  if(sig){
    metric->gLAI=gaussLAIprofile(denoised,z,nBins,dimage->laiRes,dimage->rhoRatio,metric->gHeight,dimage->maxLAIh,&metric->laiBins,mu[gInd],A[gInd],sig[gInd],dimage->gediIO.res);
  }else{
    metric->gLAI=trueLAIprofile(data->wave[data->useType],NULL,z,nBins,dimage->laiRes,dimage->rhoRatio,data->gElev,dimage->maxLAIh,&metric->laiBins);
  }
  metric->hgLAI=halfEnergyLAIprofile(denoised,z,nBins,dimage->laiRes,dimage->rhoRatio,metric->gHeight,dimage->maxLAIh,&metric->laiBins);
  metric->hiLAI=halfEnergyLAIprofile(denoised,z,nBins,dimage->laiRes,dimage->rhoRatio,metric->inflGround,dimage->maxLAIh,&metric->laiBins);
  metric->hmLAI=halfEnergyLAIprofile(denoised,z,nBins,dimage->laiRes,dimage->rhoRatio,metric->maxGround,dimage->maxLAIh,&metric->laiBins);


  /*signal start and end*/
  findSignalBounds(denoised,z,nBins,&metric->tElev,&metric->bElev,dimage);

  /*Lefsky's leading and trailing edge extents*/
  findWaveExtents(denoised,z,nBins,metric->tElev,metric->bElev,&metric->leExt,&metric->teExt);

  /*L moments*/
  /*metric->nLm=4;
  metric->LmomGau=waveLmoments(metric->rh,metric->nRH,dimage->rhRes,metric->nLm);
  metric->LmomRea=waveLmoments(metric->rhReal,metric->nRH,dimage->rhRes,metric->nLm);
  metric->LmomInf=waveLmoments(metric->rhInfl,metric->nRH,dimage->rhRes,metric->nLm);
  metric->LmomMax=waveLmoments(metric->rhMax,metric->nRH,dimage->rhRes,metric->nLm);*/

  /*cover estimates using Bryan's half feature method*/
  metric->covHalfG=halfCover(denoised,z,nBins,metric->gHeight,dimage->rhoRatio);
  metric->covHalfM=halfCover(denoised,z,nBins,metric->maxGround,dimage->rhoRatio);
  metric->covHalfI=halfCover(denoised,z,nBins,metric->inflGround,dimage->rhoRatio);

  /*Wenge Ni's metrics*/
  metric->niM2=niMetric(denoised,z,nBins,dimage->gediIO.res,metric->gHeight,2.0);
  metric->niM21=niMetric(denoised,z,nBins,dimage->gediIO.res,metric->gHeight,2.1);

  /*bayesian ground finding*/
  if(dimage->bayesGround){
    metric->bayGround=bayesGround(wave,nBins,dimage,metric,z,data);
    metric->covHalfB=halfCover(denoised,z,nBins,metric->bayGround,dimage->rhoRatio);
  }


  /*tidy up arrays*/
  TIDY(A);
  TIDY(mu);
  TIDY(sig);
  TIDY(energy);
  TIDY(smoothed);
  return;
}/*findMetrics*/


/*####################################################*/
/*LAI profile using reflected half energy*/

float *halfEnergyLAIprofile(float *denoised,double *z,int nBins,float laiRes,float rhoRatio,double gElev,float maxLAIh,int *laiBins)
{
  int i=0,gBin=0,mBin=0;
  float *LAI=NULL;
  float *canopy=NULL;
  float totG=0;
  double sep=0,minSep=0;

  /*Ground energy and bin*/
  minSep=1000000.0;
  totG=0.0;
  for(i=0;i<nBins;i++){
    if(z[i]<gElev)totG+=denoised[i];
    sep=fabs(z[i]-gElev);
    if(sep<=minSep){
      minSep=sep;
      gBin=i;
    }
  }
  totG*=2.0;

  /*canopy wave*/
  canopy=falloc(nBins,"canopy waveform",0);
  for(i=0;i<nBins;i++)canopy[i]=0.0;
  for(i=gBin;i>=0;i--){
    canopy[i]=denoised[i];
    /*subtract mirrored ground*/
    mBin=2*gBin-i;
    if((mBin>=0)&&(mBin<nBins)){
      canopy[i]-=denoised[mBin];
      if(canopy[i]<0.0)canopy[i]=0.0;
    }
  }

  /*lai profile*/
  LAI=findLAIprofile(canopy,totG,nBins,laiRes,laiBins,gElev,rhoRatio,z,maxLAIh);

  /*tidy up*/
  TIDY(canopy);
  return(LAI);
}/*halfEnergyLAIprofile*/


/*####################################################*/
/*LAI profile subtracting fitted Gaussian*/

float *gaussLAIprofile(float *denoised,double *z,int nBins,float laiRes,float rhoRatio,double gHeight,float maxLAIh,int *laiBins,float mu,float A,float sig,float res)
{
  int i=0;
  float *LAI=NULL,totG=0;
  float *canopy=NULL;

  /*canopy waveform*/
  canopy=falloc(nBins,"canopy only wave",0);
  for(i=0;i<nBins;i++){
    if(z[i]>=gHeight){
      canopy[i]=denoised[i]-A*(float)gaussian(z[i],(double)sig,(double)mu);
      if(canopy[i]<0.0)canopy[i]=0.0;
    }else canopy[i]=0.0;
  }

  /*ground energy*/
  totG=A*sig*sqrt(2.0*M_PI)/res;

  /*lai profile*/
  LAI=findLAIprofile(canopy,totG,nBins,laiRes,laiBins,gHeight,rhoRatio,z,maxLAIh);

  /*tidy up*/
  TIDY(canopy);
  return(LAI);
}/*gaussLAIprofile*/


/*####################################################*/
/*true LAI profile*/

float *trueLAIprofile(float *wave,float *ground,double *z,int nBins,float laiRes,float rhoRatio,double gElev,float maxLAIh,int *laiBins)
{
  int i=0;
  float *tLAI=NULL;
  float *canopy=NULL,totG=0;

  /*is there a ground to do this with?*/
  if(ground){
    /*Calculate canopy part*/
    canopy=falloc(nBins,"canopy part",0);
    totG=0.0;
    for(i=0;i<nBins;i++){
      totG+=ground[i];
      canopy[i]=wave[i]-ground[i];
    }

    /*lai profile*/
    tLAI=findLAIprofile(canopy,totG,nBins,laiRes,laiBins,gElev,rhoRatio,z,maxLAIh);
  }else{  /*no ground, leave blank*/
    *laiBins=(int)(maxLAIh/laiRes+0.5)+1;
    tLAI=falloc(*laiBins,"True LAI profile",0);
    for(i=0;i<*laiBins;i++)tLAI[i]=-1.0;
  }

  /*tidy up*/
  TIDY(canopy);
  return(tLAI);
}/*trueLAIprofile*/


/*####################################################*/
/*LAI profile from canopy waveform*/

float *findLAIprofile(float *canopy,float totG,int nBins,float laiRes,int *laiBins,double gElev,float rhoRatio,double *z,float maxLAIh)
{
  int i=0,bin=0;
  float *tLAI=NULL;
  float gap=0,cumul=0;
  float *lngap=NULL;
  float *laiProf=NULL;
  float totC=0;
  float G=0;           /*Ross-G function*/

  /*leaf angle distribution*/
  G=0.6;   /*spherical distribution*/

  /*total energies*/
  totC=0.0;
  for(i=0;i<nBins;i++)totC+=canopy[i];

  /*make gap profile*/
  lngap=falloc(nBins,"apparent foliage profile",0);
  cumul=0.0;
  for(i=0;i<nBins;i++){
    cumul+=canopy[i];
    gap=1.0-(cumul/(totC+rhoRatio*totG));
    lngap[i]=log(gap);
  }

  /*lai profile*/
  laiProf=falloc(nBins,"apparent foliage profile",0);
  for(i=0;i<nBins;i++){
    if((i>0)&&(i<(nBins-1)))laiProf[i]=(-2.0*(lngap[i+1]-lngap[i-1])/2.0)/G;
    else                    laiProf[i]=0.0;
  }
  TIDY(lngap);

  /*allocate and set blank*/
  *laiBins=(int)(maxLAIh/laiRes+0.5)+1;
  tLAI=falloc(*laiBins,"True LAI profile",0);
  for(i=0;i<(*laiBins);i++)tLAI[i]=0.0;

  /*bin up*/
  for(i=0;i<nBins;i++){
    bin=(int)((float)(z[i]-gElev)/laiRes);
    if(bin>=*laiBins)bin=(*laiBins)-1;
    else if(bin<0)bin=0;
    tLAI[bin]+=laiProf[i];
  }

  /*tidy up*/
  TIDY(laiProf);
  return(tLAI);
}/*findLAIprofile*/


/*####################################################*/
/*canopy cover, taking ground as double energy beneath*/

float halfCover(float *wave,double *z,int nBins,double gElev,float rhoRatio)
{
  int i=0;
  float canE=0,grE=0;
  float cov=0;

  canE=grE=0.0;
  for(i=0;i<nBins;i++){
    if(wave[i]>0.0){
      if(z[i]<=gElev)grE+=wave[i];
      else           canE+=wave[i];
    }
  }
  canE-=grE;  /*as it's half the energy*/
  grE*=2.0;   /*as it's the half energy*/
  if(grE<=(canE+grE))cov=canE/(canE+grE*rhoRatio);
  else        cov=0.0;

  return(cov);
}/*halfCover*/


/*####################################################*/
/*Bayseian ground finding*/

double bayesGround(float *wave,int nBins,control *dimage,metStruct *metric,double *z,dataStruct *data)
{
  int i=0,start=0,end=0,dir=0;
  float contN=0,prob=0;
  float *processed=NULL;
  float groundProb(float,float,float);
  float height=0;
  double bayGround=0;
  denPar *den=NULL;    /*denoising structure*/
  void gaussProps(float *,int,float,float,double *,float *,float *,float *);
  void  setDenoiseDefault(denPar *);
  void alignElevation(double,double,float *,int);

  /*allocate*/
  metric->nBgr=6;
  if(!(metric->bGr=(bGround *)calloc(metric->nBgr,sizeof(bGround)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*for first estimate, last return above noise*/
  if(z[0]>z[nBins-1]){
    start=nBins-2;
    end=2;
    dir=-1;
  }else{
    start=2;
    end=nBins-2;
    dir=1;
  }
  i=start;
  while(((dir==1)&&(i<=end))||((dir==-1)&&(i>=end))){
    if(wave[i]>0.0){
      metric->bGr[0].gHeight=z[i]+(double)(data->pSigma*2.0);
      metric->bGr[0].slope=0.0;
      metric->bGr[0].cov=1.0;
      break;
    }
    i+=dir;
  }

  /*for others, make an array*/
  if(!(den=(denPar *)calloc(metric->nBgr,sizeof(denPar)))){
    fprintf(stderr,"error waveStruct allocation.\n");
    exit(1);
  }
  for(i=1;i<metric->nBgr;i++){
    setDenoiseDefault(&(den[i]));
    den[i].varNoise=1;
    den[i].gWidth=dimage->gediIO.den->sWidth;
    if(i<=3)den[i].threshScale=3.0;
    else    den[i].threshScale=5.0;
    den[i].noiseTrack=1;
    den[i].fitGauss=1;
    if(i<=3)den[i].sWidth=(float)i*data->pSigma/1.6;
    else    den[i].sWidth=(float)(i-3)*data->pSigma/1.6;
  }

  /*loop over smoothing widths*/
  for(i=1;i<metric->nBgr;i++){
    /*fit Gaussians*/
    processed=processFloWave(wave,nBins,&den[i],1.0);
    alignElevation(z[0],z[nBins-1],den[i].gPar,den[i].nGauss);

    /*get parameters*/
    gaussProps(den[i].gPar,den[i].nGauss,dimage->gediIO.fSigma,sqrt(data->pSigma*data->pSigma+den[i].sWidth*den[i].sWidth),&(metric->bGr[i].gHeight),&(metric->bGr[i].slope),&(metric->bGr[i].cov),&height);

    TIDY(den[i].gPar);
    TIDY(processed);
  }/*smoothing width loop*/

  contN=0.0;
  bayGround=0.0;
  for(i=0;i<metric->nBgr;i++){
    prob=groundProb(metric->bGr[i].slope,metric->bGr[i].cov,height);
    bayGround+=(double)(prob*metric->bGr[i].gHeight);
    contN+=prob;
  }
  if(contN>0.0)bayGround/=(double)contN;
  TIDY(den);

  return(bayGround);
}/*bayesGround*/


/*####################################################*/
/*ground elevation probability*/

float groundProb(float slope,float cov,float height)
{
  float prob=0;

  if(slope<5.0)      prob=0.5+0.1*slope;
  else if(slope<30.0)prob=1.0;
  else               prob=slope/-30.0+2.0;
  /*do not let it be zero or negative*/
  if(prob<=0.0)prob=0.0001;

  if(height>80.0)     prob*=0.000001;
  else if(height>60.0)prob*=3.0-height/30.0;

  /*do not let it be zero or negative*/
  if(prob<=0.0)prob=0.0001;

  return(prob);
}/*groundProb*/


/*####################################################*/
/*Gaussian ground properties*/

void gaussProps(float *gPar,int nGauss,float fSigma,float pSigma,double *gHeight,float *slope,float *cov,float *height)
{
  int i=0,use=0;
  float totalE=0;
  float minH=0,maxH=0;
  float sig=0;

  totalE=0.0;
  minH=1000000.0;
  maxH=-1000000.0;
  for(i=0;i<nGauss;i++)totalE+=gPar[3*i+1]*gPar[3*i+2];
  for(i=0;i<nGauss;i++){
    sig=gPar[3*i+2];

    if(((gPar[3*i+1]*gPar[3*i+2])>=(totalE*0.001))&&(gPar[3*i]<minH)){
      minH=gPar[3*i];
      use=i;
    }
    if(((gPar[3*i+1]*gPar[3*i+2])>=(totalE*0.001))&&(gPar[3*i]>(maxH-sig))){
      maxH=gPar[3*i]+sig;
      use=i;
    }
  }

  *gHeight=minH;
  *height=maxH-minH;
  sig=gPar[3*use+2];
  if(sig>pSigma)*slope=atan2(sqrt(sig*sig-pSigma*pSigma),fSigma)*180.0/M_PI;
  else          *slope=0.0;
  *cov=1.0-gPar[3*use+1]*gPar[3*use+2]/totalE;

  return;
}/*gaussProps*/


/*####################################################*/
/*Lefsky's leading and trailing edge extents*/

void findWaveExtents(float *processed,double *z,int nBins,double tElev,double bElev,float *leExt,float *teExt)
{
  int i=0,contN=0;
  float mean=0;

  /*determine mean*/
  contN=0;
  mean=0.0;
  for(i=0;i<nBins;i++){
    if((z[i]>=bElev)&&(z[i]<=tElev)){
      mean+=processed[i];
      contN++;
    }
  }
  if(contN>0)mean/=(float)contN;

  /*leading edge*/
  for(i=0;i<nBins;i++){
    if((z[i]>=bElev)&&(z[i]<=tElev)){
      if(processed[i]>=mean){
        *leExt=(float)(tElev-z[i]);
        break;
      }
    }
  }

  /*trailing edge*/
  for(i=nBins-1;i>=0;i--){
    if((z[i]>=bElev)&&(z[i]<=tElev)){
      if(processed[i]>=mean){
        *teExt=(float)(z[i]-bElev);
        break;
      }
    }
  }
  return;
}/*findWaveExtents*/


/*####################################################*/
/*signal top and bottom*/

void findSignalBounds(float *processed,double *z,int nBins,double *tElev,double *bElev,control *dimage)
{
  int i=0;
  float cumul=0,totE=0;

  /*total energy*/
  totE=0.0;
  for(i=nBins-1;i>=0;i--)totE+=processed[i];

  /*find top*/
  cumul=0.0;
  for(i=0;i<nBins;i++){
    cumul+=processed[i]/totE;
    if(processed[i]>dimage->bThresh){
      *tElev=z[i];
      break;
    }
  }/*bin loop*/

  /*find bottom*/
  cumul=0.0;
  for(i=nBins-1;i>=0;i--){
    cumul+=processed[i]/totE;
    if(processed[i]>dimage->bThresh){
      *bElev=z[i];
      break;
    }
  }/*bin loop*/

  return;
}/*findSignalBounds*/


/*####################################################*/
/*blank RH metric array*/

float *blankRH(float rhRes,int *nRH)
{
  int i=0;
  float *rh=NULL;

  *nRH=(int)(100.0/rhRes)+1;
  rh=falloc(*nRH,"rh metrics",0);
  for(i=0;i<*nRH;i++)rh[i]=-1.0;

  return(rh);
}/*blankRH*/


/*####################################################*/
/*Gaussian canopy cover*/

float gaussCover(float *wave,int nBins,float *mu,float *energy,int nGauss,int gInd)
{
  int i=0;
  float cov=0,totE=0;


  totE=0.0;
  for(i=0;i<nGauss;i++)totE+=energy[i];

  if((gInd>0)&&(gInd<nGauss))cov=1.0-energy[gInd]/totE;
  return(cov);
}/*gaussCover*/


/*####################################################*/
/*Gaussian ground*/

double gaussianGround(float *energy,float *mu,int *gInd,int nGauss,float tot)
{
  int i=0;
  float thresh=0;
  double gHeight=0;

  (*gInd)=nGauss-1;
  thresh=0.001;  /*1% of energy*/
  gHeight=100000000.0;
  for(i=0;i<nGauss;i++){
    if((energy[i]>=thresh)&&(mu[i]<gHeight)){
      gHeight=mu[i];
      *gInd=i;
    }
  }

  return(gHeight);
}/*gaussianGround*/


/*####################################################*/
/*ground by inflection point*/

double inflGround(float *smoothed,double *z,int nBins)
{
  int i=0,dir=0;
  int start=0,end=0;
  float contN=0;
  float *d2x=NULL;
  double CofG=0; 
  char inFeat=0;

  /*determine direction of wave*/
  if(z[0]>z[nBins-1]){
    start=nBins-2;
    end=2;
    dir=-1;
  }else{
    start=2;
    end=nBins-2;
    dir=1;
  }

  /*determine derivatives*/
  d2x=falloc(nBins,"d2x",0);
  for(i=1;i<(nBins-1);i++)d2x[i]=2.0*smoothed[i]-(smoothed[i+1]+smoothed[i-1]);

  /*loop through looking for first two inflection points*/
  CofG=0.0;
  contN=0.0;
  inFeat=0;

  i=start;
  while(((dir==1)&&(i<=end))||((dir==-1)&&(i>=end))){
    if((d2x[i]<=d2x[i-1])&&(d2x[i]<d2x[i+1])){  /*minimum of the second derivative*/
      if(inFeat)break;
      else inFeat=1;
    }

    if(inFeat){
      CofG+=z[i]*(double)smoothed[i];
      contN+=smoothed[i];
    }
    i+=dir;
  }
  if(contN>0.0)CofG/=(double)contN;
  else         CofG=-1000.0;
  TIDY(d2x);

  return(CofG);
}/*inflGround*/


/*####################################################*/
/*ground by maximum*/

double maxGround(float *smoothed,double *z,int nBins)
{
  int i=0;
  double maxGround=0;

  maxGround=1000000.0;
  if(z[0]>z[nBins-1]){
    for(i=1;i<nBins-1;i++){
      if((smoothed[i]>=smoothed[i-1])&&(smoothed[i]>smoothed[i+1])){
        if(z[i]<maxGround)maxGround=z[i];
      }
    }
  }else{
    for(i=nBins-2;i>0;i--){
      if((smoothed[i]>=smoothed[i-1])&&(smoothed[i]>smoothed[i+1])){
        if(z[i]<maxGround)maxGround=z[i];
      }
    }
  }

  return(maxGround);
}/*maxGround*/


/*####################################################*/
/*read or copy LVIS L2 data*/

void setL2ground(dataStruct *data,int numb,control *dimage)
{
  uint64_t i=0;
  lvisL2struct *readLvisL2(char *);
  lvisL2struct *t=NULL;

  if(numb==0)dimage->lvisL2=readLvisL2(dimage->l2namen);

  /*for shorthand, set to a short pointer*/
  t=dimage->lvisL2;

  /*find matching L2 value*/
  for(i=0;i<t->numb;i++){
    if((t->lfid[i]==data->lfid)&&(t->shotN[i]==data->shotN)){
      data->gElev=t->zG[i];
      data->demGround=1;
      break;
    }
  }

  t=NULL;
  return;
}/*setL2ground*/

/*####################################################*/
/*read LVIS L2 data*/

lvisL2struct *readLvisL2(char *namen)
{
  int j=0,zGcol=0;
  uint64_t i=0;
  lvisL2struct *lvisL2=NULL;
  char line[5000];
  char *token=NULL;
  FILE *ipoo=NULL;

  /*allocate structures*/
  if(!(lvisL2=(lvisL2struct *)calloc(1,sizeof(lvisL2struct)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  lvisL2->numb=0;
  lvisL2->lfid=NULL;
  lvisL2->shotN=NULL;
  lvisL2->zG=NULL;

  /*open file*/
  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  /*count number of lines*/
  while(fgets(line,5000,ipoo)!=NULL)if(strncasecmp(line,"#",1))lvisL2->numb++;

  /*allocate arrays*/
  if(!(lvisL2->lfid=(uint32_t *)calloc(lvisL2->numb,sizeof(uint32_t)))){
    fprintf(stderr,"error in L2 lfid allocation.\n");
    exit(1);
  }
  if(!(lvisL2->shotN=(uint32_t *)calloc(lvisL2->numb,sizeof(uint32_t)))){
    fprintf(stderr,"error in L2 shotN allocation.\n");
    exit(1);
  }
  lvisL2->zG=falloc((int)lvisL2->numb,"L2 zG",0);

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  i=0;
  while(fgets(line,5000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      j=0;
      token=strtok(line," ");
      while(token){
        if(j==0)lvisL2->lfid[i]=atoi(token);
        else if(j==1)lvisL2->shotN[i]=atoi(token);
        else if(j==zGcol){
          lvisL2->zG[i]=atof(token);
          break;
        }
        token=strtok(NULL," ");
        j++;
      }
      i++;
    }else if(!strncasecmp(line,"# LFID",6)){
      j=0;
      token=strtok(line," ");
      while(token){
        if(!strncasecmp(token,"ZG",2))zGcol=j-1;  /*minus one to account for # at line start*/
        token=strtok(NULL," ");
        j++;
      }
    }
  }

  /*close*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(lvisL2);
}/*readLvisL2*/


/*####################################################*/
/*dummy function if photon counting not included*/

#ifndef USEPHOTON
void photonCountCloud(float *denoised,dataStruct *data,photonStruct *photonCount,char *outRoot,int numb)

{
  fprintf(stderr,"This has been compiled without photon counting functions\n");
  exit(1);
  return;
}/*photonCountCloud*/
#endif


/*####################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  void setDenoiseDefault(denPar *);
  void readPulse(denPar *);
  void writeHelp();

  /*allocate structures*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  if(!(dimage->gediIO.den=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  if(!(dimage->gediIO.gFit=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }


  /*defaults*/
  /*input/output*/
  dimage->gediIO.nFiles=1;
  dimage->gediIO.inList=chChalloc(dimage->gediIO.nFiles,"inList",0);
  dimage->gediIO.inList[0]=challoc(200,"inList",0);
  strcpy(&(dimage->gediIO.inList[0][0]),"/Users/stevenhancock/data/gedi/analysis/side_lobe/laselva/contLaselva.0.12.wave");
  strcpy(dimage->outRoot,"teastMetric");
  dimage->maxGauss=20;
  dimage->opooMet=NULL;
  dimage->opooGauss=NULL;
  dimage->hdfGedi=NULL;
  /*scan settings*/
  dimage->gediIO.pSigma=0.764331; /*pulse length*/
  dimage->gediIO.fSigma=5.5;      /*footprint width*/
  dimage->gediIO.res=0.15;

  /*switches*/
  dimage->writeFit=0;
  dimage->gediIO.ground=0;
  dimage->gediIO.useInt=0;
  dimage->gediIO.useCount=1;
  dimage->gediIO.useFrac=0;
  dimage->rhRes=5.0;
  dimage->laiRes=10.0;
  dimage->maxLAIh=30.0;
  dimage->bayesGround=0;
  dimage->noise.missGround=0;
  dimage->noise.linkNoise=0;
  dimage->noise.driftFact=0.0;
  dimage->gediIO.linkPsig=0.764331; /*pulse length*/
  dimage->gediIO.linkFsig=5.5;      /*footprint width*/
  dimage->noise.trueSig=5.0;
  dimage->noise.deSig=0.0; //0.1; //4.0*0.15/2.355;
  dimage->noise.bitRate=12;
  dimage->noise.maxDN=4096.0; //1.0/(dimage->pSigma*sqrt(2.0*M_PI));
  dimage->noise.minGap=0.0;
  dimage->noRHgauss=0;    /*do find RH metrics by Gaussian fitting*/
  dimage->renoiseWave=0;  /*do not denoise "truth"*/
  dimage->noise.newPsig=-1.0;   /*leave blank*/
  dimage->gediIO.dontTrustGround=0;  /*do trust ground in waveforms, if there*/
  dimage->readBinLVIS=0;      /*read ASCII rather than binary LVIS*/
  dimage->readHDFlvis=0;      /*read ASCII rather than HDF5 LVIS*/
  dimage->readHDFgedi=0;      /*read ASCII rather than HDF5 GEDI*/
  dimage->gediIO.readPsigma=1;       /*read pSigma from file*/
  dimage->coord2dp=1;         /*round up coords in output*/
  dimage->useBounds=0;        /*process all data provided*/
  dimage->writeGauss=0;       /*do not write Gaussian parameters*/

  /*set default denoising parameters*/
  setDenoiseDefault(dimage->gediIO.den);
  dimage->gediIO.den->meanN=0.0;  /*we haven't added noise yet*/
  dimage->gediIO.den->thresh=0.00000001;  /*tiny number as no noise yet*/
  dimage->gediIO.den->noiseTrack=0;
  dimage->gediIO.den->minWidth=0;
  dimage->gediIO.den->varNoise=0;
  dimage->gediIO.den->threshScale=1.5;
  dimage->gediIO.den->fitGauss=0;
  dimage->gediIO.den->psWidth=0.0;

  /*set default Gaussian fitting  parameters*/
  setDenoiseDefault(dimage->gediIO.gFit);
  dimage->gediIO.gFit->meanN=0.0;    /*no denoising here*/
  dimage->gediIO.gFit->thresh=0.000000005;    /*no denoising here*/
  dimage->gediIO.gFit->noiseTrack=0;    /*no denoising here*/
  dimage->gediIO.gFit->minWidth=0;    /*no denoising here*/
  dimage->gediIO.gFit->varNoise=0;    /*no denoising here*/
  dimage->gediIO.gFit->gWidth=1.2;
  dimage->gediIO.gFit->sWidth=0.0;
  dimage->gediIO.gFit->fitGauss=1;
  dimage->gediIO.gFit->minGsig=0.764331;
  /*noise parameters*/
  dimage->noise.meanN=0.0;
  dimage->noise.nSig=0.0;
  dimage->bThresh=0.001;
  dimage->noise.hNoise=0.0;
  dimage->noise.offset=94.0;
  /*LVIS data*/
  dimage->lvis.data=NULL;
  dimage->hdfLvis=NULL;
  /*LVIS level2 data*/
  dimage->readL2=0;   /*do not read L2*/
  /*photon counting*/
  dimage->ice2=0;             /*GEDI mode, rather than ICESat-2*/
  #ifdef USEPHOTON
  dimage->photonCount.designval=2.1;
  dimage->photonCount.prob=NULL;
  dimage->photonCount.pBins=0;
  dimage->photonCount.H=200.0;
  dimage->photonCount.noise_mult=0.1;
  #endif
  /*others*/
  rhoG=0.4;
  rhoC=0.57;
  dimage->rhoRatio=rhoC/rhoG;
  dimage->gTol=0.0;
  dimage->gediIO.nMessages=200;
  dimage->fhdHistRes=0.001;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        TTIDY((void **)dimage->gediIO.inList,dimage->gediIO.nFiles);
        dimage->gediIO.inList=NULL;
        dimage->gediIO.nFiles=1;
        dimage->gediIO.inList=chChalloc(dimage->gediIO.nFiles,"input name list",0);
        dimage->gediIO.inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->gediIO.inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->gediIO.inList,dimage->gediIO.nFiles);
        dimage->gediIO.inList=readInList(&dimage->gediIO.nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(dimage->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-level2",7)){
        checkArguments(1,i,argc,"-level2");
        dimage->readL2=1;
        strcpy(dimage->l2namen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-writeFit",9)){
        dimage->writeFit=1;
      }else if(!strncasecmp(argv[i],"-writeGauss",11)){
        dimage->writeGauss=1;
      }else if(!strncasecmp(argv[i],"-dcBias",7)){
        checkArguments(1,i,argc,"-dcBias");
        dimage->noise.meanN=dimage->noise.offset=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-hNoise",7)){
        checkArguments(1,i,argc,"-hNoise");
        dimage->noise.hNoise=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-nSig",5)){
        checkArguments(1,i,argc,"-nSig");
        dimage->noise.nSig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-seed",5)){
        checkArguments(1,i,argc,"-seed");
        srand(atoi(argv[++i]));
      }else if(!strncasecmp(argv[i],"-meanN",5)){
        checkArguments(1,i,argc,"-meanN");
        dimage->gediIO.den->meanN=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-thresh",6)){
        checkArguments(1,i,argc,"-thresh");
        dimage->gediIO.den->thresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-sWidth",7)){
        checkArguments(1,i,argc,"-sWidth");
        dimage->gediIO.den->sWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-psWidth",8)){
        checkArguments(1,i,argc,"-psWidth");
        dimage->gediIO.den->psWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-msWidth",8)){
        checkArguments(1,i,argc,"-msWidth");
        dimage->gediIO.den->msWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gWidth",7)){
        checkArguments(1,i,argc,"-gWidth");
        dimage->gediIO.gFit->gWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minGsig",8)){
        checkArguments(1,i,argc,"-minGsig");
        dimage->gediIO.gFit->minGsig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minWidth",9)){
        checkArguments(1,i,argc,"-minWidth");
        dimage->gediIO.den->minWidth=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-ground",7)){
        dimage->gediIO.ground=1;
      }else if(!strncasecmp(argv[i],"-varNoise",9)){
        dimage->gediIO.den->varNoise=1;
      }else if(!strncasecmp(argv[i],"-medNoise",9)){
        dimage->gediIO.den->medStats=1;
      }else if(!strncasecmp(argv[i],"-statsLen",9)){
        checkArguments(1,i,argc,"-statsLen");
        dimage->gediIO.den->statsLen=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-varScale",9)){
        checkArguments(1,i,argc,"-varScale");
        dimage->gediIO.den->varNoise=1;
        dimage->gediIO.den->threshScale=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noiseTrack",11)){
        dimage->gediIO.den->noiseTrack=1;
      }else if(!strncasecmp(argv[i],"-pFile",6)){
        checkArguments(1,i,argc,"-pFile");
        strcpy(dimage->gediIO.den->pNamen,argv[++i]);
        dimage->gediIO.den->deconGauss=0;
      }else if(!strncasecmp(argv[i],"-pSigma",7)){
        dimage->gediIO.den->pSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gold",5)){
        dimage->gediIO.den->deconMeth=0;
      }else if(!strncasecmp(argv[i],"-deconTol",9)){
        checkArguments(1,i,argc,"-deconTol");
        dimage->gediIO.den->deChang=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-preMatchF",10)){
        dimage->gediIO.den->preMatchF=1;
      }else if(!strncasecmp(argv[i],"-postMatchF",11)){
        dimage->gediIO.den->posMatchF=1;
      }else if(!strncasecmp(argv[i],"-useInt",7)){
        dimage->gediIO.useInt=1;
        dimage->gediIO.useCount=0;
        dimage->gediIO.useFrac=0;
      }else if(!strncasecmp(argv[i],"-useFrac",8)){
        dimage->gediIO.useInt=0;
        dimage->gediIO.useCount=0;
        dimage->gediIO.useFrac=1;
      }else if(!strncasecmp(argv[i],"-rhRes",6)){
        checkArguments(1,i,argc,"-rhRes");
        dimage->rhRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-laiRes",7)){
        checkArguments(1,i,argc,"-laiRes");
        dimage->laiRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-laiH",5)){
        checkArguments(1,i,argc,"-laiH");
        dimage->maxLAIh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-linkNoise",10)){
        checkArguments(2,i,argc,"-linkNoise");
        dimage->noise.linkNoise=1;
        dimage->noise.linkM=atof(argv[++i]);
        dimage->noise.linkCov=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-trueSig",8)){
        checkArguments(1,i,argc,"-trueSig");
        dimage->noise.trueSig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-missGround",11)){
        dimage->noise.missGround=1;
      }else if(!strncasecmp(argv[i],"-minGap",7)){
        checkArguments(1,i,argc,"-minGap");
        dimage->noise.missGround=1;
        dimage->noise.minGap=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-bayesGround",13)){
        dimage->bayesGround=1;
      }else if(!strncasecmp(argv[i],"-rhoG",5)){
        checkArguments(1,i,argc,"-rhoG");
        rhoG=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-rhoC",5)){
        checkArguments(1,i,argc,"-rhoC");
        rhoC=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxDN",6)){
        checkArguments(1,i,argc,"-maxDN");
        dimage->noise.maxDN=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-bitRate",8)){
        checkArguments(1,i,argc,"-bitRate");
        dimage->noise.bitRate=(char)atoi(argv[++i]);
        dimage->noise.maxDN=pow(2.0,(float)dimage->noise.bitRate);
      }else if(!strncasecmp(argv[i],"-gTol",5)){
        checkArguments(1,i,argc,"-gTol");
        dimage->gTol=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noRHgauss",10)){
        dimage->noRHgauss=1;
      }else if(!strncasecmp(argv[i],"-renoise",8)){
       dimage->renoiseWave=1;
      }else if(!strncasecmp(argv[i],"-newPsig",8)){
        checkArguments(1,i,argc,"-newPsig");
        dimage->noise.newPsig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-oldPsig",8)){
        checkArguments(1,i,argc,"-oldPsig");
        dimage->gediIO.pSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-dontTrustGround",16)){
        dimage->gediIO.dontTrustGround=1;
      }else if(!strncasecmp(argv[i],"-readBinLVIS",12)){
        dimage->readBinLVIS=1;
      }else if(!strncasecmp(argv[i],"-readHDFlvis",12)){
        dimage->readHDFlvis=1;
        dimage->gediIO.readPsigma=0;
      }else if(!strncasecmp(argv[i],"-readHDFgedi",12)){
        dimage->readHDFgedi=1;
      }else if(!strncasecmp(argv[i],"-forcePsigma",12)){
        dimage->gediIO.readPsigma=0;
      }else if(!strncasecmp(argv[i],"-noRoundCoord",13)){
        dimage->coord2dp=0;
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(4,i,argc,"-bounds");
        dimage->useBounds=1;
        dimage->minX=atof(argv[++i]);
        dimage->minY=atof(argv[++i]);
        dimage->maxX=atof(argv[++i]);
        dimage->maxY=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-linkPsig",9)){
        checkArguments(1,i,argc,"-linkPsig");
        dimage->gediIO.linkPsig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-linkFsig",9)){
        checkArguments(1,i,argc,"-linkFsig");
        dimage->gediIO.linkFsig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fhdHistRes",11)){
        checkArguments(1,i,argc,"-fhdHistRes");
        dimage->fhdHistRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-addDrift",9)){
        checkArguments(1,i,argc,"-addDrift");
        dimage->noise.driftFact=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-varDrift",9)){
        dimage->gediIO.den->corrDrift=1;
        dimage->gediIO.den->varDrift=1;
      }else if(!strncasecmp(argv[i],"-driftFac",9)){
        checkArguments(1,i,argc,"-driftFac");
        dimage->gediIO.den->corrDrift=1;
        dimage->gediIO.den->varDrift=0;
        dimage->gediIO.den->fixedDrift=atof(argv[++i]);
      #ifdef USEPHOTON
      }else if(!strncasecmp(argv[i],"-photonCount",12)){
        dimage->ice2=1;
      }else if(!strncasecmp(argv[i],"-nPhotons",9)){
        checkArguments(1,i,argc,"-nPhotons");
        dimage->photonCount.designval=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-photonWind",11)){
        checkArguments(1,i,argc,"-photonWind");
        dimage->photonCount.H=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noiseMult",10)){
        checkArguments(1,i,argc,"-noiseMult");
        dimage->photonCount.noise_mult=atof(argv[++i]);
      #endif
      }else if(!strncasecmp(argv[i],"-help",5)){
        writeHelp();
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry gediRat -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }

  /*read deconvolution pulse if needed*/
  if(dimage->gediIO.den->preMatchF||dimage->gediIO.den->preMatchF||dimage->gediIO.den->deconMeth>=0)readPulse(dimage->gediIO.den); 
  if((!dimage->gediIO.ground)&&(dimage->noise.missGround)){
    fprintf(stderr,"Noise option conflict. Cannot use missGround without ground\n");
    exit(1);
  }

  return(dimage);
}/*readCommands*/


/*###########################################################*/
/*write help statement*/

void writeHelp()
{
  fprintf(stdout,"\n#########################\nProgram to calculate GEDI waveform metrics\n#########################\n\
\nInput output\n\
-input name;      waveform  input filename\n\
-outRoot name;    output filename root\n\
-inList list;     input file list for multiple files\n\
-writeFit;        write fitted waveform\n\
-writeGauss;      write Gaussian parameters\n\
-readBinLVIS;     input is an LVIS binary file\n\
-readHDFlvis;     read LVIS HDF5 input\n\
-readHDFgedi;     read GEDI simulator HDF5 input\n\
-level2 name;     level2 filename for LVIS ZG\n\
-bounds minX minY maxX maxY;    only analyse data within bounds\n\
\nSwitches\n\
-ground;          read true ground from file\n\
-useInt;          use discrete intensity instead of count\n\
-useFrac;         use fractional hits rather than counts\n\
-rhRes r;         percentage energy resolution of RH metrics\n\
-laiRes res;      lai profile resolution in metres\n\
-laiH h;          height to calculate LAI to\n\
-noRHgauss;       do not fit Gaussians\n\
-gTol tol;        ALS ground tolerance. Used to calculate slope.\n\
-fhdHistRes res;  waveform intensity resolution to use when calculating FHD from histograms\n\
-forcePsigma;     do not read pulse sigma from file\n\
-bayesGround;     use Bayseian ground finding\n\
-dontTrustGround; don't trust ground in waveforms, if included\n\
-noRoundCoord;    do not round up coords when outputting\n\
\nAdding noise:\n\
-dcBias n;        mean noise level\n\
-nSig sig;        noise sigma\n\
-seed n;          random number seed\n\
-hNoise n;        hard threshold noise as a fraction of integral\n\
-linkNoise linkM cov;     apply Gaussian noise based on link margin at a cover\n\
-linkFsig sig;    footprint width to use when calculating and applying signal noise\n\
-linkPsig sig;    pulse width to use when calculating and applying signal noise\n\
-trueSig sig;     true sigma of background noise\n\
-bitRate n;       digitisation bit rate\n\
-maxDN max;       maximum DN\n\
-renoise;         remove noise from truth before applying new noise level\n\
-newPsig sig;     new value for pulse width, when lengthening pulse\n\
-oldPsig sig;     old value for pulse width if not defined in waveform file, when lengthening pulse\n\
-addDrift xi;     apply detector background drift\n\
-missGround;      assume ground is missed to assess RH metrics\n\
-minGap gap;      delete signal beneath min detectable gap fraction\n");
  #ifdef USEPHOTON
  fprintf(stdout,"\nPhoton counting\n\
-photonCount;     output point cloud from photon counting\n\
-nPhotons n;      mean number of photons\n\
-photonWind x;    window length for photon counting search, metres\n\
-noiseMult x;     noise multiplier for photon\
-counting\n");
  #endif
  fprintf(stdout,"\nDenoising:\n\
-meanN n;         mean noise level, if using a predefined mean level\n\
-thresh n;        noise threshold, if using a predefined noise threshold\n\
-varNoise;        use a variable noise threshold\n\
-varScale x;      variable noise threshold scale (multiple of stdev above mean to set threshold)\n\
-statsLen len;    length to calculate noise stats over for varNoise\n\
-noiseTrack;      use noise tracking\n\
-sWidth sig;      smoothing width, after densoising\n\
-psWidth sigma;   smoothing width, before denoising\n\
-msWidth sig;     smoothing width, after noise stats, before denoising\n\
-preMatchF;       matched filter before denoising\n\
-postMatchF;      matched filter after denoising\n\
-pFile file;      read pulse file, for deconvoltuion and matched filters\n\
-gWidth sig;      Gaussian parameter selection smoothing width\n\
-minGsig sig;     minimum Gaussian sigma to fit\n\
-minWidth n;      minimum feature width in bins\n\
-medNoise;        use median stats rather than mean\n\
-varDrift;        correct detector drift with variable factor\n\
-driftFac xi;     fix drift with constant drift factor\n\
-rhoG rho;        ground reflectance\n\
-rhoC rho;        canopy reflectance\n\
-pSigma sig;      pulse width to smooth by if using Gaussian pulse\n\
-gold;            deconvolve with Gold's method\n\
-deconTol;        deconvolution tolerance\n\
\nQuestions to svenhancock@gmail.com\n\n");

  return;
}/*writeHelp*/

/*the end*/
/*###########################################################*/

