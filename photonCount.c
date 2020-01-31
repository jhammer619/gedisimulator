#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "hdf5.h"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "libOctree.h"
#include "tools.h"
#include "gediIO.h"
#include "gediNoise.h"
#include "photonCount.h"



/*##############################*/
/*# Generates photon count from#*/
/*# simulated GEDI waveforms   #*/
/*# 2019 svenhancock@gmail.com #*/
/*##############################*/

/*#######################################*/
/*# Copyright 2015-2019, Steven Hancock #*/
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


/*####################################################*/
/*turn a waveform in to a photon-count pseudo-wave*/

void uncompressPhotons(float *denoised,dataStruct *data,photonStruct *photonCount,denPar *den,noisePar *noise)
{
  int i=0,nPhot=0;
  int bin=0;
  float *temp=NULL;
  float **phots=NULL;

  /*allocate space*/
  temp=falloc(data->nBins,"temp pcl photon",0);

  /*extract photon coords along with their flags*/
  phots=countPhotons(denoised,data,photonCount,&nPhot,den,noise);

  /*bin up in to new wave*/
  for(i=0;i<nPhot;i++){
    bin=(int)(((float)data->z[0]-phots[0][i])/data->res);
    temp[bin]+=1.0;
  }
  TTIDY((void **)phots,3);

  /*copy over results*/
  TIDY(denoised);
  denoised=temp;
  temp=NULL;

  return;
}/*extractPhotons*/


/*####################################################*/
/*count photons*/

float **countPhotons(float *denoised,dataStruct *data,photonStruct *photonCount,int *nPhot,denPar *den,noisePar *noise)
{
  int i=0;
  int nPhotons=0,nNoise=0;
  int setNumberNoise(float,float,float);
  float **phots=NULL;  /*arrray with z, isSignal and isGround*/
  float photThresh=0,d=0,thisZ=0;
  float pickArrayElement(float,float *,int,char);
  float noiseInt=0;
  float photonNoiseIntensity(float);
  float minZ=0,maxZ=0;
  float *thisGr=NULL;
  float *wave=NULL;
  float *adjustPhotonProb(float *,dataStruct *,denPar *,noisePar *,int,photonStruct *);
  void setPhotonProb(photonStruct *);
  void setPhotonGround(float *,float *,float,double,float *,float *,double *,int);
  char testPhotonGround(dataStruct *,float);


  /*rescale waveform for reflectance*/
  wave=adjustPhotonProb(denoised,data,den,noise,data->useType,photonCount);

  /*do we need to set up the probability array?*/
  if(photonCount->prob==NULL)setPhotonProb(photonCount);

  /*choose a number of signal photons to use*/
  photThresh=(float)rand()/(float)RAND_MAX;
  nPhotons=(int)pickArrayElement(photThresh,photonCount->prob,photonCount->pBins,0);

 /*determine reflectance for noise intensity*/
  noiseInt=photonNoiseIntensity(data->cov);
  /*generate noise photons*/
  nNoise=setNumberNoise(data->cov,photonCount->noise_mult,photonCount->H);
  *nPhot=nPhotons+nNoise;

  /*allocate space*/
  phots=fFalloc(3,"photon coordinates",0);
  for(i=0;i<3;i++)phots[i]=falloc(*nPhot,"photon coordinates",i+1);

  /*generate signal photons*/
  for(i=0;i<nPhotons;i++){
    /*pick a point along the waveform*/
    photThresh=(float)rand()/(float)RAND_MAX;
    d=pickArrayElement(photThresh,wave,data->nBins,1);

    phots[0][i]=(float)data->z[0]-d*data->res;   /*determine range*/
    phots[1][i]=1.0;                             /*is signal*/
    phots[2][i]=(float)testPhotonGround(data,d); /*is this ground or canopy?*/
  }/*signal photon loop*/

  /*Noise*/
  /*set bounds of search window*/
  if(data->ground)thisGr=data->ground[data->useType];
  else            thisGr=NULL;
  setPhotonGround(&minZ,&maxZ,photonCount->H,data->gElev,data->wave[data->useType],thisGr,data->z,data->nBins);
  thisGr=NULL;

  /*add noise photons*/
  for(i=0;i<nNoise;i++){
    thisZ=(maxZ-minZ)*((float)rand()/(float)RAND_MAX)+minZ;

    phots[0][i+nPhotons]=thisZ;   /*determine range*/
    phots[1][i+nPhotons]=0.0;     /*is signal*/
    phots[2][i+nPhotons]=0.0;     /*is this ground or canopy?*/
  }/*noise loop*/

  /*tidy up*/
  if(wave!=denoised){
    TIDY(wave);
  }else wave=NULL;
  return(phots);
}/*countPhotons*/


/*####################################################*/
/*select photons for photon counting*/

void photonCountCloud(float *denoised,dataStruct *data,photonStruct *photonCount,char *outRoot,int numb,denPar *den,noisePar *noise)
{
  int i=0,nRH=0,nPhot=0;
  float *rhReal=NULL,noiseInt=0;
  float photonNoiseIntensity(float);
  float **phots=NULL;

  /*open file if needed*/
  if(photonCount->opoo==NULL){
    sprintf(photonCount->outNamen,"%s.pts",outRoot);
    if((photonCount->opoo=fopen(photonCount->outNamen,"w"))==NULL){
      fprintf(stderr,"Error opening input file %s\n",photonCount->outNamen);
      exit(1);
    }
    fprintf(photonCount->opoo,"# 1 X, 2 Y, 3 Z, 4 minht, 5 WFGroundZ, 6 RH50, 7 RH60, 8 RH75, 9 RH90, 10 RH95, 11 CanopyZ, 12 canopycover, 13 shot#, 14 photon#, 15 iteration#, 16 refdem, 17 noiseInt, 18 signal, 19 ground\n");
  }

  /*generate photons*/
  phots=countPhotons(denoised,data,photonCount,&nPhot,den,noise);

  /*get true RH metrics*/
  rhReal=findRH(data->wave[data->useType],data->z,data->nBins,data->gElev,5.0,&nRH);

  /*determine reflectance for noise intensity*/
  noiseInt=photonNoiseIntensity(data->cov);

  /*write out photons*/
  for(i=0;i<nPhot;i++){
    fprintf(photonCount->opoo,"%.2f %.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %d %d 1 %.3f %.3f %d %d\n",data->lon,data->lat,phots[0][i],rhReal[0],data->gElev,rhReal[10],rhReal[12],rhReal[15],rhReal[18],rhReal[19],rhReal[nRH-1],data->cov,numb,i,data->gElev,noiseInt,(int)phots[1][i],(int)phots[2][i]);
  }/*mutiple photon loop*/

  TTIDY((void **)phots,3);
  TIDY(rhReal);
  return;
}/*photonCountCloud*/


/*########################################################*/
/*is this photon from the ground?*/

char testPhotonGround(dataStruct *data,float d)
{
  int bin=0;
  float gFrac=0;
  char isGround=0;

  if(data->ground==NULL)isGround=0;
  else{
    bin=(int)d;
    if(data->wave[data->useType][bin]>0.0){
      gFrac=data->ground[data->useType][bin]/data->wave[data->useType][bin];
      isGround=(gFrac>=0.5)?1:0;
    }else isGround=0;
  }

  return(isGround);
}/*testPhotonGround*/


/*########################################################*/
/*adjust waveform to account for refl difference*/

float *adjustPhotonProb(float *denoised,dataStruct *data,denPar *den,noisePar *noise,int numb,photonStruct *phot)
{
  int i=0;
  float tot=0;
  float *wave=NULL,*canopy=NULL;
  float *smooGr=NULL,*smooCan=NULL;

  /*is any adjustment needed*/
  if(fabs(1.0-phot->rhoVrhoG)<TOL)wave=denoised;
  else{
    if(den->varNoise||noise->linkNoise){
      fprintf(stderr,"Not able to readjust denoised waveforms just yet\n");
      exit(1);
    }else{
      /*find canopy portion*/
      canopy=falloc((uint64_t)data->nBins,"canopy wave",0);
      for(i=0;i<data->nBins;i++)canopy[i]=data->wave[numb][i]-data->ground[numb][i];
      /*smooth ground if needed*/
      if((den->sWidth>0.001)||(den->psWidth>0.001)||(den->msWidth>0.001)){
        /*smooth if needed*/
        smooCan=processFloWave(canopy,data->nBins,den,1.0);
        smooGr=processFloWave(data->ground[numb],data->nBins,den,1.0);
      }else{
        smooCan=canopy;
        smooGr=data->ground[numb];
      }
      /*add up and normalise*/
      wave=falloc((uint64_t)data->nBins,"rescaled erflectance wave",0);
      tot=0.0;
      for(i=0;i<data->nBins;i++){
        wave[i]=smooCan[i]*phot->rhoVrhoG+smooGr[i]/phot->rhoVrhoG;
        tot+=wave[i];
      }
      if(fabs(1.0-tot)>TOL){
        for(i=0;i<data->nBins;i++)wave[i]/=tot;
      }
    }
  }

  /*tidy up*/
  if(smooCan!=canopy){
    TIDY(smooCan);
  }
  TIDY(canopy);
  if(smooGr!=data->ground[numb]){
    TIDY(smooGr);
  }
  return(wave);
}/*adjustPhotonProb*/


/*########################################################*/
/*determine ground for photon counting noise*/

void setPhotonGround(float *minZ,float *maxZ,float H,double gElev,float *wave,float *ground,double *z,int nBins)
{
  int i=0;
  float tot=0;
  float CofG=0;

  /*do we have any ground energy?*/
  tot=0.0;
  if(ground){
    for(i=0;i<nBins;i++)tot+=ground[i];
  }

  /*if no ground return, use CofG*/
  if(tot<=TOL){
    tot=0.0;
    for(i=0;i<nBins;i++){
      tot+=wave[i];
      CofG+=(float)z[i]*wave[i];
    }
    CofG/=tot;
    *maxZ=(float)gElev+H/2.0;
    *minZ=(float)gElev-H/2.0;
  }else CofG=gElev;  /*otehrwise use the ground elevation*/
  *maxZ=(float)CofG+H/2.0;
  *minZ=(float)CofG-H/2.0;

  return;
}/*setPhotonGround*/


/*########################################################*/
/*set number of noise photons for photon-counting*/

int setNumberNoise(float cov,float noise_mult,float H)
{
  int nNoise=0;
  float refl=0;
  float noiseRate=0;
  float c=299792458.0;

  /*surface reflectance*/
  if((cov<0.0)||(cov>1.0))cov=0.5;
  refl=cov*0.15+(1.0-cov)*0.22;  /*assuming ground and canopy reflectance values*/
  noiseRate=noise_mult*refl*pow(10.0,6.0);
  nNoise=(int)(50.0*(H/c)*noiseRate+0.5);

  return(nNoise);
}/*setNumberNoise*/


/*########################################################*/
/*determine solar background intensity for photon count*/

float photonNoiseIntensity(float cov)
{
  float noiseInt=0;

  /*it should be a backwards average of this, using the FOV (1D slice thereof)*/
  /*noiseInt=1.0+(1.0-cov)*1.5;*/

  /*this is what is in the ICEsat-2 code, just about*/
  if(cov<0.25)noiseInt=2.5;
  else if(cov>0.75)noiseInt=1.0;
  else noiseInt=1.5+(float)rand()/(float)RAND_MAX;

  return(noiseInt);
}/*photonNoiseIntensity*/


/*####################################################*/
/*set photon probability*/

void setPhotonProb(photonStruct *photonCount)
{
  int i=0;
  float y=0;
  float poissonPDF(float,float);

  /*determine number of steps*/
  photonCount->pBins=0;
  do{
    y=poissonPDF((float)photonCount->pBins,photonCount->designval);
    photonCount->pBins++;
  }while((y>0.00001)||((float)photonCount->pBins<photonCount->designval));  /*at least to mean and then to low prob*/

  /*allocate space*/
  photonCount->prob=falloc((uint64_t)photonCount->pBins,"photon prob",0);

  /*set probabilities*/
  for(i=0;i<photonCount->pBins;i++)photonCount->prob[i]=poissonPDF((float)i,photonCount->designval);

  return;
}/*setPhotonProb*/


/*####################################################*/
/*Amy's Poission function*/

float poissonPDF(float n,float lambda)
{
  float y=0;
  y=exp(-1.0*lambda+n*log(lambda)-lgamma(n+1.0));
  return(y);
}/*poissonPDF*/


/*####################################################*/
/*determine point at which array exceeds threshold*/

float pickArrayElement(float photThresh,float *jimlad,int nBins,char interpolate)
{
  int i=0;
  float x=0,y0=0;
  float tot=0,*cumul=NULL;

  /*determine total energy and adjust threshold*/
  tot=0.0;
  cumul=falloc((uint64_t)nBins,"cumul",0);
  for(i=0;i<nBins;i++){
    tot+=jimlad[i];
    if(i>0)cumul[i]=cumul[i-1]+jimlad[i];
    else   cumul[i]=jimlad[i];
  }
  photThresh*=tot;

  /*determine point above*/
  for(i=0;i<nBins;i++)if(cumul[i]>=photThresh)break;

  /*extrapolate between two elements*/
  if(interpolate){
    if(i>0)y0=cumul[i-1];
    else   y0=0.0;
    if(i<(nBins-1))x=(cumul[i]-photThresh)/(cumul[i]-y0)+(float)i;
    else           x=(float)(nBins-1);
  }else x=(float)i;
  TIDY(cumul);


  return(x);
}/*pickArrayElement*/


/*########################################################*/
/*set photon rates*/

void setPhotonRates(photonStruct *photonCount)
{
  /*do we need to?*/
  if(photonCount->nPhotG>0.0){
    photonCount->designval=(photonCount->nPhotC+photonCount->nPhotG)/2.0;
    photonCount->rhoVrhoG=photonCount->nPhotC/photonCount->nPhotG;
  }
  return;
}/*setPhotonRates*/


/*the end*/
/*########################################################*/

