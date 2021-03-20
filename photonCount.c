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
#include "functionWrappers.h"
#include "gediIO.h"
#include "gediNoise.h"
#include "photonCount.h"
#include "gsl/gsl_fft_complex.h"

//#define DEBUG


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

float *uncompressPhotons(float *wave,dataStruct *data,photonStruct *photonCount,noisePar *noise,gediIOstruct *gediIO)
{
  float *photWave=NULL;
  float *corrWave=NULL;
  float *crossCorrelateWaves(float *,float,int,pulseStruct *,float);
  float *crossCorrelateTime(float *,float,int,pulseStruct *,float);

  #ifdef DEBUG
  int i=0;
  static int c=0;
  #endif

  /*do we have a usable pulse?*/
  if(gediIO->pulse==NULL){
    errorf("No pulse. Cannot use PCL\n");
    return(NULL);
  }

  /*first perform photon counting, if needed*/
  if(gediIO->pclPhoton) { ASSIGN_CHECKNULL_RETNULL(photWave,countWaveform(wave,data,photonCount,gediIO->den,noise)); }
  else                 photWave=wave;

  /*perform cross-correlation*/
  ASSIGN_CHECKNULL_RETNULL(corrWave,crossCorrelateTime(photWave,data->res,data->nBins,gediIO->pulse,gediIO->pRes));
  //corrWave=crossCorrelateWaves(photWave,data->res,data->nBins,gediIO->pulse,gediIO->pRes);

  #ifdef DEBUG
  for(i=0;i<data->nBins;i++)msgf("%d %f %f %f %f debug2\n",c,data->z[i],wave[i],photWave[i],corrWave[i]);
  c++;
  #endif


  /*tidy up*/
  if(photWave!=wave){
    TIDY(photWave);
  }else photWave=NULL;

  return(corrWave);
}/*uncompressPhotons*/


/*####################################################*/
/*resample a pulse for PCL*/

int resamplePclPulse(pulseStruct *pulse,float res,float pRes)
{
  int i=0,*contN=NULL;
  int bin=0;

  /*allocate space and zero*/
  pulse->rBins=(int)((float)pulse->nBins*pRes/res);
  ASSIGN_CHECKNULL_RETINT(pulse->resamp,falloc(pulse->rBins,"",0));
  ASSIGN_CHECKNULL_RETINT(contN,ialloc(pulse->rBins,"",0));
  for(i=0;i<pulse->rBins;i++){
    pulse->resamp[i]=0.0;
    contN[i]=0;
  }


  /*resample pulse*/
  for(i=0;i<pulse->nBins;i++){
    bin=(int)((float)i*pRes/res);
    if((bin>=0)&&(bin<pulse->rBins)){
      pulse->resamp[bin]+=pulse->y[i];
      contN[bin]++;
    }
  }

  /*normalise resampled*/
  for(i=0;i<pulse->rBins;i++){
    if(contN[i]>0)pulse->resamp[i]/=(float)contN[i];
  }
  TIDY(contN);
  pulse->rCent=(int)((float)pulse->centBin*pRes/res);

  return 0;
}/*resamplePclPulse*/


/*####################################################*/
/*perform a cross-correlation in the time domain*/

float *crossCorrelateTime(float *photWave,float res,int nBins,pulseStruct *pulse,float pRes)
{
  int i=0,j=0,bin=0;
  int thisCont=0;
  float *compCorr=NULL;
  float meanP=0,meanW=0;
  float stdevP=0,stdevW=0;
  int resamplePclPulse(pulseStruct *,float,float);

  /*allocate space*/
  ASSIGN_CHECKNULL_RETNULL(compCorr,falloc(nBins,"compCorr",0));

  /*allocate resampled pulse if needed*/
  if(pulse->resamp==NULL) {ISINTRETNULL(resamplePclPulse(pulse,res,pRes)); }

  /*find the average of the pulse*/
  meanP=0.0;
  for(i=0;i<pulse->rBins;i++)meanP+=pulse->resamp[i];
  meanP/=(float)pulse->rBins;
  //meanP=singleMedian(pulse->resamp,pulse->rBins);

  /*find the stdev of the pulse*/
  stdevP=0.0;
  for(i=0;i<pulse->rBins;i++)stdevP+=(pulse->resamp[i]-meanP)*(pulse->resamp[i]-meanP);
  stdevP=sqrt(stdevP/(float)pulse->rBins);

  /*find the average of the wave*/
  meanW=0.0;
  for(i=0;i<nBins;i++)meanW+=photWave[i];
  meanW/=(float)nBins;
  //meanW=singleMedian(photWave,nBins);

  /*find the stdev of the wave*/
  stdevW=0.0;
  for(i=0;i<nBins;i++)stdevW+=(photWave[i]-meanW)*(photWave[i]-meanW);
  stdevW=sqrt(stdevW/(float)nBins);

  /*loop over bins in time domain*/
  //for(i=pulse->rCent;i<(nBins-pulse->rCent);i++){ /*step in by the cent bins of the pulse*/
  for(i=0;i<nBins;i++){
    compCorr[i]=0.0;
    thisCont=0;

    /*loop over pulse to convolve*/
    for(j=0;j<pulse->rBins;j++){
      bin=i+j-pulse->rCent;  /*bin on the pulse*/

      /*are we within the pulse array?*/
      if((bin>=0)&&(bin<nBins)){
        compCorr[i]+=(photWave[bin]-meanW)*(pulse->resamp[pulse->rBins-(j+1)]-meanP)/(stdevP*stdevW);
        //compCorr[i]+=photWave[bin]*pulse->resamp[pulse->rBins-(j+1)]/(stdevP*stdevW);
        thisCont++;
      }
    }/*pulse bin loop*/

    /*normalise*/
    //if(thisCont>0)compCorr[i]/=(float)thisCont;
    compCorr[i]/=(float)pulse->rBins;
  }/*wave bin loop*/

  return(compCorr);
}/*crossCorrelateTime*/


/*####################################################*/
/*perform a cross-correlation using Fourier*/

float *crossCorrelateWaves(float *photWave,float res,int nBins,pulseStruct *pulse,float pRes)
{
  int i=0,bin=0;
  int numb=0;
  int *contN=NULL;
  float *corrWave=NULL;
  double *compPulse=NULL;
  double *compWave=NULL;
  double *compCorr=NULL;
  int gsl_fft_complex_radix2_forward(gsl_complex_packed_array,size_t,size_t);
  int gsl_fft_complex_radix2_backward(gsl_complex_packed_array, size_t,size_t);

  /*FFT requires that array is a power of 2 long*/
  numb=pow(2.0,(float)(int)(log((double)nBins)/log(2.0)+1.0));

  /*allocate space for resampled complex pulse*/
  ASSIGN_CHECKNULL_RETNULL(compPulse,dalloc(2*numb,"complex pulse",0));
  ASSIGN_CHECKNULL_RETNULL(contN,ialloc(numb,"contribution counter",0));
  for(i=0;i<numb;i++){
    compPulse[2*i]=compPulse[2*i+1]=0.0;
    contN[i]=0;
  }

  /*split pulse over 0*/
  /*for(i=pulse->centBin;i<pulse->nBins;i++){
    bin=(int)((float)(i-pulse->centBin)*pRes/res+0.5);
    if((bin<0)||(bin>=numb)){
      errorf("Out of bounds when splitting pulse, %d of %d\n",bin,numb);
      continue;
    }
    compPulse[2*bin]+=(double)pulse->y[i];
    contN[bin]++;
  }
  for(i=0;i<pulse->centBin;i++){
    bin=numb-((int)((float)(pulse->centBin-i)*pRes/res+0.5)+1);
    if((bin<0)||(bin>=numb)){
      errorf("Out of bounds when splitting pulse, %d of %d\n",bin,numb);
      continue;
    }
    compPulse[2*bin]+=(double)pulse->y[i];
    contN[bin]++;
  }*/


  /*resample pulse to match waveform*/
  for(i=0;i<pulse->nBins;i++){
    bin=(int)((float)i*pRes/res+0.5);
    if((bin<0)||(bin>=numb))continue;
    compPulse[2*bin]+=(double)pulse->y[i];
    contN[bin]++;
  }

  /*normalise*/
  for(i=0;i<numb;i++){
    if(contN[i]>0)compPulse[2*i]/=(float)contN[i];
  }
  TIDY(contN);

  /*median of pulse*/
  //meanP=singleMedian(pulse->y,pulse->nBins);

  /*make waveform complex*/
  ASSIGN_CHECKNULL_RETNULL(compWave,dalloc(2*numb,"complex wave",0));
  for(i=0;i<nBins;i++){
    compWave[2*i]=(double)photWave[i];
    compWave[2*i+1]=0.0;  /*imaginary part*/
  }
  //meanW=singleMedian(photWave,nBins);
  for(i=2*nBins;i<2*numb;i++)compWave[i]=0.0;

  /*subtract means*/
  /*for(i=numb-1;i>=0;i--){
    compWave[2*i]-=meanW;
    compPulse[2*i]-=meanP;
  }*/

  /*remove assymmetry of signal*/
  //removeAsymmetryPCL(compWave,numb);

  /*fourier transform both*/
  gsl_fft_complex_radix2_forward((gsl_complex_packed_array)compPulse,1,numb);
  gsl_fft_complex_radix2_forward((gsl_complex_packed_array)compWave,1,numb);

  /*correlate*/
  ASSIGN_CHECKNULL_RETNULL(compCorr,dalloc(2*numb,"complex correlation",0));
  for(i=0;i<numb;i++){
    compCorr[2*i]=compPulse[2*i]*compWave[2*i]+compPulse[2*i+1]*compWave[2*i+1];
    compCorr[2*i+1]=compPulse[2*i+1]*compWave[2*i]-compPulse[2*i]*compWave[2*i+1];
  }

  /*inverse fourier*/
  gsl_fft_complex_radix2_backward((gsl_complex_packed_array)compCorr,1,numb);

  /*make real*/
  ASSIGN_CHECKNULL_RETNULL(corrWave,falloc(nBins,"correlated wave",0));
  for(i=0;i<nBins;i++){
    corrWave[i]=(float)compCorr[2*i]; //sqrt(compCorr[2*i]*compCorr[2*i]+compCorr[2*i+1]*compCorr[2*i+1]);
  }

  /*tidy up*/
  TIDY(compPulse);
  TIDY(compWave);
  TIDY(compCorr);

  return(corrWave);
}/*crossCorrelateWaves*/


/*####################################################*/
/*truncate the assymmetry of a waveform for PCL*/

void removeAsymmetryPCL(float *wave,int numb)
{
  int i=0;
  float medianW=0;

  /*find the median*/
  medianW=singleMedian(wave,numb);

  /*find start point*/
  for(i=0;i<numb;i++){
    /*if less than mean, set to mean*/
    if(wave[i]>=medianW){
      for(;i>=0;i--)wave[i]=medianW;
      break;
    }
  }

  /*find end point*/
  for(i=numb-1;i>=0;i--){
    /*if less then mean, set to mean*/
    if(wave[i]>=medianW){
      for(;i<numb;i++)wave[i]=medianW;
      break;
    }
  }

  return;
}/*removeAsymmetryPCL*/


/*####################################################*/
/*produce a photon-counting pseudo-waveform*/

float *countWaveform(float *denoised,dataStruct *data,photonStruct *photonCount,denPar *den,noisePar *noise)
{
  int i=0,nPhot=0;
  int bin=0;
  float minI=0;
  float shotSig=0;
  float shotNoise=0;
  float *temp=NULL;
  float **phots=NULL;

  #ifdef DEBUG
  static int count=0;
  msgf("Photon counting\n");
  #endif

  /*set minimum to zero*/
  minI=10000.0;
  for(i=0;i<data->nBins;i++){
    if(denoised[i]<minI)minI=denoised[i];
  }
  ASSIGN_CHECKNULL_RETNULL(temp,falloc(data->nBins,"temp pcl waveform",0));
  for(i=0;i<data->nBins;i++)temp[i]=denoised[i]-minI;

  /*set window size for background noise photons*/
  photonCount->H=data->res*(float)data->nBins*2.0;  /*this is the two way distance*/

  /*extract photon coords along with their flags*/
  ASSIGN_CHECKNULL_RETNULL(phots,countPhotons(temp,data,photonCount,&nPhot,den,noise));

  /*reset temp array*/
  for(i=0;i<data->nBins;i++)temp[i]=0.0;

  /*bin up in to new wave*/
  for(i=0;i<nPhot;i++){
    bin=(int)(((float)data->z[0]-phots[0][i])/data->res);
    if((bin>=0)&&(bin<data->nBins))temp[bin]+=1.0;
  }
  TTIDY((void **)phots,3);

  /*apply shot noise if neeed*/
  if(noise->shotNoise){
    /*loop over bins*/
    for(i=0;i<data->nBins;i++){
      /*only if there are photons*/
      if(temp[i]>0.0){
        /*set sigma*/
        shotSig=sqrt(temp[i]);
  
        /*draw Gaussian random number*/
        shotNoise=GaussNoise()*shotSig;

        /*add noise and truncate negative*/
        temp[i]+=round(shotNoise);
        if(temp[i]<0.0)temp[i]=0.0;
      }/*return check*/
    }/*bin loop*/
  }/*apply shot noise check*/

  #ifdef DEBUG
  for(i=0;i<data->nBins;i++)msgf("%d %f %f %f ta\n",count,data->z[i],temp[i],data->wave[0][i]);
  count++;
  #endif

  return(temp);
}/*countWaveform*/


/*####################################################*/
/*produce photon counting photons*/

float **countPhotons(float *denoised,dataStruct *data,photonStruct *photonCount,int *nPhot,denPar *den,noisePar *noise)
{
  int i=0;
  int nPhotons=0,nNoise=0;
  int setNumberNoise(float,float,float);
  float **phots=NULL;  /*arrray with z, isSignal and isGround*/
  float photThresh=0,d=0,thisZ=0;
  float photonNoiseIntensity(float);
  float minZ=0,maxZ=0;
  float *thisGr=NULL;
  float *wave=NULL;
  float *adjustPhotonProb(float *,dataStruct *,denPar *,noisePar *,int,photonStruct *);
  void knockOffNegativeWaves(float *,dataStruct *);
  void adjustTotalPhotRate(photonStruct *,float);

  void setPhotonGround(float *,float *,float,double,float *,float *,double *,int);
  char testPhotonGround(dataStruct *,float);


  /*remove negagives if needed*/
  knockOffNegativeWaves(denoised,data);

  /*rescale waveform for reflectance*/
  ASSIGN_CHECKNULL_RETNULL(wave,adjustPhotonProb(denoised,data,den,noise,data->useType,photonCount));

  /*do we need to set up the probability array?*/
  if(photonCount->prob==NULL){
    if(photonCount->reflDiff)adjustTotalPhotRate(photonCount,data->cov);
    ISINTRETNULL(setPhotonProb(photonCount));
  }

  /*choose a number of signal photons to use*/
  photThresh=frand();
  ASSIGN_CHECKINT_RETNULL(nPhotons,(int)pickArrayElement(photThresh,photonCount->prob,photonCount->pBins,0));

  /*generate noise photons*/
  nNoise=setNumberNoise(data->cov,photonCount->noise_mult,photonCount->H);
  *nPhot=nPhotons+nNoise;

  /*allocate space*/
  ASSIGN_CHECKNULL_RETNULL(phots,fFalloc(3,"photon coordinates",0));
  for(i=0;i<3;i++) {ASSIGN_CHECKNULL_RETNULL(phots[i],falloc(*nPhot,"photon coordinates",i+1));}

  /*generate signal photons*/
  for(i=0;i<nPhotons;i++){
    /*pick a point along the waveform*/
    photThresh=frand();
    d=pickArrayElement(photThresh,wave,data->nBins,1);
    if(fabs(d+1.0)< TOL)
      return(NULL);

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
  #ifdef DEBUG
  fprintf(stdout,"Adding %d noise over %f signal %d\n",nNoise,photonCount->H,nPhotons);
  #endif
  for(i=0;i<nNoise;i++){
    thisZ=(maxZ-minZ)*frand()+minZ;

    phots[0][i+nPhotons]=thisZ;   /*determine range*/
    phots[1][i+nPhotons]=0.0;     /*is signal*/
    phots[2][i+nPhotons]=0.0;     /*is this ground or canopy?*/
  }/*noise loop*/

  /*clear out prob if needed*/
  if(photonCount->reflDiff)TIDY(photonCount->prob);

  /*tidy up*/
  if(wave!=denoised){
    TIDY(wave);
  }else wave=NULL;

  return(phots);
}/*countPhotons*/


/*####################################################*/
/*select photons for photon counting*/

int photonCountCloud(float *denoised,dataStruct *data,photonStruct *photonCount,char *outRoot,int numb,denPar *den,noisePar *noise)
{
  int i=0,nRH=0,nPhot=0;
  float *rhReal=NULL,noiseInt=0;
  float photonNoiseIntensity(float);
  float **phots=NULL;

  /*open file if needed*/
  if(photonCount->opoo==NULL){
    sprintf(photonCount->outNamen,"%s.pts",outRoot);
    if((photonCount->opoo=fopen(photonCount->outNamen,"w"))==NULL){
      errorf("Error opening input file %s\n",photonCount->outNamen);
      return(-1);
    }
    fprintf(photonCount->opoo,"# 1 X, 2 Y, 3 Z, 4 minht, 5 WFGroundZ, 6 RH50, 7 RH60, 8 RH75, 9 RH90, 10 RH95, 11 CanopyZ, 12 canopycover, 13 shot#, 14 photon#, 15 iteration#, 16 refdem, 17 noiseInt, 18 signal, 19 ground\n");
  }

  /*generate photons*/
  phots=countPhotons(denoised,data,photonCount,&nPhot,den,noise);

  /*get true RH metrics*/
  ASSIGN_CHECKNULL_RETINT(rhReal,findRH(data->wave[data->useType],data->z,data->nBins,data->gElev,5.0,&nRH));

  /*determine reflectance for noise intensity*/
  noiseInt=photonNoiseIntensity(data->cov);

  /*write out photons*/
  for(i=0;i<nPhot;i++){
    fprintf(photonCount->opoo,"%.2f %.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %d %d 1 %.3f %.3f %d %d\n",data->lon,data->lat,phots[0][i],rhReal[0],data->gElev,rhReal[10],rhReal[12],rhReal[15],rhReal[18],rhReal[19],rhReal[nRH-1],data->cov,numb,i,data->gElev,noiseInt,(int)phots[1][i],(int)phots[2][i]);
  }/*mutiple photon loop*/

  TTIDY((void **)phots,3);
  TIDY(rhReal);
  return(0);
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
/*remove negative values for photon counting*/

void knockOffNegativeWaves(float *denoised,dataStruct *data)
{
  int i=0;
  float min=0;

  /*find the minimum*/
  min=100000.0;
  for(i=0;i<data->nBins;i++){
    if(denoised[i]<min)min=denoised[i];
  }


  /*translate waves if needed*/
  if(min<0.0){
    for(i=0;i<data->nBins;i++){
      denoised[i]-=min;
      if(data->ground)data->ground[data->useType][i]-=min;
    }
  }

  return;
}/*knockOffNegativeWaves*/


/*########################################################*/
/*adjust waveform to account for refl difference*/

float *adjustPhotonProb(float *denoised,dataStruct *data,denPar *den,noisePar *noise,int numb,photonStruct *phot)
{
  int i=0;
  float tot=0;
  float *wave=NULL,*canopy=NULL;
  float *smooGr=NULL,*smooCan=NULL;

  /*do we have a ground*/
  if(data->ground==NULL){  /*no ground*/
    wave=denoised;
  }else{   /*there is a ground*/
    /*is any adjustment needed*/
    if(fabs(1.0-phot->rhoVrhoG)<TOL)wave=denoised;
    else{  /*if it is needed, do it*/
      if(den->varNoise||noise->linkNoise){
        errorf("Not able to readjust denoised waveforms just yet\n");
        return NULL;
      }else{
        /*find canopy portion*/
        ASSIGN_CHECKNULL_RETNULL(canopy,falloc((uint64_t)data->nBins,"canopy wave",0));
        for(i=0;i<data->nBins;i++)canopy[i]=data->wave[numb][i]-data->ground[numb][i];

        /*smooth ground if needed*/
        if((den->sWidth>0.001)||(den->psWidth>0.001)||(den->msWidth>0.001)){
          /*smooth if needed*/
          ASSIGN_CHECKNULL_RETNULL(smooCan,processFloWave(canopy,data->nBins,den,1.0));
          ASSIGN_CHECKNULL_RETNULL(smooGr,processFloWave(data->ground[numb],data->nBins,den,1.0));
        }else{
          smooCan=canopy;
          smooGr=data->ground[numb];
        }

        /*add up and normalise*/
        ASSIGN_CHECKNULL_RETNULL(wave,falloc((uint64_t)data->nBins,"rescaled erflectance wave",0));
        tot=0.0;
        for(i=0;i<data->nBins;i++){
          wave[i]=smooCan[i]*phot->nPhotC+smooGr[i]*phot->nPhotG;
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
  }/*is there a ground check?*/

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
  /*float refl=0; OLD, for varying reflectance*/
  float photThresh=0;
  float noiseRate=0;
  float c=299792458.0;
  photonStruct tempPhot;

  /*is any noise being added?*/
  if(noise_mult>TOL){
    /*surface reflectance*/
    /*if((cov<0.0)||(cov>1.0))cov=0.5;   OLD, for varying ground reflectance
    refl=cov*0.15+(1.0-cov)*0.22;*/  /*assuming ground and canopy reflectance values*/

    /*noise rate in photons per window*/
    noiseRate=noise_mult*pow(10.0,6.0); //*refl     used to be *refl to account for xchanging surface reflectance
    /*nNoise=(int)(50.0*(H/c)*noiseRate+0.5);  This is to match Kaitlin's matlab code, but unsure where the 50 came from*/
    tempPhot.designval=(H/c)*noiseRate;

    /*pick from a Poisson*/
    ISINTRETINT(setPhotonProb(&tempPhot));
    photThresh=(float)frand()/(float)RAND_MAX;
    nNoise=(int)pickArrayElement(photThresh,tempPhot.prob,tempPhot.pBins,0);
    TIDY(tempPhot.prob);
  }else nNoise=0;

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
  else noiseInt=1.5+frand();

  return(noiseInt);
}/*photonNoiseIntensity*/


/*####################################################*/
/*adjust photon rate for changing reflectance*/

void adjustTotalPhotRate(photonStruct *photonCount,float cov)
{
  if(cov>=0.0)photonCount->designval=cov*photonCount->nPhotC+(1.0-cov)*photonCount->nPhotG;

  return;
}/*adjustTotalPhotRate*/


/*####################################################*/
/*set photon probability*/

int setPhotonProb(photonStruct *photonCount)
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
  ASSIGN_CHECKNULL_RETINT(photonCount->prob,falloc((uint64_t)photonCount->pBins,"photon prob",0));

  /*set probabilities*/
  for(i=0;i<photonCount->pBins;i++)photonCount->prob[i]=poissonPDF((float)i,photonCount->designval);

  return(0);
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
  if(cumul==NULL)
    return(-1.1);
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

  /*what mode has been used to define photon rates?*/
  if(fabs(photonCount->rhoVrhoG-1.0)>TOL){                   /*adjusted rhoV/rhoG*/
    photonCount->nPhotG=2.0*photonCount->designval/photonCount->rhoVrhoG-1.0;
    photonCount->nPhotC=photonCount->nPhotG*photonCount->rhoVrhoG;
    photonCount->reflDiff=1;
    /*prevent negative rates. only needed for black soil etc.*/
    if(photonCount->nPhotC<0.0)photonCount->nPhotC=0.0;
    if(photonCount->nPhotG<0.0)photonCount->nPhotG=0.0;

  }else if((photonCount->nPhotG+photonCount->nPhotC)>0.0){   /*separately defined photon rates*/
    photonCount->designval=(photonCount->nPhotC+photonCount->nPhotG)/2.0;
    photonCount->rhoVrhoG=photonCount->nPhotC/photonCount->nPhotG;
    photonCount->reflDiff=1;

  }else if(fabs(photonCount->designval)<TOL){                /*no reflectance correction. seperate rates*/
    photonCount->designval=(photonCount->nPhotC+photonCount->nPhotG)/2.0;
    photonCount->rhoVrhoG=1.0;
    photonCount->reflDiff=0;

  }else{                                                     /*no reflectance correction. total rate*/
    photonCount->nPhotG=photonCount->nPhotC=photonCount->designval;
    photonCount->rhoVrhoG=1.0;
    photonCount->reflDiff=0;

  }

  return;
}/*setPhotonRates*/


/*the end*/
/*########################################################*/

