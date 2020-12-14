#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "stdint.h"
#include "hdf5.h"
#include "tools.h"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "libOctree.h"
#include "gediIO.h"
#include "gediNoise.h"


/*##############################*/
/*# Adds noise to simulated    #*/
/*# GEDI waveforms             #*/
/*# 2018 svenhancock@gmail.com #*/
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


#define DRIFTTOL 0.000001

/*####################################################*/
/*add noise to waveform*/

void addNoise(dataStruct *data,noisePar *gNoise,float fSigma,float pSigma,float res,float rhoc,float rhog)
{
  int i=0;
  float noise=0;
  float tot=0.0,thresh=0;
  float minE=0;
  float *tempNoise=NULL;
  float GaussNoise();
  float *smooNoise=NULL,*tempWave=NULL;
  float *digitiseWave(float *,int,char,float,float);
  float reflScale=0;
  void deleteGround(float *,float *,float *,int,float,float,float,float,float,float,float);
  void scaleNoiseDN(float *,int,float,float,float);
  float *detectorDrift(float *,int,float,float);

  /*allocate*/
  data->noised=falloc((uint64_t)data->nBins,"noised wave",0);

  /*apply detecor drift if using*/
  tempWave=detectorDrift(data->wave[data->useType],data->nBins,gNoise->driftFact,res);

  if(gNoise->missGround){        /*Delete all signal beneath ground peak*/
    if(gNoise->minGap==0.0){
      fprintf(stderr,"Cannot delete ground without a defined minimum gap\n");
      exit(1);
    }
    deleteGround(data->noised,tempWave,data->ground[data->useType],data->nBins,gNoise->minGap,pSigma,fSigma,res,data->cov,rhoc,rhog);
  }else if(gNoise->linkNoise){   /*link margin based noise*/
    /*Gaussian noise*/
    tempNoise=falloc((uint64_t)data->nBins,"temp noised",0);
    /*in case of PCL, subtract min*/
    minE=1000000.0;
    for(i=0;i<data->nBins;i++)if(data->wave[data->useType][i]<minE)minE=data->wave[data->useType][i];
    if(minE>0.0)minE=0.0;
    tot=0.0;
    for(i=0;i<data->nBins;i++)tot+=(data->wave[data->useType][i]-minE)*res;
    reflScale=(data->cov*rhoc+(1.0-data->cov)*rhog)*tot/(gNoise->linkCov*rhoc+(1.0-gNoise->linkCov)*rhog);
    for(i=0;i<data->nBins;i++)tempNoise[i]=gNoise->linkSig*GaussNoise()*reflScale;
    /*smooth noise by detector response*/
    smooNoise=smooth(gNoise->deSig,data->nBins,tempNoise,res);
    for(i=0;i<data->nBins;i++)tempNoise[i]=tempWave[i]+smooNoise[i];
    TIDY(smooNoise);
    /*scale to match sigma*/
    scaleNoiseDN(tempNoise,data->nBins,gNoise->linkSig*reflScale,gNoise->trueSig,gNoise->offset);
    /*digitise*/
    TIDY(data->noised);
    data->noised=digitiseWave(tempNoise,data->nBins,gNoise->bitRate,gNoise->maxDN,tot);
    TIDY(tempNoise);
  }else if((gNoise->nSig>0.0)||(gNoise->meanN>0.0)){   /*mean and stdev based noise*/
    for(i=0;i<data->nBins;i++){
      noise=gNoise->nSig*GaussNoise();
      if((float)rand()/(float)RAND_MAX<0.5)noise*=-1.0; /*to allow negative numbers*/
      data->noised[i]=tempWave[i]+gNoise->meanN+noise;
    }/*bin loop*/
  }else if(gNoise->hNoise>0.0){  /*hard threshold noise*/
    tot=0.0;
    for(i=0;i<data->nBins;i++)tot+=data->wave[data->useType][i];
    thresh=gNoise->hNoise*tot;
    for(i=0;i<data->nBins;i++){
      data->noised[i]=tempWave[i]-thresh;
      if(data->noised[i]<0.0)data->noised[i]=0.0;
    }
  }else{  /*no noise*/
    for(i=0;i<data->nBins;i++)data->noised[i]=tempWave[i];
  }

  TIDY(tempWave);
  return;
}/*addNoise*/

/*####################################################*/
/*calculate sigma for link noise*/

float setNoiseSigma(float linkM,float cov,float pSigma,float fSigma,float rhoc,float rhog)
{
  float sig=0;
  float groundAmp=0;
  float slope=0;
  float gRefl=0;
  float tanSlope=0;
  float sigEff=0;          /*effective ground return width*/
  float probNoise=0,probMiss=0;
  float findSigma(float,float,float,float);

  slope=2.0*M_PI/180.0;

  gRefl=(1.0-cov)*rhoc;

  tanSlope=sin(slope)/cos(slope);
  sigEff=sqrt(pSigma*pSigma+fSigma*fSigma*tanSlope*tanSlope);
  groundAmp=(gRefl/(gRefl+rhoc*cov))/(sigEff*sqrt(2.0*M_PI));  /*normalise by total waveform reflectance*/

  probNoise=0.05;
  probMiss=0.1;
  sig=findSigma(probNoise,probMiss,groundAmp,linkM);

  return(sig);
}/*setNoiseSigma*/


/*####################################################*/
/*find sigma for this combination, from Xioali's notes*/

float findSigma(float probNoise,float probMiss,float groundAmp,float linkM)
{
  float sig=0;
  float step=0;
  float err=0,minErr=0;
  float thisLink=0;
  double nNsig=0,nSsig=0;
  float threshS=0,threshN=0;
  void gaussThresholds(double,double,double,double,double *,double *,noiseThreshStruct *);
  char direction=0;
  noiseThreshStruct noiseSigs;

  /*make a blank noise threshold structure*/
  noiseSigs.nThreshes=0;
  noiseSigs.threshN=noiseSigs.threshS=noiseSigs.probNoise=noiseSigs.probMiss=NULL;

  /*initial guess*/
  sig=groundAmp/2.0;
  step=groundAmp/10.0;

  /*determine threshold in terms of standard deviations*/
  gaussThresholds(1.0,XRES,(double)probNoise,(double)probMiss,&nNsig,&nSsig,&noiseSigs);

  direction=0;
  minErr=0.00015;
  do{
    /*scale thresholds by sigma*/
    threshN=(float)nNsig*sig;
    threshS=(float)nSsig*sig;

    thisLink=10.0*log((groundAmp-threshS)/threshN)/log(10.0);
    if(thisLink<linkM){
      if(direction==1)step/=2.0;
      sig-=step;
      direction=-1;
      if(sig<=0.0){
        step/=2.0;
        sig+=step;
        direction=0;
      }
    }else if(thisLink>linkM){
      if(direction==-1)step/=2.0;
      sig+=step;
      direction=1;
    }
    err=fabs(thisLink-linkM);
  }while(err>minErr);

  return(sig);
}/*findSigma*/


/*####################################################################################################*/
/*calculate Gaussian thresholds for normaly distributed noise*/

void gaussThresholds(double sig,double res,double probNoise,double probMiss,double *threshN,double *threshS,noiseThreshStruct *noiseSigs)
{
  int i=0;
  double x=0,y=0;
  double cumul=0,tot=0;
  double probN=0,probS=0;
  double *markDo(int,double *,double);
  char foundS=0,foundN=0;
  char beenCalc=0;         /*has this level already been calculated or not?*/

  /*do we need to calculate, or has it already been done?*/
  for(i=0;i<noiseSigs->nThreshes;i++){
    if((fabs(probNoise-noiseSigs->probNoise[i])<DRIFTTOL)&&(fabs(probMiss-noiseSigs->probMiss[i])<DRIFTTOL)){
      beenCalc=1;
      *threshN=noiseSigs->threshN[i];
      *threshS=noiseSigs->threshS[i];
      break;
    }
  }

  /*does this need calculating*/
  if(beenCalc==0){
    probN=1.0-probNoise/(30.0/0.15);
    probS=1.0-probMiss;

    /*determine start*/
    x=0.0;
    tot=0.0;
    do{
      y=gaussian(x,sig,0.0);
      if(fabs(x)>res)tot+=2.0*y;
      else           tot+=y;
      x-=res;
    }while(y>=YTOL);
    tot*=res;

    /*numerical integration and look for threshold crossings*/
    do{
      y=gaussian(x,sig,0.0);
      cumul+=y*res/tot;

      if(foundS==0){
        if(cumul>=probS){
          foundS=1;
          *threshS=x;
        }
      }
      if(foundN==0){
        if(cumul>=probN){
          foundN=1;
          *threshN=x;
        }
      }
      x+=res;
    }while((foundS==0)||(foundN==0));

    /*save the new thresholds*/
    noiseSigs->probNoise=markDo(noiseSigs->nThreshes,noiseSigs->probNoise,probNoise);
    noiseSigs->probMiss=markDo(noiseSigs->nThreshes,noiseSigs->probMiss,probMiss);
    noiseSigs->threshN=markDo(noiseSigs->nThreshes,noiseSigs->threshN,*threshN);
    noiseSigs->threshS=markDo(noiseSigs->nThreshes,noiseSigs->threshS,*threshS);
    noiseSigs->nThreshes++;
  }
  return;
}/*gaussThresholds*/


/*####################################################*/
/*digitse*/

float *digitiseWave(float *wave,int nBins,char bitRate,float maxDN,float tot)
{
  int i=0;
  int nDN=0;
  float *sampled=NULL;
  float resDN=0;

  sampled=falloc((uint64_t)nBins,"sampled wave",0);

  /*number of bins*/
  nDN=1;
  for(i=0;i<bitRate;i++)nDN*=2;

  resDN=maxDN/(float)nDN;
  for(i=0;i<nBins;i++)sampled[i]=floor(wave[i]/resDN)*resDN;

  return(sampled);
}/*digitiseWave*/


/*####################################################*/
/*scale noise to match Bryan's numbers*/

void scaleNoiseDN(float *noised,int nBins,float noiseSig,float trueSig,float offset)
{
  int i=0;
  float sigScale=0;

  sigScale=trueSig/noiseSig;

  for(i=0;i<nBins;i++)noised[i]=noised[i]*sigScale+offset;

  return;
}/*scaleNoiseDN*/


/*####################################################*/
/*noise from a Gaussian*/

float GaussNoise()
{
  float noise=0,max=0;
  float x1=0,x2=0,w=0;

  if(RAND_MAX>0)max=(float)RAND_MAX;
  else          max=-1.0*(float)RAND_MAX;

  /*Box approximation to Gaussian random number*/
  w=0.0;
  do{
    x1=2.0*(float)rand()/max-1.0;
    x2=2.0*(float)rand()/max-1.0;
    w=x1*x1+x2*x2;
  }while(w>=1.0);
  w=sqrt((-2.0*log(w))/w);

  noise=x1*w;

  return(noise);
}/*GaussNoise*/


/*####################################################*/
/*Remove noise on truth. Conservative*/

float *denoiseTruth(float *wave,int nBins)
{
  float *tempWave=NULL;
  void setDenoiseDefault(denPar *);
  denPar den;

  /*pick some denoising parameters, very conservative*/
  setDenoiseDefault(&den);
  den.varNoise=1;
  den.statsLen=15.0;
  den.noiseTrack=1;
  den.threshScale=4.0;

  /*denoise*/
  tempWave=processFloWave(wave,nBins,&den,1.0);

  return(tempWave);
}/*denoiseTruth*/


/*####################################################*/
/*modify the truth in terms of noise and pulse width*/

void modifyTruth(dataStruct *data,noisePar *gNoise)
{
  float *tempWave=NULL;
  float *denoiseTruth(float *,int);
  float sigDiff=0;
  char doNothing=0;

  /*remove noise*/
  tempWave=denoiseTruth(data->wave[data->useType],data->nBins);

  /*change pulse width*/
  if(gNoise->newPsig>0.0){
    if(gNoise->newPsig<data->pSigma){   /*reduce pulse width*/
      fprintf(stderr,"Can't deconvolve for new pulse length just yet. Old sigma %f new sigma %f\n",data->pSigma,gNoise->newPsig);
      exit(1);
    }else if(gNoise->newPsig>data->pSigma){  /*increase pulse width*/
      sigDiff=sqrt(gNoise->newPsig*gNoise->newPsig-data->pSigma*data->pSigma);
      TIDY(data->wave[data->useType]);
      data->wave[data->useType]=smooth(sigDiff,data->nBins,tempWave,data->res);
    }else{  /*do not change*/
      doNothing=1;
    }
  }else doNothing=1;


  if(doNothing==1){
    TIDY(data->wave[data->useType]);
    data->wave[data->useType]=tempWave;
    tempWave=NULL;
  }

  TIDY(tempWave);
  return;
}/*modifyTruth*/


/*####################################################*/
/*delete all signal beneath ground*/

void deleteGround(float *noised,float *wave,float *ground,int nBins,float minGap,float pSigma,float fSigma,float res,float trueCov,float rhoc,float rhog)
{
  int i=0;
  float maxGr=0;
  float tot=0,thresh=0;
  float rhoTot=0;
  float groundAmp=0,gRefl=0;
  float slope=0,tanSlope=0;
  float sigEff=0;


  if(minGap>0.0){  /*delete min detectable ground intensity*/
    tot=0.0;
    for(i=0;i<nBins;i++)tot+=wave[i];
    tot*=res;

    slope=2.0*M_PI/180.0;
    gRefl=minGap*rhog;
    tanSlope=sin(slope)/cos(slope);
    sigEff=sqrt(pSigma*pSigma+fSigma*fSigma*tanSlope*tanSlope);
    groundAmp=gRefl/(sigEff*sqrt(2.0*M_PI));

    if(trueCov<=0.0)rhoTot=rhog*minGap+rhoc*(1.0-minGap);
    else            rhoTot=rhog*(1.0-trueCov)+rhoc*trueCov;
    thresh=groundAmp*tot/rhoTot;

    for(i=0;i<nBins;i++){
      noised[i]=wave[i]-thresh;
      if(noised[i]<0.0)noised[i]=0.0;
    }
  }else{   /*delete based on ground intensity*/
    /*determine max ground intensity*/
    maxGr=-100.0;
    for(i=0;i<nBins;i++)if(ground[i]>maxGr)maxGr=ground[i];
    if(maxGr<0.0)maxGr=0.0;
    /*delete that from waveform*/
    for(i=0;i<nBins;i++){
      noised[i]=wave[i]-maxGr;
      if(noised[i]<0.0)noised[i]=0.0;
    }
  }

  return;
}/*deleteGround*/


/*####################################################*/
/*apply detector mean background drift*/

float *detectorDrift(float *wave,int nBins,float driftFact,float res)
{
  int i=0;
  float *tempWave=NULL;
  float cumul=0,tot=0;

  /*allocate*/
  tempWave=falloc((uint64_t)nBins,"drifted wave",0);

  /*do we need to apply drift?*/
  if(fabs(driftFact)>DRIFTTOL){
    /*total energy*/
    tot=0.0;
    for(i=0;i<nBins;i++)tot+=wave[i]*res;
    /*loop along wave*/
    cumul=0.0;
    for(i=0;i<nBins;i++){
      cumul+=wave[i]*res/tot;
      tempWave[i]=wave[i]-cumul*driftFact;
    }/*wave loop*/
  }else{  /*if not, copy waveform*/
    for(i=0;i<nBins;i++)tempWave[i]=wave[i];
  }/*apply or not switch*/

  return(tempWave);
}/*detectorDrift*/


/*the end*/
/*####################################################*/

