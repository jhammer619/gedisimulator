#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "libReadLVIS.h"


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
#define XRES 0.00005
#define YTOL 0.0000001
#define MINERR 0.0000001

/*element reflectance*/
float rhoG;
float rhoC;

/*####################################*/
/*control structure*/

typedef struct{
  /*input/output*/
  int nFiles;   /*number of waveforms*/
  char **inList;
  char listNamen[200];
  char outRoot[200];
  FILE *opooGauss;  /*Gaussian parameter output*/
  FILE *opooMet;    /*waveform metric output*/
  int maxGauss;     /*maximum number of Gaussians for output*/

  /*switches*/
  char ground;      /*read separateground wave or not*/
  char writeFit;    /*write fitted wave switch*/
  char useInt;      /*use discrete intensity instead of count*/
  char useFrac;     /*use fraction of hits per beam for weighting*/
  float rhRes;      /*rh resolution*/
  char bayesGround; /*Bayseian ground finding*/
  char noRHgauss;   /*do not do Gaussian fitting*/
  char renoiseWave; /*remove noise before adding*/
  char dontTrustGround; /*don't trust ground included with waveforms*/
  char readBinLVIS;  /*read binary LVIS rather than a list of ASCII files*/

  /*denoising parameters*/
  denPar *den;   /*for denoising*/
  denPar *gFit;  /*for Gaussian fitting*/

  /*noise parameters*/
  float meanN;
  float nSig;
  float bThresh;   /*bounds threshold*/
  float hNoise;    /*hard threshold noise as a fraction of integral*/
  char missGround; /*force to miss ground to get RH errors*/
  float minGap;    /*minimum detectavle gap fraction*/
  char linkNoise;  /*use link noise or not*/
  float linkM;     /*link margin*/
  float linkCov;   /*cover at which link margin is defined*/
  float linkSig;   /*link noise sigma*/
  float linkFsig;  /*footprint sigma used for link margin*/
  float linkPsig;  /*pulse sigma used for link margin*/
  float trueSig;   /*true noise sigma in DN*/
  float deSig;     /*detector sigma*/
  char bitRate;   /*digitiser bit rate*/
  float maxDN;    /*maximum DN we need to digitise*/
  float offset;   /*waveform DN offset*/

  /*pulse parameters*/
  float pSigma;    /*pulse length*/
  float fSigma;    /*footprint width*/
  float newPsig;   /*new pulse sigma*/

  /*LVIS data*/
  int verMaj;      /*major version*/
  int verMin;      /*minor version*/
  lvisStruct lvis;

  /*others*/
  float rhoRatio; /*ration of canopy to ground reflectance*/
  float res;      /*range resolution*/
  float gTol;     /*toleranve used to label ALS ground finding*/
  float zen;      /*zenith angle*/
  int nMessages;  /*number of progress messages*/
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
  float FHD;        /*foliage height diversity*/
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

  int nBgr;         /*number of ground estimates*/
  bGround *bGr;     /*Bayesian ground structure*/
  double bayGround; /*Bayesian ground elevation*/
}metStruct;


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
}dataStruct;


/*###########################################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct *data=NULL;
  dataStruct *readASCIIdata(char *,control *);
  dataStruct *readBinaryLVIS(char *,lvisStruct *,int,control *);
  metStruct *metric=NULL;
  void findMetrics(metStruct *,float *,int,float *,float *,int,double *,control *,dataStruct *);
  void tidySMoothPulse();
  void alignElevation(double,double,float *,int);
  void writeResults(dataStruct *,control *,metStruct *,int,float *,float *,char *);
  void addNoise(dataStruct *,control *);
  void determineTruth(dataStruct *,control *);
  void modifyTruth(dataStruct *,control *);
  float *processed=NULL,*denoised=NULL;;
  float setNoiseSigma(float,float,float,float);


  /*read command Line*/
  dimage=readCommands(argc,argv);

  /*set link noise if needed*/
  dimage->linkSig=setNoiseSigma(dimage->linkM,dimage->linkCov,dimage->pSigma,dimage->fSigma);

  /*allocate metric array*/
  if(!(metric=(metStruct *)calloc(1,sizeof(metStruct)))){
    fprintf(stderr,"error metric structure allocation.\n");
    exit(1);
  }


  /*loop over files*/
  for(i=0;i<dimage->nFiles;i++){
    if((i%dimage->nMessages)==0)fprintf(stdout,"Wave %d of %d\n",i+1,dimage->nFiles);

    /*read waveform*/
    if(dimage->readBinLVIS)data=readBinaryLVIS(dimage->inList[0],&dimage->lvis,i,dimage);
    else                   data=readASCIIdata(dimage->inList[i],dimage);

    /*is the data usable*/
    if(data->usable){
      /*adjust link noise if needed*/
      if(dimage->linkNoise&&((dimage->pSigma!=data->pSigma)||(dimage->fSigma!=data->fSigma))){
        dimage->linkSig=setNoiseSigma(dimage->linkM,dimage->linkCov,dimage->linkPsig,dimage->linkFsig);
        dimage->pSigma=data->pSigma;
        dimage->fSigma=data->fSigma;
      }

      /*denoise and change pulse if needed*/
      if(dimage->renoiseWave)modifyTruth(data,dimage);

      /*determine truths before noising*/
      determineTruth(data,dimage);

      /*add noise if needed*/
      addNoise(data,dimage);

      /*process waveform*/
      /*denoise*/
      denoised=processFloWave(data->noised,data->nBins,dimage->den,1.0);

      /*Gaussian fit*/
      if(dimage->noRHgauss==0)processed=processFloWave(denoised,data->nBins,dimage->gFit,1.0);

      /*shift Gaussian centres to align to absolute elevation*/
      alignElevation(data->z[0],data->z[data->nBins-1],dimage->gFit->gPar,dimage->gFit->nGauss);

      /*determine metrics*/
      findMetrics(metric,dimage->gFit->gPar,dimage->gFit->nGauss,denoised,data->noised,data->nBins,data->z,dimage,data);

      /*write results*/
      if(dimage->readBinLVIS)writeResults(data,dimage,metric,i,denoised,processed,dimage->inList[0]);
      else                   writeResults(data,dimage,metric,i,denoised,processed,dimage->inList[i]);
    }/*is the data usable*/


    /*tidy as we go along*/
    TIDY(processed);
    TIDY(denoised);
    if(data){
      TIDY(data->noised);
      TIDY(data->ground);
      TIDY(data->wave);
      TIDY(data->z);
      TIDY(data);
    }
    TIDY(dimage->gFit->gPar);
    TIDY(dimage->den->gPar);
    dimage->den->nGauss=0;
    dimage->gFit->nGauss=0;
    TIDY(metric->rhMax);
    TIDY(metric->rhInfl);
    TIDY(metric->rhReal);
    TIDY(metric->rh);
    TIDY(metric->bGr);
    //TIDY(metric->LmomGau);
    //TIDY(metric->LmomRea);
    //TIDY(metric->LmomInf);
    //TIDY(metric->LmomMax);
  }/*file loop*/

  /*TIDY LVIS data if it was read*/
  if(dimage->readBinLVIS)TIDY(dimage->lvis.data);

  fprintf(stdout,"Written to %s.gauss.txt\n",dimage->outRoot);
  fprintf(stdout,"Written to %s.metric.txt\n",dimage->outRoot);


  /*tidy up arrays*/
  tidySMoothPulse();
  TIDY(metric);
  if(dimage){
    if(dimage->readBinLVIS)TTIDY((void **)dimage->inList,1);
    else                   TTIDY((void **)dimage->inList,dimage->nFiles);
    dimage->inList=NULL;
    if(dimage->opooMet){
      fclose(dimage->opooMet);
      dimage->opooMet=NULL;
    }
    if(dimage->opooGauss){
      fclose(dimage->opooGauss);
      dimage->opooGauss=NULL;
    }
    if(dimage->den){
      TTIDY((void **)dimage->den->pulse,2);
      TIDY(dimage->den->matchPulse);
      TIDY(dimage->den->hardPulse);
      TIDY(dimage->den);
    }
    if(dimage->gFit){
      TTIDY((void **)dimage->gFit->pulse,2);
      TIDY(dimage->gFit);
    }
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################################*/
/*calculate sigma for link noise*/

float setNoiseSigma(float linkM,float cov,float pSigma,float fSigma)
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

  gRefl=(1.0-cov)*rhoG;

  tanSlope=sin(slope)/cos(slope);
  sigEff=sqrt(pSigma*pSigma+fSigma*fSigma*tanSlope*tanSlope);
  groundAmp=(gRefl/(gRefl+rhoC*cov))/(sigEff*sqrt(2.0*M_PI));  /*normalise by total waveform reflectance*/

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
  float nNsig=0,nSsig=0;
  float threshS=0,threshN=0;
  void gaussThresholds(float,float,float,float,float *,float *);

  char direction=0;

  /*initial guess*/
  sig=groundAmp/2.0;
  step=groundAmp/10.0;

  /*determine threshold in terms of standard deviations*/
  gaussThresholds(1.0,XRES,probNoise,probMiss,&nNsig,&nSsig);

  direction=0;
  minErr=0.00015;
  do{
    /*scale thresholds by sigma*/
    threshN=nNsig*sig;
    threshS=nSsig*sig;

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


/*####################################################*/
/*calculate Gaussian hresholds*/

void gaussThresholds(float sig,float res,float probNoise,float probMiss,float *threshN,float *threshS)
{
  float x=0,y=0;
  float cumul=0;
  char foundS=0,foundN=0;
  float probN=0,probS=0;

  probN=1.0-probNoise/(30.0/0.15);
  probS=1.0-probMiss;

  /*determine start*/
  x=0.0;
  do{
    y=(float)gaussian((double)x,(double)sig,0.0);
    x-=res;
  }while(y>=YTOL);

  do{
    y=(float)gaussian((double)x,(double)sig,0.0);
    cumul+=y*res;

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

  return;
}/*gaussThresholds*/


/*####################################################*/
/*determine true variables*/

void determineTruth(dataStruct *data,control *dimage)
{
  int i=0;
  float totE=0,cumul=0;
  float totC=0,totG=0;
  float groundOverlap(float *,float *,int);
  float groundMinAmp(float *,float *,int);
  float groundInflection(float *,float *,int);

  /*determine ground*/
  if(!data->demGround){  /*unless it's already been calculalated from the DEM*/
    if(dimage->ground){
      totE=0.0;
      data->gElev=0.0;
      for(i=0;i<data->nBins;i++){
        totE+=data->ground[i];
        data->gElev+=(double)data->ground[i]*data->z[i];
      }
      if(totE>0.0)data->gElev/=(double)totE;
      else        data->gElev=-1000000.0;

      /*standard deviation as a measure of slope*/
      data->gStdev=0.0;
      for(i=0;i<data->nBins;i++){
        data->gStdev+=(float)((data->z[i]-data->gElev)*(data->z[i]-data->gElev))*data->ground[i];
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
  for(i=0;i<data->nBins;i++)totE+=data->wave[i];
  cumul=0.0;
  for(i=0;i<data->nBins;i++){
    cumul+=data->wave[i]/totE;
    if(cumul>=dimage->bThresh){
      data->tElev=data->z[i];
      break;
    }
  }/*top finding*/

  /*canopy cover*/
  if(dimage->ground){
    totG=totC=0.0;
    for(i=0;i<data->nBins;i++){
     totG+=data->ground[i];
     totC+=data->wave[i]-data->ground[i];
    }
    if((totG+totC)>0.0)data->cov=totC/(totC+totG*dimage->rhoRatio);
    else               data->cov=-1.0;
  }else{
    data->cov=-1.0;
  }

  /*understorey metrics*/
  if(dimage->ground){
    data->gLap=groundOverlap(data->wave,data->ground,data->nBins);
    data->gMinimum=groundMinAmp(data->wave,data->ground,data->nBins);
    data->gInfl=groundInflection(data->wave,data->ground,data->nBins);
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
/*add noise to waveform*/

void addNoise(dataStruct *data,control *dimage)
{
  int i=0;
  float noise=0;
  float tot=0.0,thresh=0;
  float *tempNoise=NULL;
  float GaussNoise();
  float *smooNoise=NULL;
  float *digitiseWave(float *,int,char,float,float);
  float reflScale=0;
  void deleteGround(float *,float *,float *,int,float,float,float,float,float);
  void scaleNoiseDN(float *,int,float,float,float);

  /*allocate*/
  data->noised=falloc(data->nBins,"noised wave",0);

  if(dimage->missGround){        /*Delete all signal beneath ground peak*/
    if((dimage->ground==0)&&(dimage->minGap==0.0)){
      fprintf(stderr,"Cannot delete ground when we do not know it\n");
      exit(1);
    }
    deleteGround(data->noised,data->wave,data->ground,data->nBins,dimage->minGap,dimage->pSigma,dimage->fSigma,dimage->res,data->cov);
  }else if(dimage->linkNoise){   /*link margin based noise*/
    /*Gaussian noise*/
    tempNoise=falloc(data->nBins,"temp noised",0);
    tot=0.0;
    for(i=0;i<data->nBins;i++)tot+=data->wave[i]*dimage->res;  
    reflScale=(data->cov*rhoC+(1.0-data->cov)*rhoG)*tot/(dimage->linkCov*rhoC+(1.0-dimage->linkCov)*rhoG);
    for(i=0;i<data->nBins;i++)tempNoise[i]=dimage->linkSig*GaussNoise()*reflScale;
    /*smooth noise by detector response*/
    smooNoise=smooth(dimage->deSig,data->nBins,tempNoise,dimage->res);
    for(i=0;i<data->nBins;i++)tempNoise[i]=data->wave[i]+smooNoise[i];
    TIDY(smooNoise);
    /*scale to match sigma*/
    scaleNoiseDN(tempNoise,data->nBins,dimage->linkSig*reflScale,dimage->trueSig,dimage->offset);
    /*digitise*/
    TIDY(data->noised);
    data->noised=digitiseWave(tempNoise,data->nBins,dimage->bitRate,dimage->maxDN,tot);
    TIDY(tempNoise);
  }else if((dimage->nSig>0.0)||(dimage->meanN>0.0)){   /*mean and stdev based noise*/
    for(i=0;i<data->nBins;i++){
      noise=dimage->nSig*GaussNoise();
      if((float)rand()/(float)RAND_MAX<0.5)noise*=-1.0; /*to allow negative numbers*/
      data->noised[i]=data->wave[i]+dimage->meanN+noise;
    }/*bin loop*/
  }else if(dimage->hNoise>0.0){  /*hard threshold noise*/
    tot=0.0;
    for(i=0;i<data->nBins;i++)tot+=data->wave[i];
    thresh=dimage->hNoise*tot;
    for(i=0;i<data->nBins;i++){
      data->noised[i]=data->wave[i]-thresh;
      if(data->noised[i]<0.0)data->noised[i]=0.0;
    }
  }else{  /*no noise*/
    for(i=0;i<data->nBins;i++)data->noised[i]=data->wave[i];
  }

  return;
}/*addNoise*/


/*####################################################*/
/*modify the truth in terms of noise and pulse width*/

void modifyTruth(dataStruct *data,control *dimage)
{
  float *tempWave=NULL;
  float *denoiseTruth(float *,int,control *);
  float sigDiff=0;
  char doNothing=0;

  /*remove noise*/
  tempWave=denoiseTruth(data->wave,data->nBins,dimage);

  /*change pulse width*/
  if(dimage->newPsig>0.0){
    if(dimage->newPsig<data->pSigma){   /*reduce pulse width*/
      fprintf(stderr,"Can't deconvolve for new pulse length just yet\n");
      exit(1);
    }else if(dimage->newPsig>data->pSigma){  /*increase pulse width*/
      sigDiff=sqrt(dimage->newPsig*dimage->newPsig-data->pSigma*data->pSigma);
      TIDY(data->wave);
      data->wave=smooth(sigDiff,data->nBins,tempWave,data->res);
    }else{  /*do not change*/
      doNothing=1;
    }
  }else doNothing=1;


  if(doNothing==1){
    TIDY(data->wave);
    data->wave=tempWave;
    tempWave=NULL;
  }

  TIDY(tempWave);
  return;
}/*modifyTruth*/


/*####################################################*/
/*Remove noise on truth. Conservative*/

float *denoiseTruth(float *wave,int nBins,control *dimage)
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
/*delete all signal beneath ground*/

void deleteGround(float *noised,float *wave,float *ground,int nBins,float minGap,float pSigma,float fSigma,float res,float trueCov)
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
    gRefl=minGap*rhoG;
    tanSlope=sin(slope)/cos(slope);
    sigEff=sqrt(pSigma*pSigma+fSigma*fSigma*tanSlope*tanSlope);
    groundAmp=gRefl/(sigEff*sqrt(2.0*M_PI));

    if(trueCov<=0.0)rhoTot=rhoG*minGap+rhoC*(1.0-minGap);
    else            rhoTot=rhoG*(1.0-trueCov)+rhoC*trueCov;
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
/*digitse*/

float *digitiseWave(float *wave,int nBins,char bitRate,float maxDN,float tot)
{
  int i=0;
  int nDN=0;
  float *sampled=NULL;
  float resDN=0;

  sampled=falloc(nBins,"sampled wave",0);

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
/*write results*/

void writeResults(dataStruct *data,control *dimage,metStruct *metric,int numb,float *denoised,float *processed,char *inNamen)
{
  int i=0,j=0;
  char waveNamen[200];
  char namen[200];
  float gauss(float,float,float);
  FILE *opoo=NULL;


  /*open file if needed*/
  if(dimage->opooGauss==NULL){
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
    //for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet,", %d LmomGauss%d",14+4*metric->nRH+1+dimage->bayesGround+21+i,i+1);
    //for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet,", %d LmomInfl%d",14+4*metric->nRH+1+dimage->bayesGround+21+metric->nLm+i,i+1);
    //for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet,", %d LmomMax%d",14+4*metric->nRH+1+dimage->bayesGround+21+2*metric->nLm+i,i+1);
    //for(i=0;i<metric->nLm;i++)fprintf(dimage->opooMet,", %d LmomReal%d",14+4*metric->nRH+1+dimage->bayesGround+21+3*metric->nLm+i,i+1);
    fprintf(dimage->opooMet,",\n");
  }

  if(dimage->gFit->nGauss>dimage->maxGauss)fprintf(stderr,"More Gaussians than header entries %d\n",dimage->gFit->nGauss);

  /*fitted Gaussians*/
  fprintf(dimage->opooGauss,"%d %d",numb,dimage->gFit->nGauss);
  for(i=0;i<dimage->gFit->nGauss;i++){
    if((dimage->gFit->gPar[3*i]>=0.0)&&(dimage->gFit->gPar[3*i+1]>=0.0)&&(dimage->gFit->gPar[3*i+2]>=0.0)){
      fprintf(dimage->opooGauss," %f %f %f",dimage->gFit->gPar[3*i],dimage->gFit->gPar[3*i+1],dimage->gFit->gPar[3*i+2]);
    }else fprintf(dimage->opooGauss," ? ? ?");
  }
  for(i=dimage->gFit->nGauss;i<dimage->maxGauss;i++)fprintf(dimage->opooGauss," ? ? ?");
  fprintf(dimage->opooGauss," %s\n",inNamen);


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
  if(dimage->linkNoise)fprintf(dimage->opooMet," %f %f",dimage->linkM,dimage->linkCov);
  else                 fprintf(dimage->opooMet," ? ?");
  fprintf(dimage->opooMet," %.2f %.2f",data->lon,data->lat);
  fprintf(dimage->opooMet," %f %f %f",data->gLap,data->gMinimum,data->gInfl); 
  fprintf(dimage->opooMet," %f %f",metric->totE,metric->blairSense);
  fprintf(dimage->opooMet," %f %f",data->pointDense,data->beamDense);
  fprintf(dimage->opooMet," %f %f",data->zen,metric->FHD);
  fprintf(dimage->opooMet," %f %f",metric->niM2,metric->niM21);
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
    for(i=0;i<dimage->gFit->nGauss;i++)fprintf(opoo,", %d Gauss %d",i+8,i+1);
    fprintf(opoo,"\n");
    fprintf(opoo,"# fSigma %f pSigma %f res %f\n",data->fSigma,data->pSigma,dimage->res);
    fprintf(opoo,"# coord %.2f %.2f\n",data->lon,data->lat);
    fprintf(opoo,"# cover %f rhoG %f rhoC %f\n",data->cov,rhoG,rhoC);
    fprintf(opoo,"# ground %.2f slope %f\n",data->gElev,data->slope);
    for(i=0;i<data->nBins;i++){
      fprintf(opoo,"%f %f %f %f %f",data->z[i],data->noised[i],denoised[i],processed[i],data->wave[i]);
      if(dimage->ground)fprintf(opoo," %f %f",data->ground[i],data->wave[i]-data->ground[i]);
      else              fprintf(opoo," 0 0");
      for(j=0;j<dimage->gFit->nGauss;j++)fprintf(opoo," %f",dimage->gFit->gPar[j*3+1]*gauss((float)data->z[i],dimage->gFit->gPar[3*j+2],dimage->gFit->gPar[3*j]));
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

void findMetrics(metStruct *metric,float *gPar,int nGauss,float *processed,float *wave,int nBins,double *z,control *dimage,dataStruct *data)
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
  float *smoothed=NULL;
  float halfCover(float *,double *,int,double,float);
  float findBlairSense(metStruct *,dataStruct *,control *);
  void findSignalBounds(float *,double *,int,double *,double *,control *);
  void findWaveExtents(float *,double *,int,double,double,float *,float *);
  void setDenoiseDefault(denPar *);
  denPar den;


  /*total energy. The sqrt(2.pi) can be normalised out*/
  tot=0.0;
  energy=falloc(nGauss,"energy",0);
  mu=falloc(nGauss,"mu",0);
  A=falloc(nGauss,"A",0);
  sig=falloc(nGauss,"sigma",0);
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
  for(i=0;i<nBins;i++)metric->totE+=processed[i]*dimage->res;

  /*BLair sensitivity*/
  metric->blairSense=findBlairSense(metric,data,dimage);

  /*smooth waveform for fonding ground by max and inflection*/
  setDenoiseDefault(&den);
  den.varNoise=0;
  den.meanN=0;
  den.thresh=0;
  den.sWidth=0.76*3.0/4.0;  /*according to Bryan Blair*/
  den.noiseTrack=0;
  smoothed=processFloWave(processed,nBins,&den,1.0);

  /*ground by Gaussian fit*/
  if(dimage->noRHgauss==0)metric->gHeight=gaussianGround(energy,mu,&gInd,nGauss,tot);
  else                    metric->gHeight=-1.0;

  /*canopy cover*/
  metric->cov=gaussCover(processed,nBins,mu,energy,nGauss,gInd);

  /*ground by maximum*/
  metric->maxGround=maxGround(smoothed,z,nBins);

  /*ground by inflection*/
  metric->inflGround=inflGround(smoothed,z,nBins);

  /*rh metrics with Gaussian ground*/
  if(dimage->noRHgauss==0)metric->rh=findRH(processed,z,nBins,metric->gHeight,dimage->rhRes,&metric->nRH);
  else                    metric->rh=blankRH(dimage->rhRes,&metric->nRH);

  /*rh metrics with maximum ground*/
  metric->rhMax=findRH(processed,z,nBins,metric->maxGround,dimage->rhRes,&metric->nRH);

  /*rh metrics with inflection ground*/
  metric->rhInfl=findRH(processed,z,nBins,metric->inflGround,dimage->rhRes,&metric->nRH);

  /*rh metrics with real ground, if we have the ground*/
  if(dimage->ground)metric->rhReal=findRH(data->wave,z,nBins,data->gElev,dimage->rhRes,&metric->nRH);
  else{
    metric->rhReal=falloc(metric->nRH,"rhReal",0);
    for(i=0;i<metric->nRH;i++)metric->rhReal[i]=-1.0;
  }

  /*foliage height diversity*/
  metric->FHD=foliageHeightDiversity(processed,nBins);

  /*signal start and end*/
  findSignalBounds(processed,z,nBins,&metric->tElev,&metric->bElev,dimage);

  /*Lefsky's leading and trailing edge extents*/
  findWaveExtents(processed,z,nBins,metric->tElev,metric->bElev,&metric->leExt,&metric->teExt);

  /*L moments*/
  /*metric->nLm=4;
  metric->LmomGau=waveLmoments(metric->rh,metric->nRH,dimage->rhRes,metric->nLm);
  metric->LmomRea=waveLmoments(metric->rhReal,metric->nRH,dimage->rhRes,metric->nLm);
  metric->LmomInf=waveLmoments(metric->rhInfl,metric->nRH,dimage->rhRes,metric->nLm);
  metric->LmomMax=waveLmoments(metric->rhMax,metric->nRH,dimage->rhRes,metric->nLm);*/

  /*cover estimates using Bryan's half feature method*/
  metric->covHalfG=halfCover(processed,z,nBins,metric->gHeight,dimage->rhoRatio);
  metric->covHalfM=halfCover(processed,z,nBins,metric->maxGround,dimage->rhoRatio);
  metric->covHalfI=halfCover(processed,z,nBins,metric->inflGround,dimage->rhoRatio);

  /*Wenge Ni's metrics*/
  metric->niM2=niMetric(processed,z,nBins,dimage->res,metric->gHeight,2.0);
  metric->niM21=niMetric(processed,z,nBins,dimage->res,metric->gHeight,2.1);

  /*bayesian ground finding*/
  if(dimage->bayesGround){
    metric->bayGround=bayesGround(wave,nBins,dimage,metric,z,data);
    metric->covHalfB=halfCover(processed,z,nBins,metric->bayGround,dimage->rhoRatio);
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
/*determine Blair sensitivity metric*/

float findBlairSense(metStruct *metric,dataStruct *data,control *dimage)
{
  float gAmp=0;
  float sigEff=0,gArea=0;
  float slope=0,tanSlope=0;
  float blairSense=0;
  float meanN=0,stdev=0;
  float notNeeded=0;
  float nNsig=0,nSsig=0;
  float probNoise=0,probMiss=0;
  void meanNoiseStats(float *,uint32_t,float *,float *,float *,float,float,int);
  void gaussThresholds(float,float,float,float,float *,float *);

  /*determine noise stats for sensitivity metric*/
  meanNoiseStats(data->noised,(uint32_t)data->nBins,&meanN,&stdev,&notNeeded,-1.0,1.0,(int)(dimage->den->statsLen/dimage->res));
  stdev-=meanN;

  if(stdev>0.0){
    probNoise=0.05;
    probMiss=0.1;
    gaussThresholds(1.0,XRES,probNoise,probMiss,&nNsig,&nSsig);

    slope=2.0*M_PI/180.0;
    tanSlope=sin(slope)/cos(slope);
    gAmp=(nNsig+nSsig)*stdev;
    sigEff=sqrt(dimage->linkPsig*dimage->linkPsig+dimage->linkFsig*dimage->linkFsig*tanSlope*tanSlope);
    gArea=gAmp*sigEff*sqrt(2.0*M_PI)/metric->totE;

    if(gArea>0.0)blairSense=1.0-gArea;
    else         blairSense=0.0;
  }else blairSense=1.0;

  return(blairSense);
}/*findBlairSense*/


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
    den[i].gWidth=dimage->den->sWidth;
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
    gaussProps(den[i].gPar,den[i].nGauss,dimage->fSigma,sqrt(data->pSigma*data->pSigma+den[i].sWidth*den[i].sWidth),&(metric->bGr[i].gHeight),&(metric->bGr[i].slope),&(metric->bGr[i].cov),&height);

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
/*read LVIS binary data*/

dataStruct *readBinaryLVIS(char *namen,lvisStruct *lvis,int numb,control *dimage)
{
  int i=0;
  dataStruct *data=NULL;
  float pulseLenFromTX(unsigned char *,int);

  /*do we need to read all the data*/
  if(lvis->data==NULL){
    lvis->verMaj=dimage->verMaj;
    lvis->verMin=dimage->verMin;
    /*check version number*/
    if((lvis->verMaj!=1)||(lvis->verMin!=3)){
      fprintf(stderr,"Currently only works for version 1.3, not %d.%d\n",lvis->verMaj,lvis->verMin);
      exit(1);
    }
    /*read data*/
    lvis->nBins=432;
    lvis->data=readLVISdata(namen,&lvis->nWaves);
    dimage->nFiles=lvis->nWaves;
  }


  /*allocate space*/
  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  data->useID=0;
  data->nBins=lvis->nBins;
  data->wave=falloc(data->nBins,"waveform",0);
  data->z=dalloc(data->nBins,"z",0);
  data->ground=NULL;
  data->useID=0;
  data->demGround=0;
  data->pSigma=-1.0;    /*nonesense pulse length*/
  data->fSigma=-1.0;    /*nonesense footprint width*/
  data->usable=1;

  /*copy data to structure*/
  data->zen=lvis->data[numb].zen;
  data->res=dimage->den->res=dimage->gFit->res=(lvis->data[numb].z0-lvis->data[numb].z431)/(float)lvis->nBins;
  data->totE=0.0;
  for(i=0;i<lvis->nBins;i++){
    data->wave[i]=(float)lvis->data[numb].rxwave[i];
    data->z[i]=(double)(lvis->data[numb].z0-(float)i*data->res);
    data->totE+=data->wave[i];
  }
  if(dimage->den->res<TOL)data->usable=0;
  if(data->totE<=0.0)data->usable=0;
  if(dimage->ground==0){   /*set to blank*/
    data->cov=-1.0;
    data->gLap=-1.0;
    data->gMinimum=-1.0;
    data->gInfl=-1.0;
  }
  data->lon=(lvis->data[numb].lon0+lvis->data[numb].lon431)/2.0;
  data->lat=(lvis->data[numb].lat0+lvis->data[numb].lat431)/2.0;

  /*analyse pulse*/
  data->pSigma=pulseLenFromTX(lvis->data[numb].rxwave,80);

  if(data->fSigma<0.0)data->fSigma=dimage->fSigma;
  if(data->pSigma<0.0)data->pSigma=dimage->pSigma;


  /*set up number of messages*/
  if(lvis->nWaves>dimage->nMessages)dimage->nMessages=(int)(lvis->nWaves/dimage->nMessages);
  else                              dimage->nMessages=1;

  return(data);
}/*readBinaryLVIS*/


/*####################################################*/
/*determine pulse width from TXwave*/

float pulseLenFromTX(unsigned char *rxwave,int nBins)
{
  int i=0;
  float pSigma=0;
  float *temp=NULL;
  float *denoised=NULL;
  float tot=0,CofG=0;
  denPar den;
  void setDenoiseDefault(denPar *);

  /*copy to float*/
  temp=falloc(nBins,"temp pulse",0);
  for(i=0;i<nBins;i++)temp[i]=(float)rxwave[i];

  /*denoise*/
  setDenoiseDefault(&den);
  den.varNoise=1;
  den.threshScale=5.0;
  den.noiseTrack=1;
  den.minWidth=5;
  den.statsLen=3.0;
  den.res=0.15;
  denoised=processFloWave(temp,nBins,&den,1.0);

  /*CofG*/
  CofG=tot=0.0;
  for(i=0;i<nBins;i++){
    CofG+=(float)i*den.res*denoised[i];
    tot+=denoised[i];
  }
  if(tot>0.0){
    CofG/=tot;

    pSigma=0.0;
    for(i=0;i<nBins;i++)pSigma+=((float)i*den.res*denoised[i]-CofG)*((float)i*den.res*denoised[i]-CofG);
    pSigma=sqrt(pSigma/tot);

  }else pSigma=-1.0;

  TIDY(denoised);
  TIDY(temp);
  return(pSigma);
}/*pulseLenFromTX*/


/*####################################################*/
/*read ASCII data*/

dataStruct *readASCIIdata(char *namen,control *dimage)
{
  int i=0;
  dataStruct *data=NULL;
  char line[400],temp1[100],temp2[100];
  char temp3[100],temp4[100],temp5[100];
  char temp6[100],temp7[100],temp8[100];
  char temp9[100],temp10[100];
  FILE *ipoo=NULL;

  /*open input*/
  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  data->pSigma=-1.0;    /*nonesense pulse length*/
  data->fSigma=-1.0;    /*nonesense footprint width*/
  data->usable=1;
  data->zen=0.0;        /*straight down*/

  /*count number of wavebins*/
  data->nBins=0;
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))data->nBins++;

  data->wave=falloc(data->nBins,"waveform",0);
  data->z=dalloc(data->nBins,"z",0);
  if(dimage->ground)data->ground=falloc(data->nBins,"ground",0);
  data->useID=0;
  data->demGround=0;

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  i=0;
  while(fgets(line,400,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(dimage->useInt){  /*read intensity*/
        if(dimage->ground==0){   /*don't read ground*/
          if(sscanf(line,"%s %s",temp1,temp2)==2){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp2);
            i++;
          }
        }else{                   /*read ground*/
          if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp2);
            data->ground[i]=atof(temp4);
            i++;
          }
        }/*ground switch*/
      }else if(dimage->useFrac){ /*read fraction*/
        if(dimage->ground==0){   /*don't read ground*/
         if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp3);
            i++;
          }
        }else{                   /*read ground*/
          if(sscanf(line,"%s %s %s %s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10)==10){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp8);
            data->ground[i]=atof(temp10);
            i++;
          }
        }/*ground switch*/
      }else{ /*read count*/
        if(dimage->ground==0){   /*don't read ground*/
          if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp3);
            i++;
          }
        }else{                   /*read ground*/
          if(sscanf(line,"%s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7)==7){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp5);
            data->ground[i]=atof(temp7);
            i++;
          }
        }/*ground switch*/
      }/*intensity or count switch*/
    }else{  /*read the header*/
      if(sscanf(line,"%s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7)==7){
        if(!strncasecmp(temp2,"fSigma",6)){
          data->fSigma=atof(temp3);
          data->pSigma=atof(temp5);
        }
      }
      if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
        if(!strncasecmp(temp2,"waveID",6)){
          data->useID=1;
          strcpy(&(data->waveID[0]),temp3);
        }else if(!strncasecmp(temp2,"meanScanAng",11)){
          data->zen=atof(temp3);
        }
      }
      if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
        if(!strncasecmp(temp2,"coord",5)){
          data->lon=atof(temp3);
          data->lat=atof(temp4);
        }else if(!strncasecmp(temp2,"ground",6)){
          if(!dimage->dontTrustGround){
            data->gElev=atof(temp3);
            data->slope=atof(temp4);
            data->demGround=1;
          }else data->demGround=0;
        }
      }
      if(sscanf(line,"%s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6)==6){
        if(!strncasecmp(temp2,"density",7)){
          data->pointDense=atof(temp4);
          data->beamDense=atof(temp6);
        }else if(!strncasecmp(temp2,"lvis",4)){
          data->res=atof(temp4);
          data->zen=atof(temp6);
        }
      }
    }
  }/*line loop*/

  if(data->res<=0.0)data->res=dimage->res;

  /*add up energy*/
  data->totE=0.0;
  for(i=0;i<data->nBins;i++)data->totE+=data->wave[i];


  dimage->res=dimage->den->res=dimage->gFit->res=fabs(data->z[1]-data->z[0]);
  if(dimage->den->res<TOL)data->usable=0;
  if(data->totE<=0.0)data->usable=0;
  if(dimage->ground==0){   /*set to blank*/
    data->cov=-1.0;
    data->gLap=-1.0;
    data->gMinimum=-1.0;
    data->gInfl=-1.0;
  }
  if(data->fSigma<0.0)data->fSigma=dimage->fSigma;
  if(data->pSigma<0.0)data->pSigma=dimage->pSigma;

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  /*set up number of messages*/
  if(dimage->nFiles>dimage->nMessages)dimage->nMessages=(int)(dimage->nFiles/dimage->nMessages);
  else                                dimage->nMessages=1;

  return(data);
}/*readASCIIdata*/


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
  if(!(dimage->den=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  if(!(dimage->gFit=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }


  /*defaults*/
  /*input/output*/
  dimage->nFiles=1;
  dimage->inList=chChalloc(dimage->nFiles,"inList",0);
  dimage->inList[0]=challoc(200,"inList",0);
  strcpy(&(dimage->inList[0][0]),"/Users/stevenhancock/data/gedi/analysis/side_lobe/laselva/contLaselva.0.12.wave");
  strcpy(dimage->outRoot,"teastMetric");
  dimage->maxGauss=20;
  dimage->opooMet=NULL;
  dimage->opooGauss=NULL;
  /*scan settings*/
  dimage->pSigma=0.764331; /*pulse length*/
  dimage->fSigma=5.5;      /*footprint width*/
  dimage->res=0.15;

  /*switches*/
  dimage->writeFit=0;
  dimage->ground=0;
  dimage->useInt=0;
  dimage->useFrac=0;
  dimage->rhRes=10.0;
  dimage->bayesGround=0;
  dimage->missGround=0;
  dimage->linkNoise=0;
  dimage->linkPsig=0.764331; /*pulse length*/
  dimage->linkFsig=5.5;      /*footprint width*/
  dimage->trueSig=5.0;
  dimage->deSig=0.0; //0.1; //4.0*0.15/2.355;
  dimage->bitRate=12;
  dimage->maxDN=4096.0; //1.0/(dimage->pSigma*sqrt(2.0*M_PI));
  dimage->minGap=0.0;
  dimage->noRHgauss=0;    /*do find RH metrics by Gaussian fitting*/
  dimage->renoiseWave=0;  /*do not denoise "truth"*/
  dimage->newPsig=-1.0;   /*leave blank*/
  dimage->dontTrustGround=0;  /*do trust ground in waveforms, if there*/
  dimage->readBinLVIS=0;          /*read ASCII rather than binary LVIS*/

  /*set default denoising parameters*/
  setDenoiseDefault(dimage->den);
  dimage->den->meanN=0.0;  /*we haven't added noise yet*/
  dimage->den->thresh=0.00000001;  /*tiny number as no noise yet*/
  dimage->den->noiseTrack=0;
  dimage->den->minWidth=0;
  dimage->den->varNoise=0;
  dimage->den->threshScale=1.5;
  dimage->den->fitGauss=0;
  dimage->den->psWidth=0.0;

  /*set default Gaussian fitting  parameters*/
  setDenoiseDefault(dimage->gFit);
  dimage->gFit->meanN=0.0;    /*no denoising here*/
  dimage->gFit->thresh=0.000000005;    /*no denoising here*/
  dimage->gFit->noiseTrack=0;    /*no denoising here*/
  dimage->gFit->minWidth=0;    /*no denoising here*/
  dimage->gFit->varNoise=0;    /*no denoising here*/
  dimage->gFit->gWidth=1.2;
  dimage->gFit->sWidth=0.0;
  dimage->gFit->fitGauss=1;
  dimage->gFit->minGsig=0.9;
  /*noise parameters*/
  dimage->meanN=0.0;
  dimage->nSig=0.0;
  dimage->bThresh=0.001;
  dimage->hNoise=0.0;
  dimage->offset=94.0;
  /*LVIS data*/
  dimage->verMaj=1;  /*major version*/
  dimage->verMin=3;  /*minor version*/
  dimage->lvis.data=NULL;
  /*others*/
  dimage->rhoRatio=0.57/0.4;
  rhoG=0.4;
  rhoC=0.57;
  dimage->gTol=0.0;
  dimage->nMessages=200;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->nFiles=1;
        dimage->inList=chChalloc(dimage->nFiles,"input name list",0);
        dimage->inList[0]=challoc(strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=readInList(&dimage->nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(dimage->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-writeFit",9)){
        dimage->writeFit=1;
      }else if(!strncasecmp(argv[i],"-dcBias",7)){
        checkArguments(1,i,argc,"-dcBias");
        dimage->meanN=dimage->offset=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-hNoise",7)){
        checkArguments(1,i,argc,"-hNoise");
        dimage->hNoise=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-nSig",5)){
        checkArguments(1,i,argc,"-nSig");
        dimage->nSig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-seed",5)){
        checkArguments(1,i,argc,"-seed");
        srand(atoi(argv[++i]));
      }else if(!strncasecmp(argv[i],"-meanN",5)){
        checkArguments(1,i,argc,"-meanN");
        dimage->den->meanN=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-thresh",6)){
        checkArguments(1,i,argc,"-thresh");
        dimage->den->thresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-sWidth",7)){
        checkArguments(1,i,argc,"-sWidth");
        dimage->den->sWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-psWidth",8)){
        checkArguments(1,i,argc,"-psWidth");
        dimage->den->psWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gWidth",7)){
        checkArguments(1,i,argc,"-gWidth");
        dimage->gFit->gWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minGsig",8)){
        checkArguments(1,i,argc,"-minGsig");
        dimage->gFit->minGsig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minWidth",9)){
        checkArguments(1,i,argc,"-minWidth");
        dimage->den->minWidth=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-ground",7)){
        dimage->ground=1;
      }else if(!strncasecmp(argv[i],"-varNoise",9)){
        dimage->den->varNoise=1;
      }else if(!strncasecmp(argv[i],"-medNoise",9)){
        dimage->den->medStats=1;
      }else if(!strncasecmp(argv[i],"-statsLen",9)){
        checkArguments(1,i,argc,"-statsLen");
        dimage->den->statsLen=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-varScale",9)){
        checkArguments(1,i,argc,"-varScale");
        dimage->den->threshScale=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noiseTrack",11)){
        dimage->den->noiseTrack=1;
      }else if(!strncasecmp(argv[i],"-pFile",6)){
        checkArguments(1,i,argc,"-pFile");
        strcpy(dimage->den->pNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-gold",5)){
        dimage->den->deconMeth=0;
      }else if(!strncasecmp(argv[i],"-deconTol",9)){
        checkArguments(1,i,argc,"-deconTol");
        dimage->den->deChang=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-preMatchF",10)){
        dimage->den->preMatchF=1;
      }else if(!strncasecmp(argv[i],"-postMatchF",11)){
        dimage->den->posMatchF=1;
      }else if(!strncasecmp(argv[i],"-useInt",7)){
        dimage->useInt=1;
      }else if(!strncasecmp(argv[i],"-useFrac",8)){
        dimage->useFrac=1;
      }else if(!strncasecmp(argv[i],"-rhRes",6)){
        checkArguments(1,i,argc,"-rhRes");
        dimage->rhRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-linkNoise",10)){
        checkArguments(2,i,argc,"-linkNoise");
        dimage->linkNoise=1;
        dimage->linkM=atof(argv[++i]);
        dimage->linkCov=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-trueSig",8)){
        checkArguments(1,i,argc,"-trueSig");
        dimage->trueSig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-missGround",11)){
        dimage->missGround=1;
      }else if(!strncasecmp(argv[i],"-minGap",7)){
        checkArguments(1,i,argc,"-minGap");
        dimage->missGround=1;
        dimage->minGap=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-bayesGround",13)){
        dimage->bayesGround=1;
      }else if(!strncasecmp(argv[i],"-rhoG",5)){
        checkArguments(1,i,argc,"-rhoG");
        rhoG=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-rhoC",5)){
        checkArguments(1,i,argc,"-rhoC");
        rhoC=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-offset",7)){
        checkArguments(1,i,argc,"-offset");
        dimage->offset=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxDN",6)){
        checkArguments(1,i,argc,"-maxDN");
        dimage->maxDN=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-bitRate",8)){
        checkArguments(1,i,argc,"-bitRate");
        dimage->bitRate=(char)atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gTol",5)){
        checkArguments(1,i,argc,"-gTol");
        dimage->gTol=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noRHgauss",10)){
        dimage->noRHgauss=1;
      }else if(!strncasecmp(argv[i],"-renoise",8)){
       dimage->renoiseWave=1;
      }else if(!strncasecmp(argv[i],"-newPsig",8)){
        checkArguments(1,i,argc,"-newPsig");
        dimage->newPsig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-oldPsig",8)){
        checkArguments(1,i,argc,"-oldPsig");
        dimage->pSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-dontTrustGround",16)){
        dimage->dontTrustGround=1;
      }else if(!strncasecmp(argv[i],"-readBinLVIS",12)){
        dimage->readBinLVIS=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to calculate GEDI waveform metrics\n#####\n\n-input name;     waveform  input filename\n-outRoot name;   output filename root\n-inList list;    input file list for multiple files\n-writeFit;       write fitted waveform\n-ground;         read true ground from file\n-useInt;         use discrete intensity instead of count\n-useFrac;        use fractional hits rather than counts\n-readBinLVIS;    input is an LVIS binary file\n-rhRes r;        percentage energy resolution of RH metrics\n-bayesGround;    use Bayseian ground finding\n-gTol tol;       ALS ground tolerance. Used to calculate slope.\n-noRHgauss; do not fit Gaussians\n-dontTrustGround;     don't trust ground in waveforms, if included\n\nAdding noise:\n-dcBias n;       mean noise level\n-nSig sig;       noise sigma\n-seed n;         random number seed\n-hNoise n;       hard threshold noise as a fraction of integral\n-linkNoise linkM cov;     apply Gaussian noise based on link margin at a cover\n-trueSig sig;    true sigma of background noise\n-renoise;        remove noise feom truth\n-newPsig sig;    new value for pulse width\n-oldPsig sig;    old value for pulse width if not defined in waveform file\n-missGround;     assume ground is missed to assess RH metrics\n-minGap gap;     delete signal beneath min detectable gap fraction\n-maxDN max;      maximum DN\n-bitRate n;      DN bit rate\n\nDenoising:\n-meanN n;        mean noise level\n-thresh n;       noise threshold\n-sWidth sig;     smoothing width\n-psWidth sigma;  pre-smoothing width\n-gWidth sig;     Gaussian paremter selection width\n-minGsig sig;    minimum Gaussian sigma to fit\n-minWidth n;     minimum feature width in bins\n-varNoise;       variable noise threshold\n-varScale x;     variable noise threshold scale\n-statsLen len;   length to calculate noise stats over\n-medNoise;       use median stats rather than mean\n-noiseTrack;     use noise tracking\n-rhoG rho;       ground reflectance\n-rhoC rho;       canopy reflectance\n-offset y;       waveform DN offset\n-pFile file;     read pulse file, for deconvoltuion and matched filters\n-preMatchF;      matched filter before denoising\n-postMatchF;     matched filter after denoising\n-gold;           deconvolve with Gold's method\n-deconTol;       deconvolution tolerance\n\nQuestions to svenhancock@gmail.com\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry gediRat -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }

  /*read deconvolution pulse if needed*/
  if(dimage->den->preMatchF||dimage->den->preMatchF||dimage->den->deconMeth>=0)readPulse(dimage->den); 

  return(dimage);
}/*readCommands*/

/*the end*/
/*###########################################################*/

