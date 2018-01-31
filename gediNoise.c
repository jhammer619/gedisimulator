#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
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
  void gaussThresholds(double,double,double,double,double *,double *);
  char direction=0;

  /*initial guess*/
  sig=groundAmp/2.0;
  step=groundAmp/10.0;

  /*determine threshold in terms of standard deviations*/
  gaussThresholds(1.0,XRES,(double)probNoise,(double)probMiss,&nNsig,&nSsig);

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

void gaussThresholds(double sig,double res,double probNoise,double probMiss,double *threshN,double *threshS)
{
  double x=0,y=0;
  double cumul=0,tot=0;
  char foundS=0,foundN=0;
  double probN=0,probS=0;

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

/*the end*/
/*####################################################*/

