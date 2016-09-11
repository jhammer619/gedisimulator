#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "tools.h"
#include "tools.c"

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


#define TOL 0.00000001

/*#############################*/
/*main*/

int main()
{
  int nBins;
  float res=0;
  float pSig=0,fSig=0;
  float A=0,slope=0,fSlope=0;
  float maxSlope=0,step=0;
  float *wave=NULL;
  float *makeWave(float,float,float,float,int *,float);
  float stdev=0;
  float waveStdev(float *,int,float);
  float findSlope(float,float,float);

  /*instrument parameters*/
  res=0.15;
  A=0.5;
  pSig=0.764331;
  fSig=11.0;

  /*slope range*/
  maxSlope=60.0*M_PI/180.0;
  step=0.5*M_PI/180.0;
  slope=0.0;

  /*loop over slopes*/
  while(slope<=maxSlope){
    /*make the waveform*/
    wave=makeWave(slope,pSig,fSig,A,&nBins,res);

    /*calculate stdev*/
    stdev=waveStdev(wave,nBins,res);
    TIDY(wave);

    /*calculate slope*/
    fSlope=findSlope(stdev,pSig,fSig);

    fprintf(stdout,"%f %f %f %f %f\n",slope*180.0/M_PI,stdev,fSlope,pSig,fSig);
    slope+=step;
  }

  return(0);
}/*main*/


/*#############################*/
/*slope from stdev*/

float findSlope(float stdev,float pSig,float fSig)
{
  float slope=0;

  slope=atan2(sqrt(stdev*stdev-pSig*pSig),fSig)*180.0/M_PI;

  return(slope);
}/*findSlope*/


/*#############################*/
/*determine standard deviation*/

float waveStdev(float *wave,int nBins,float res)
{
  int i=0;
  float x=0,tot=0;
  float mean=0,stdev=0;

  mean=tot=0.0;
  for(i=0;i<nBins;i++){
    x=(float)i*res;
    mean+=wave[i]*x;
    tot+=wave[i];
  }
  mean/=tot;

  stdev=0.0;
  for(i=0;i<nBins;i++){
    x=(float)i*res;
    stdev+=(x-mean)*(x-mean)*wave[i];
  }
  stdev=sqrt(stdev/tot);

  return(stdev);
}/*waveStdev*/


/*#############################*/
/*make waveform*/

float *makeWave(float slope,float pSig,float fSig,float A,int *nBins,float res)
{
  int bin=0;
  float *wave=NULL,aRes=0;
  float footRad=0;
  float x=0,y=0,z=0;
  float tot=0,*smoothed=NULL;;
  float *smooth(float,int,float *,float);
  double gaussian(double,double,double);

  *nBins=2000;
  wave=falloc(*nBins,"wave",0);

  /*determine footprint radius*/
  footRad=0.0;
  aRes=0.05;
  do{
    y=A*(float)gaussian((double)footRad,(double)fSig,0.0);
    footRad+=aRes;
  }while(y>=TOL);

  /*calculate floor return*/
  x=-1.0*footRad;
  do{
    z=x*sin(slope)/cos(slope);
    bin=(int)(z/res)+(*nBins)/2;
    if((bin<0)||(bin>=(*nBins))){
      fprintf(stderr,"Not enough bins\n");
      exit(1);
    }

    tot=0.0;
    z=0.0;
    do{
      y=(float)gaussian((double)z,(double)fSig,0.0);
      tot+=y*aRes;
      z+=aRes;
    }while(y>=TOL);
    tot*=2.0;
    wave[bin]+=tot*(float)gaussian((double)x,(double)fSig,0.0)*aRes;
    x+=aRes;
  }while(x<=footRad);

  /*smooth by pulse blurring*/
  smoothed=smooth(pSig,*nBins,wave,res);
  TIDY(wave);

  return(smoothed);
}/*makeWave*/

/*the end*/
/*#############################*/

