#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"



/*#################################*/
/*# Fit a Gaussian to transmitted #*/
/*# LVIS pulse to find centre     #*/
/*# 2016    svenhancock@gmail.com #*/
/*#################################*/

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



/*####################################*/
/*control structure*/

typedef struct{
  char inNamen[200];
  char outNamen[200];
  int meanBins;
  float oRes;
  int minN;
}control;


/*####################################*/
/*data structure*/

typedef struct{
  int nWaves;   /*number of waveforms*/
  float **wave; /*waveforms*/
  int nBins;   /*number of bins per waveform*/
}dataStruct;


/*############################################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct *data=NULL;
  dataStruct *readData(char *);
  float **fitPulseGauss(dataStruct *,int *,float,int);
  float **meanWaves=NULL;
  void writeResults(float **,int,float,char *);


  /*read command Line*/
  dimage=readCommands(argc,argv);

  /*read data*/
  data=readData(dimage->inNamen);

  /*perform fits*/
  meanWaves=fitPulseGauss(data,&dimage->meanBins,dimage->oRes,dimage->minN);

  /*write results*/
  writeResults(meanWaves,dimage->meanBins,dimage->oRes,dimage->outNamen);

  /*tidy up*/
  TTIDY((void **)meanWaves,2);
  if(data){
    TTIDY((void **)data->wave,data->nWaves);
    TIDY(data);
  }
  if(dimage){
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*############################################################*/
/*write results*/

void writeResults(float **meanWaves,int nBins,float res,char *outNamen)
{
  int i=0;
  FILE *opoo=NULL;

  if((opoo=fopen(outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",outNamen);
    exit(1);
  }

  fprintf(opoo,"# 1 x, 2 CofG, 3 Gaussian\n");
  for(i=0;i<nBins;i++)fprintf(opoo,"%f %f %f\n",(float)(i-nBins/2)*res,meanWaves[0][i],meanWaves[1][i]);

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",outNamen);
  return;
}/*writeResults*/


/*############################################################*/
/*fit Gaussian to pulse*/

float **fitPulseGauss(dataStruct *data,int *meanBins,float oRes,int minN)
{
  int i=0,numb=0,nGauss=0;
  int **nIn=NULL,bin=0;
  float *temp=NULL,*denoise=NULL;
  float *copyLastFeature(float *,int);
  float *fitSingleGauss(float *,float *,int,float,int *,float **);
  float *x=NULL,*fitWave=NULL;
  float *gaussPar=NULL,CofG=0;
  float total=0;
  float findCofG(float *,float *,int);
  float **meanWaves=NULL;
  denPar den;
  void setDenoiseDefault(denPar *);

  /*denoising parameters*/
  setDenoiseDefault(&den);
  den.varNoise=1;
  den.threshScale=5.0;
  den.noiseTrack=1;
  den.statsLen=3.0;

  (*meanBins)=(int)((float)data->nBins*den.res/oRes);
  meanWaves=fFalloc(2,"meanWaves",0);
  nIn=iIalloc(2,"nIn",0);
  for(i=0;i<2;i++){
    meanWaves[i]=falloc((*meanBins),"meanWaves",i+1);
    nIn[i]=ialloc((*meanBins),"nIn",i+1);
    for(numb=0;numb<(*meanBins);numb++){
      meanWaves[i][numb]=0.0;
      nIn[i][numb]=0;
    }
  }

  x=falloc(data->nBins,"x",0);
  for(i=0;i<data->nBins;i++)x[i]=(float)i*den.res;

  /*loop over waveforms*/
  for(numb=0;numb<data->nWaves;numb++){
    /*reverse waveform to ignore early reflectaion*/
    temp=falloc(data->nBins,"temp wave",0);
    for(i=0;i<data->nBins;i++)temp[i]=data->wave[numb][data->nBins-(i+1)];

    /*denoise*/
    denoise=processFloWave(temp,data->nBins,&den,1.0);
    TIDY(temp);

    /*keep last waveform*/
    temp=copyLastFeature(denoise,data->nBins);

    /*fit a single Gaussian*/
    fitWave=fitSingleGauss(x,temp,data->nBins,0.5,&nGauss,&gaussPar);
    CofG=findCofG(x,temp,data->nBins);

    /*load into mean arrays*/
    for(i=0;i<data->nBins;i++){
      bin=(int)((x[i]-CofG)/oRes+0.5)+(*meanBins)/2;
      if((bin>=0)&&(bin<(*meanBins))){
        meanWaves[0][bin]+=temp[i];
        nIn[0][bin]++;
      }
      bin=(int)((x[i]-gaussPar[0])/oRes+0.5)+(*meanBins)/2;
      if((bin>=0)&&(bin<(*meanBins))){
        meanWaves[1][bin]+=temp[i];
        nIn[1][bin]++;
      }
    }
    TIDY(temp);
    TIDY(denoise);
    TIDY(gaussPar);
    TIDY(fitWave);
    TIDY(gaussPar);
  }/*waveform loop*/

  /*normalise*/
  for(i=0;i<(*meanBins);i++){
    for(numb=0;numb<2;numb++){
      if(nIn[numb][i]>minN)meanWaves[numb][i]/=(float)nIn[numb][i];
      else                 meanWaves[numb][i]=0.0;
    }
  }

  /*remove outliers*/
  for(numb=0;numb<2;numb++){
    for(i=(*meanBins)/2;i<(*meanBins);i++){
      if(meanWaves[numb][i]<=0.0){
        for(;i<(*meanBins);i++)meanWaves[numb][i]=0.0;
      }
    }
    for(i=(*meanBins)/2;i>=0;i--){
      if(meanWaves[numb][i]<=0.0){
        for(;i>=0;i--)meanWaves[numb][i]=0.0;
      }
    }
  }

  /*normalise*/
  for(numb=0;numb<2;numb++){
    total=0.0;
    for(i=0;i<(*meanBins);i++)total+=meanWaves[numb][i];
    for(i=0;i<(*meanBins);i++)meanWaves[numb][i]/=total;
  }

  TIDY(x);
  TTIDY((void **)nIn,2);
  return(meanWaves);
}/*fitPulseGauss*/


/*############################################################*/
/*centre pf gravity*/

float findCofG(float *x,float *temp,int nBins)
{
  int i=0;
  float CofG=0,total=0;

  CofG=total=0.0;
  for(i=0;i<nBins;i++){
    CofG+=x[i]*temp[i];
    total+=temp[i];
  }
  if(total>0.0)CofG/=total;
  else         CofG=-9999.0;

  return(CofG);
}/*findCofG*/


/*############################################################*/
/*copy last feature to avoid reflection*/

float *copyLastFeature(float *wave,int nBins)
{
  int i=0,bin=0;
  int sBin=0,eBin=0;
  float *temp=NULL;
  float total=0,cumul=0;
  float thresh=0.0;

  /*total waveform*/
  total=0.0;
  for(i=0;i<nBins;i++)total+=wave[i];

  cumul=0.0;
  thresh=0.02*total;
  sBin=-1;
  for(i=0;i<nBins;i++){
    cumul+=wave[i];
    if((cumul>=thresh)&&(sBin<0))sBin=i;
    else if((sBin>=0)&&(wave[i]<=0.0)){
      eBin=i;
      break;
    }
  }

  temp=falloc(nBins,"temp",0);
  for(i=0;i<nBins;i++){
    bin=nBins-(i+1);
    if((i>=sBin)&&(i<=eBin))temp[bin]=wave[i];
    else                    temp[bin]=0.0;
  }

  return(temp);
}/*copyLastFeature*/


/*############################################################*/
/*read data*/

dataStruct *readData(char *namen)
{
  int i=0,j=0,bin=0;
  int sCol=0,eCol=0;
  dataStruct *data=NULL;
  char line[10000],*token=NULL;
  FILE *ipoo=NULL;

  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error metric structure allocation.\n");
    exit(1);
  }

  /*first count the number of waveforms*/
  data->nWaves=0;
  while(fgets(line,10000,ipoo)!=NULL){
    if(strncasecmp(line,"lfid",4))data->nWaves++;
    else{
      i=0;
      token=strtok(line," ");
      while(token){
        if(!strncasecmp(token,"tx00",4))sCol=i;
        else if(!strncasecmp(token,"rx000",5))eCol=i;
        token=strtok(NULL," ");
        i++;
      }
    }
  }/*count number of wavefoms*/


  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*allocate space*/
  data->nBins=eCol-sCol;
  data->wave=fFalloc(data->nWaves,"waveforms",0);
  for(i=0;i<data->nWaves;i++)data->wave[i]=falloc(data->nBins,"waveforms",i+1);

  /*read data*/
  j=0;
  while(fgets(line,10000,ipoo)!=NULL){
    if(strncasecmp(line,"lfid",4)){
      i=0;
      token=strtok(line," ");
      while(token){
        bin=i-sCol;
        if((bin>=0)&&(bin<data->nBins))data->wave[j][bin]=atof(token);
        token=strtok(NULL," ");
        i++;
      }
      j++;
    }
  }

  /*close file*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(data);
}/*readData*/


/*############################################################*/
/*read commands*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;

  /*allocate structures*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  strcpy(dimage->inNamen,"/Users/stevenhancock/data/teast/pulse/howland.waves");
  strcpy(dimage->outNamen,"teast.dat");
  dimage->oRes=0.15;
  dimage->minN=100;


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-res",4)){
        checkArguments(1,i,argc,"-res");
        dimage->oRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-minN",5)){
        checkArguments(1,i,argc,"-minN");
        dimage->minN=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to determine LVIS pulse shape\n#####\n\n-input name;   input filaname\n-output name;  output filename\n-res res;      output resolution\n-minN min;     minimum number of samples to trust\n\n");
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
/*############################################################*/

