#include "stdio.h"
#include "stdlib.h"
#include "string.h"
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

/*#######################################*/
/*control structure*/

typedef struct{
  char inNamen[200];   /*input filename*/
  char outNamen[200];  /*output filename*/
  /*input settings*/
  char readHDFgedi;    /*read GEDI HDF*/
  /*noise*/
  noisePar noise;  /*noise parameter structure*/
  float linkFsig;  /*footprint sigma used for link margin*/
  float linkPsig;  /*pulse sigma used for link margin*/
  /*parameters*/
  gediIOstruct gediIO;
}control;


/*#######################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  gediHDF *hdfGedi;     /*GEDI HDF5 structure*/
  gediHDF *noiseGEDI(control *);

  /*read command Line*/
  dimage=readCommands(argc,argv);

  /*read data and add noise*/
  hdfGedi=noiseGEDI(dimage);

  /*write data*/
  writeGEDIhdf(hdfGedi,dimage->outNamen,&dimage->gediIO);

  /*tidy up*/
  hdfGedi=tidyGediHDF(hdfGedi);
  if(dimage){
    TIDY(dimage);
    dimage=NULL;
  }
  return(0);
}/*main*/


/*####################################################*/
/*add noise to waveforms*/

gediHDF *noiseGEDI(control *dimage)
{
  int i=0,j=0;
  float rhoC=0,rhoG=0;
  float rhoRatio=0;
  gediHDF *hdfGedi=NULL;
  gediHDF *allHDF=NULL;
  dataStruct *data=NULL;
  dataStruct *tidyData(dataStruct *,char);

  /*determine noise level*/
  rhoG=0.4;
  rhoC=0.57;
  rhoRatio=rhoC/rhoG;
  dimage->noise.linkSig=setNoiseSigma(dimage->noise.linkM,dimage->noise.linkCov,dimage->linkPsig,dimage->linkFsig,rhoC,rhoG);

  /*read dummy hdf data*/
  dimage->gediIO.useInt=dimage->gediIO.useCount=dimage->gediIO.useFrac=0;
  data=unpackHDFgedi(dimage->inNamen,&dimage->gediIO,&allHDF,0);
  data=tidyData(data,dimage->readHDFgedi);

  /*loop over wave types*/
  for(j=0;j<allHDF->nTypeWaves;j++){
    if(allHDF->nTypeWaves==3){
      dimage->gediIO.useInt=(j==0)?1:0;
      dimage->gediIO.useCount=(j==1)?1:0;
      dimage->gediIO.useFrac=(j==2)?1:0;
    }else{
      dimage->gediIO.useCount=1;
      dimage->gediIO.useInt=dimage->gediIO.useFrac=0;
    }
    /*loop over waves*/
    for(i=0;i<allHDF->nWaves;i++){
      data=unpackHDFgedi(dimage->inNamen,&dimage->gediIO,&hdfGedi,i);
      /*determine true cover to use for link margin*/
      data->cov=waveformTrueCover(data,&dimage->gediIO,rhoRatio);
      /*add noise*/
      addNoise(data,&dimage->noise,dimage->linkFsig,dimage->linkPsig,data->res,rhoC,rhoG);
      /*copy back*/
      memcpy(&allHDF->wave[j][i*hdfGedi->nBins[0]],data->noised,data->nBins*sizeof(float));
      data=tidyData(data,dimage->readHDFgedi);
    }
    hdfGedi=tidyGediHDF(hdfGedi);
  }/*wave loop*/

  /*reset all wave flags*/
  if(allHDF->nTypeWaves==3){
    dimage->gediIO.useInt=1;
    dimage->gediIO.useCount=1;
    dimage->gediIO.useFrac=1;
  }else{
    dimage->gediIO.useCount=1;
    dimage->gediIO.useInt=dimage->gediIO.useFrac=0;
  }

  return(allHDF);
}/*noiseGEDI*/


/*####################################################*/
/*tidy data structure*/

dataStruct *tidyData(dataStruct *data,char readHDFgedi)
{
  if(data){
    TIDY(data->noised);
    if(readHDFgedi){  /*pointer to array. do not free*/
      data->wave[0]=NULL;
      if(data->ground)data->ground[0]=NULL;
    }
    TTIDY((void **)data->ground,data->nWaveTypes);
    TTIDY((void **)data->wave,data->nWaveTypes);
    TIDY(data->totE);
    TIDY(data->z);
    TIDY(data);
  }

  return(data);
}/*tidyData*/


/*####################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;

  /*allocate structures*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    errorf("error control allocation.\n");
    exit(1);
  }

  /*defaults*/
  dimage->readHDFgedi=1;
  dimage->noise.linkNoise=1;
  dimage->noise.shotNoise=0;
  dimage->noise.missGround=0;
  dimage->noise.trueSig=5.0;
  dimage->noise.maxDN=4096.0;
  dimage->noise.offset=94.0;
  dimage->noise.deSig=0.0;
  dimage->noise.bitRate=12;
  dimage->noise.linkCov=0.95;
  dimage->noise.linkM=3.0;
  dimage->noise.nSig=0.0;
  dimage->noise.meanN=0.0;
  dimage->noise.hNoise=0.0;
  dimage->gediIO.useInt=1;
  dimage->gediIO.useCount=1;
  dimage->gediIO.useFrac=1;
  dimage->gediIO.nMessages=200;
  dimage->gediIO.ground=1;
  dimage->linkFsig=5.5;
  dimage->linkPsig=0.9;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-linkFsig",9)){
        checkArguments(1,i,argc,"-linkFsig");
        dimage->linkFsig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-linkPsig",9)){
        checkArguments(1,i,argc,"-linkPsig");
        dimage->linkPsig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-linkNoise",10)){
        checkArguments(2,i,argc,"-linkNoise");
        dimage->noise.linkNoise=1;
        dimage->noise.linkM=atof(argv[++i]);
        dimage->noise.linkCov=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-trueSig",8)){
        checkArguments(1,i,argc,"-trueSig");
        dimage->noise.trueSig=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-bitRate",8)){
        checkArguments(1,i,argc,"-bitRate");
        dimage->noise.bitRate=(char)atoi(argv[++i]);
        dimage->noise.maxDN=pow(2.0,(float)dimage->noise.bitRate);
      }else if(!strncasecmp(argv[i],"-seed",5)){
        checkArguments(1,i,argc,"-seed");
        srand2(atoi(argv[++i]));
      }else if(!strncasecmp(argv[i],"-dcBias",7)){
        checkArguments(1,i,argc,"-dcBias");
        dimage->noise.offset=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        msgf("\n#####\nProgram to calculate GEDI waveform metrics\n#####\n\n-input name;     waveform  input filename\n-output name;   output filename\n-seed n;         random number seed\n-linkNoise linkM cov;     apply Gaussian noise based on link margin at a cover\n-dcBias dn;      mean noise level\n-linkFsig sig;       footprint width to use when calculating and applying signal noise\n-linkPsig sig;       pulse width to use when calculating and applying signal noise\n-trueSig sig;    true sigma of background noise\n-bitRate n;      DN bit rate\n\n");
        exit(1);
      }else{
        errorf("%s: unknown argument on command line: %s\nTry gediRat -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }

  return(dimage);
}/*readCommands*/


/*the end*/
/*#######################################*/

