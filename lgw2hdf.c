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


/*############################*/
/*# Converts LVIS lgw format #*/
/*# format into HDF5    2017 #*/
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
  char inNamen[200];   /*input filename*/
  char outNamen[200];  /*output filename*/
  gediIOstruct gediIO; /*input/output structure*/
  lvisLGWstruct lvis;   /*LVIS lgw structure*/
}control;


/*###########################################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct *data=NULL;

  /*read command Line*/
  dimage=readCommands(argc,argv);

  /*read LVIS data*/
  data=readBinaryLVIS(dimage->inNamen,&dimage->lvis,i,&dimage->gediIO);


  /*tidy up*/
  if(dimage){
    TIDY(dimage->gediIO.den);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*###########################################################*/
/*read command Line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;


  /*allocate structures*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  if(!(dimage->gediIO.den=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }


  /*defaults*/
  strcpy(dimage->outNamen,"teast.h5");
  /*set default denoising parameters*/
  setDenoiseDefault(dimage->gediIO.den);
  dimage->gediIO.den->meanN=0.0;  /*we haven't added noise yet*/
  dimage->gediIO.den->thresh=0.0;  /*tiny number as no noise yet*/
  dimage->gediIO.den->noiseTrack=0;
  dimage->gediIO.den->minWidth=0;
  dimage->gediIO.den->varNoise=0;
  dimage->gediIO.den->threshScale=1.5;
  dimage->gediIO.den->fitGauss=0;
  dimage->gediIO.den->psWidth=0.0;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to convert LVIS lgw files to HDF5\n#####\n\n-input name;     input LGW filename\n-output name;   output HDF5 filename\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry lgw2hdf -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }

  return(dimage);
}/*readCommands*/


/*the end*/
/*###########################################################*/

