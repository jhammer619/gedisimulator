#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"

/*##############################*/
/*# Generates GEDI trakcs      #*/
/*# from ornital simulations   #*/
/*# 2017 svenhancock@gmail.com #*/
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


/*##############################################*/
/*control structure*/

typedef struct{
  char trackNamen[200];
  char metricNamen[200];
  char outNamen[200];

  /*switches*/
  char readMetric;   /*read metric or set coordinate switch*/

  /*GEDI settings*/
  float res;         /*grid resolution*/
  float alongTrack;  /*along track footprint spacing*/
  float acrossTrack; /*across track footprint spacing*/

  /*others*/
}dimage;


/*##############################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);

  /*read command Line*/
  dimage=readCommands(argc,argv);

  /*read metric file for coords and bounds*/

  /*read track file*/

  /*choose coords to use*/

  /*output relevant lines*/

  /*tidy up*/
  if(dimage){
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*##############################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;


  /*allocate structures*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*defaults*/
  strcpy(dimage->trackNamen,"/Users/stevenhancock/data/teast/tracks/GEDI_grid1000m_dy1_395days_lat0N_parsed.csv");
  strcpy(simage->metricNamen,"");
  strcpy(dimage->outNamen,"teast.dat");

  /*switches*/
  dimage->readMetric=1;   /*read metric or set coordinate switch*/

  /*GEDI settings*/
  dimage->res=1000.0;         /*grid resolution*/
  dimage->alongTrack=60.0;    /*along track footprint spacing*/
  dimage->acrossTrack=600.0;  /*across track footprint spacing*/

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-trackFile",10)){
        checkArguments(1,i,argc,"-trackFile");
        strcpy(dimage->trackNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-metricFile",11)){
        checkArguments(1,i,argc,"-metricFile");
        strcpy(dimage->metricNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-seed",5)){
        checkArguments(1,i,argc,"-seed");
        srand(atoi(argv[++i]));
      }else if(!strncasecmp(argv[i],"-cloudFrac",10)){
        checkArguments(1,i,argc,"-cloudFrac");
        dimage->cloudFrac=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to calculate GEDI waveform metrics\n#####\n\n");
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
/*##############################################*/

