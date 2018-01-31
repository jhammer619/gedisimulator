


/*##############################*/
/*# Header for gediNoise.c    #*/
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


/*tolerances*/
#define XRES 0.0000025
#define YTOL 0.000000001   /*for determining Gaussian thresholds*/


/*#######################################*/
/*noise structure*/

typedef struct{
  /*link margin defined noise*/
  char linkNoise;  /*use link noise or not*/
  float linkCov;   /*cover at which link margin is defined*/
  float linkM;     /*link margin*/
  float linkSig;   /*link noise sigma*/
  /*general noise statistics*/
  float meanN;     /*mean noise offset*/
  float trueSig;   /*true noise sigma in DN*/
  float nSig;      /*noise sigma*/
  float offset;   /*waveform DN offset*/
  char bitRate;   /*digitiser bit rate*/
  float maxDN;    /*maximum DN we need to digitise*/
  /*others*/
  char missGround; /*force to miss ground to get RH errors*/
  float deSig;     /*detector sigma*/
}noisePar;


/*the end*/
/*#######################################*/

