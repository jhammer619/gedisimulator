#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "tiffWrite.h"


/*#########################*/
/*# Makes geotiffs from   #*/
/*# las files   2015      #*/
/*# svenhancock@gmail.com #*/
/*#########################*/

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


/*##################################################*/
/*control structure*/

typedef struct{
  /*input output*/
  char **inList;
  int nFiles;
  char outNamen[200];
  char bNamen[200];
  FILE *bFile;        /*bounds file output*/

  /*switches*/
  char drawInt;       /*intensity image switch*/
  char drawHeight;    /*height.elevation image switch*/
  char drawDens;      /*draw ensity images*/
  char findDens;      /*find point and footprint density*/
  char writeBounds;   /*write out file bounds*/
  uint64_t pBuffSize; /*point buffer rading size in bytes*/
  char printNpoint;   /*write number of points to the screen*/

  /*geotiff parts*/
  float res;
  float maxDN;
  uint16_t epsg;
}control;


/*##################################################*/
/*image structure*/

typedef struct{
  int nX;
  int nY;
  int *nIn;       /*point density*/
  int *nFoot;     /*footprint density*/
  float *jimlad;  /*FLOATING POINT IMAGE*/
  float min;      /*min intensity*/
  float max;      /*max intensity*/
  unsigned char *image;  /*IMAGE TO WRITE*/
  double minX;
  double minY;
  double maxX;
  double maxY;
  double geoL[2];
  int geoI[2];
  uint16_t epsg;   /*EPSG code*/
  int maxPoint;    /*max number of points per pixels*/
  int maxFoot;     /*max number of footprints per pixels*/
}imageStruct;


/*##################################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  lasFile **las=NULL;
  imageStruct *image=NULL;
  imageStruct *allocateImage(control *,lasFile **);
  imageStruct *tidyImage(imageStruct *);
  void collateImage(control *,lasFile **,imageStruct *);
  void writeImage(control *,imageStruct *);
  void writeFileBounds(lasFile *,char *,control *);


  /*read command line*/
  dimage=readCommands(argc,argv);

  /*allocate space*/
  if(!(las=(lasFile **)calloc(dimage->nFiles,sizeof(lasFile *)))){
    fprintf(stderr,"error contN allocation.\n");
    exit(1);
  }

  /*read file headers*/
  for(i=0;i<dimage->nFiles;i++){
    las[i]=readLasHead(dimage->inList[i],dimage->pBuffSize);
    if(dimage->epsg==0)readLasGeo(las[i]);
    else               las[i]->epsg=dimage->epsg;
    if(dimage->writeBounds)writeFileBounds(las[i],dimage->inList[i],dimage);
    if(dimage->printNpoint)fprintf(stdout,"nPoints %s %u\n",dimage->inList[i],las[i]->nPoints);
  }

  if(dimage->drawInt||dimage->drawHeight||dimage->findDens||dimage->drawDens){
    /*allocate image array*/
    image=allocateImage(dimage,las);

    /*read data and set up image*/
    collateImage(dimage,las,image);
  }
  if(las){
    for(i=0;i<dimage->nFiles;i++)las[i]=tidyLasFile(las[i]);
    TIDY(las);
  }

  /*write image*/
  if(dimage->drawInt||dimage->drawHeight)writeImage(dimage,image);
  if(dimage->writeBounds)fprintf(stdout,"Written to %s\n",dimage->bNamen);


  /*tidy up arrays*/
  image=tidyImage(image);
  if(dimage){
    TTIDY((void **)dimage->inList,dimage->nFiles);
    if(dimage->bFile){
      fclose(dimage->bFile);
      dimage->bFile=NULL;
    }
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*##################################################*/
/*tidy image structure*/

imageStruct *tidyImage(imageStruct *image)
{
  if(image){
    TIDY(image->jimlad);
    TIDY(image->nIn);
    TIDY(image->image);
    TIDY(image->nFoot);
    TIDY(image);
  }
  return(image);
}/*tidyImage*/


/*##################################################*/

void writeFileBounds(lasFile *las,char *namen,control *dimage)
{
  fprintf(dimage->bFile,"%s %.2f %.2f %.2f %.2f %.2f %.2f\n",namen,las->minB[0],las->minB[1],las->minB[2],las->maxB[0],las->maxB[1],las->maxB[2]);
  
  return;
}/*writeFileBounds*/


/*##################################################*/
/*write image to geotiff*/

void writeImage(control *dimage,imageStruct *image)
{
  drawTiff(dimage->outNamen,&(image->geoL[0]),&(image->geoI[0]),(double)dimage->res,image->image,image->nX,image->nY,1.0,image->epsg);

  return;
}/*writeImage*/


/*##################################################*/
/*collate image*/

void collateImage(control *dimage,lasFile **las,imageStruct *image)
{
  int i=0,place=0;
  int xBin=0,yBin=0;
  int nContP=0,nContF=0;
  uint32_t j=0;
  double x=0,y=0,z=0;
  float meanPoint=0,meanFoot=0;

  /*geolocation*/
  image->geoI[0]=image->geoI[1]=0;
  image->geoL[0]=image->minX+0.5*(double)dimage->res;
  image->geoL[1]=image->maxY-0.5*(double)dimage->res;
  image->epsg=las[0]->epsg;

  for(i=0;i<dimage->nFiles;i++){
    if(las[i]->epsg!=image->epsg){
      fprintf(stderr,"EPSG mismatch %d %d\n",(int)image->epsg,las[i]->epsg);
      exit(1);
    }
    for(j=0;j<las[i]->nPoints;j++){
      readLasPoint(las[i],j);
      setCoords(&x,&y,&z,las[i]);

      xBin=(int)((x-image->minX)/(double)dimage->res);
      yBin=(int)((image->maxY-y)/(double)dimage->res);

      if((xBin>=0)&&(xBin<image->nX)&&(yBin>=0)&&(yBin<image->nY)){
        place=yBin*image->nX+xBin;
        if(dimage->drawInt)image->jimlad[place]+=(float)las[i]->refl;
        else if(dimage->drawHeight)image->jimlad[place]+=(float)z;
        if(dimage->findDens&&(las[i]->field.retNumb==las[i]->field.nRet))image->nFoot[place]++;
        image->nIn[place]++;
      }else{
        fprintf(stderr,"How can we be outside %d %d? x %f %f %f y %f %f %f\n",xBin,yBin,image->minX,x,image->maxX,image->minY,y,image->maxY);
      }
    }/*point loop*/
  }/*file loop*/

  /*normalise and find bounds*/
  meanPoint=meanFoot=0.0;
  nContP=nContF=0;
  for(i=image->nX*image->nY-1;i>=0;i--){
    if(image->nIn[i]>0){
      if(dimage->drawInt||dimage->drawHeight){
        image->jimlad[i]/=(float)image->nIn[i];
        if(image->jimlad[i]<image->min)image->min=image->jimlad[i];
        if(image->jimlad[i]>image->max)image->max=image->jimlad[i];
      }
      if(dimage->findDens){
        if(image->nIn[i]>image->maxPoint)image->maxPoint=image->nIn[i];
        if(image->nFoot[i]>image->maxFoot)image->maxFoot=image->nFoot[i];
        if(image->nIn[i]>0){
          meanPoint+=(float)image->nIn[i];
          nContP++;
        }
        if(image->nFoot[i]>0){
          meanFoot+=(float)image->nFoot[i];
          nContF++;
        }
      }
    }
  }

  if(!dimage->drawDens)TIDY(image->nIn);
  if(dimage->findDens){
    if(nContF>0)meanFoot/=(float)nContF;
    if(nContP>0)meanPoint/=(float)nContP;
    fprintf(stdout,"Mean point density %f per m2\n",meanPoint/(dimage->res*dimage->res));
    fprintf(stdout,"Mean footprint density %f per m2\n",meanFoot/(dimage->res*dimage->res));
  }


  if(dimage->maxDN>0.0){
    if(image->max>dimage->maxDN)dimage->maxDN=image->max;
  }

  /*copy to uchar array*/
  if(dimage->drawInt||dimage->drawHeight){
    image->image=uchalloc(image->nX*image->nY,"image",0);
    for(i=image->nX*image->nY-1;i>=0;i--){
      if(image->jimlad[i]<image->max)image->image[i]=(unsigned char)((image->jimlad[i]-image->min)*255.0/(image->max-image->min));
      else                           image->image[i]=255;
    }
  }

  TIDY(image->jimlad);
  return;
}/*collateImage*/


/*##################################################*/
/*allocate image struture*/

imageStruct *allocateImage(control *dimage,lasFile **las)
{
  int i=0;
  /*uint32_t j=0;
  double x=0,y=0,z=0;*/
  imageStruct *image=NULL;

  if(!(image=(imageStruct *)calloc(1,sizeof(imageStruct)))){
    fprintf(stderr,"error imageStruct allocation.\n");
    exit(1);
  }

  image->maxX=image->maxY=-10000000000.0;
  image->minX=image->minY=10000000000.0;

  for(i=0;i<dimage->nFiles;i++){
    if(las[i]->minB[0]<image->minX)image->minX=las[i]->minB[0];
    if(las[i]->minB[1]<image->minY)image->minY=las[i]->minB[1];
    if(las[i]->maxB[0]>image->maxX)image->maxX=las[i]->maxB[0];
    if(las[i]->maxB[1]>image->maxY)image->maxY=las[i]->maxB[1];

    /*for(j=0;j<las[i]->nPoints;j++){
      readLasPoint(las[i],j);
      setCoords(&x,&y,&z,las[i]);
      if(x<image->minX)image->minX=x;
      if(y<image->minY)image->minY=y;
      if(x>image->maxX)image->maxX=x;
      if(y>image->maxY)image->maxY=y;
    }*/
  }

  image->nX=(int)((image->maxX-image->minX)/(double)dimage->res)+1;
  image->nY=(int)((image->maxY-image->minY)/(double)dimage->res)+1;
  fprintf(stdout,"Image will be %d by %d\n",image->nX,image->nY);

  if(dimage->drawInt||dimage->drawHeight)image->jimlad=falloc(image->nX*image->nY,"jimlad",0);
  else                                   image->jimlad=NULL;
  image->nIn=ialloc(image->nX*image->nY,"nIn",0);
  if(dimage->findDens)image->nFoot=ialloc(image->nX*image->nY,"nFoot",0);
  else                image->nFoot=NULL;
  for(i=image->nX*image->nY-1;i>=0;i--){
    if(dimage->drawInt||dimage->drawHeight)image->jimlad[i]=0.0;
    if(dimage->findDens)image->nFoot[i]=0;
    image->nIn[i]=0;
  }
  image->min=1000000.0;
  image->max=-1000000.0;
  image->maxFoot=image->maxPoint=0;

  return(image);
}/*allocateImage*/


/*##################################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  char **readInList(int *,char *);

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error contN allocation.\n");
    exit(1);
  }

  dimage->nFiles=1;
  dimage->inList=chChalloc(dimage->nFiles,"inList",0);
  dimage->inList[0]=challoc(200,"inList",0);
  strcpy(&(dimage->inList[0][0]),"/Users/stevenhancock/data/gedi/USDA_pilot/waveform/Lift02/WF_V13_-Riegl680i-HRZ-140701_130919_1-originalpoints.las");
  strcpy(dimage->outNamen,"teast.tif");
  dimage->drawInt=1;
  dimage->drawHeight=0;
  dimage->drawDens=0;
  dimage->findDens=0;
  dimage->writeBounds=0;
  dimage->printNpoint=0;
  dimage->bFile=NULL;
  dimage->epsg=0;    /*leave bank*/
  dimage->pBuffSize=(uint64_t)200000000;

  dimage->res=100.0;
  dimage->maxDN=-1.0;


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
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=readInList(&dimage->nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-res",4)){
        checkArguments(1,i,argc,"-res");
        dimage->res=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxDN",6)){
        checkArguments(1,i,argc,"-maxDN");
        dimage->maxDN=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-epsg",5)){
        checkArguments(1,i,argc,"-epsg");
        dimage->epsg=(uint16_t)atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noInt",6)){
        dimage->drawInt=0;
      }else if(!strncasecmp(argv[i],"-height",7)){
        dimage->drawHeight=1;
        dimage->drawInt=0;
      }else if(!strncasecmp(argv[i],"-findDens",9)){
        dimage->findDens=1;
      }else if(!strncasecmp(argv[i],"-writeBound",11)){
        checkArguments(1,i,argc,"-writeBound");
        dimage->writeBounds=1;
        strcpy(dimage->bNamen,argv[++i]);
        if((dimage->bFile=fopen(dimage->bNamen,"w"))==NULL){
          fprintf(stderr,"Error opening output file %s\n",dimage->bNamen);
          exit(1);
        }
      }else if(!strncasecmp(argv[i],"-pBuff",6)){
        checkArguments(1,i,argc,"-pBuff");
        dimage->pBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0);
      }else if(!strncasecmp(argv[i],"-printNpoint",12)){
        dimage->printNpoint=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to create GEDI waveforms from ALS las files\n#####\n\n-input name;     lasfile input filename\n-output name;    output filename\n-inList list;    input file list for multiple files\n-res res;        image resolution, in metres\n-height;     draw height image\n-noInt;          no image\n-findDens;       find point and footprint density\n-epsg n;         geolocation code if not read from file\n-writeBound n;   write file bounds to a file\n-pBuff s;        point reading buffer size in Gbytes\n-printNpoint;    print number of points in each file\n\nQuestions to svenhancock@gmail.com\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry gediRat -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }/*command parser*/

  return(dimage);
}/*readCommands*/

/*the end*/
/*##################################################*/

