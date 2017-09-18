#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "hdf5.h"
#include "libLasRead.h"
#include "libLidVoxel.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "gediIO.h"


/*#########################*/
/*# Generated grids of    #*/
/*# waveforms    2015     #*/
/*# svenhancock@gmail.com #*/
/*#########################*/

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

#define TOL 0.000001   /*a tolerance*/


/*####################################*/
/*pulse structure*/

typedef struct{
  int nBins;
  int centBin;  /*peak bin*/
  float *y;
  float *x;
}pulseStruct;


/*####################################*/
/*lobe structure*/

typedef struct{
  double coord[2];  /*central coordinate*/
  float E;          /*fraction of energy here*/
  float fSigma;     /*footprint sigma*/
  double maxSepSq;  /*maximum distance from footprint needed*/
}lobeStruct;


/*####################################*/

typedef struct{
  char **inList;
  int nFiles;
  char outNamen[200];
  char waveNamen[200];
  char pulseFile[200];

  /*IO structure*/
  gediIOstruct gediIO; /*generic IO options*/

  /*options*/
  char ground;         /*only use ground points*/
  char readWave;       /*read waveform switch*/
  char listFiles;      /*list waves only*/
  char checkCover;     /*check that the whole footprit is covered by data*/
  char useFootprint;   /*use footprint or not flag*/
  char overWrite;      /*overwrite old waveform switch*/
  char cleanOut;       /*clean subterranean outliers*/
  char normCover;      /*normalise for variable ALS coverage*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/
  char waveID[200];    /*wave ID if we are to use it*/
  char useID;          /*use wave ID*/
  char readPulse;      /*read pulse to simulate with*/
  char useShadow;      /*account for shadowing through voxelisation*/
  float maxScanAng;    /*maximum scan angle*/
  char polyGr;         /*fit a polynomial to the ground*/
  char nnGr;           /*ground DEM from nearest neighbour*/

  /*HDF5 output*/
  char writeHDF;       /*write output as hdf5*/
  int maxBins;         /*bins per wave for HDF5 output*/
  int hdfCount;        /*count used footprints*/

  /*read a batch of coords*/
  char readALSonce;    /*read all ALS data once*/
  char coordList[200]; /*list of coordinates*/
  double **coords;     /*list of coordinates*/
  char **waveIDlist;   /*list of waveform IDs*/

  /*grid to output multiple waveforms per run*/
  char doGrid;         /*gridded switch*/
  double gRes;         /*grid resolution*/
  double gMinX;        /*minimum x of grid*/
  double gMaxX;        /*maximum x of grid*/
  double gMinY;        /*minimum y of grid*/
  double gMaxY;        /*maximum y of grid*/
  int gNx;             /*number of x steps*/
  int gNy;             /*number of y steps*/

  /*GEDI footprint parameters, per footprint*/
  char sideLobe;     /*side lobe switch*/
  float lobeAng;     /*lobe major axis, degrees*/
  int nLobes;        /*number of side lobes*/
  lobeStruct *lobe;  /*lobe structure*/
  double minX;       /*minimum latitude of interest*/
  double maxX;       /*maximum latitude of interest*/
  double minY;       /*minimum longitude of interest*/
  double maxY;       /*maximum longitude of interest*/
  char topHat;       /*use a top hat wavefront rather than Gaussian*/
  double coord[2];

  /*global area of interest*/
  double globMinX;
  double globMaxX;
  double globMinY;
  double globMaxY;
  
  /*GEDI characteristics*/
  float res;       /*range resolution*/
  float pFWHM;     /*pulse width in ns*/
  float pSigma;    /*pulse width in metres*/
  float fWidth;    /*footprint width*/
  float fSigma;    
  float pRes;
  pulseStruct *pulse;
  double maxSep;   /*maximum acceptable separation*/

  /*tolerances*/
  float iThresh;   /*intensity threshold*/
  float meanN;
  char doDecon;    /*deconolution switch*/
  char indDecon;   /*deconolve individual ALS waves*/
  denPar *decon;   /*denoising parameters*/

  /*grid to normalise sampling density*/
  int *nGrid;    /*beam per grid cell*/
  int gX;
  int gY;
  float gridRes;
  double g0[2];  /*grid origin*/

  /*point and beam density*/
  float denseRadSq; /*radius to calculate density within*/
  float pointDense; /*point density within 2 sigma*/
  float beamDense;  /*beam density within 2 sigma*/

  /*for voxel shadows*/
  float vRes[3];   /*resolution along each axis*/
  float beamRad;   /*beam radius at ground*/
}control;


/*####################################*/
/*waveform structure*/

typedef struct{
  float **wave;  /*waveforms*/
  float **canopy;/*canopy waveform*/
  float **ground;/*ground waveform*/
  double gElev;   /*ground elevation if calculated*/
  float gSlope;    /*ground sope if calculated*/
  double gElevSimp; /*simple ground elevation if calculated*/
  float gSlopeSimp;  /*simple ground sope if calculated*/
  float meanScanAng;/*mean ALS scan angle*/
  double minZ;    /*elevation bounds*/
  double maxZ;   /*elevation bounds*/
  int nBins;     /*number of wave bins*/
  int nWaves;    /*number of different ways*/
  double groundBreakElev;  /*break in ground*/
}waveStruct;


/*####################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0,j=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  lasFile *las=NULL;
  pCloudStruct **data=NULL;
  pCloudStruct *readALSdata(lasFile *,control *);
  waveStruct *waves=NULL;
  waveStruct *makeGediWaves(control *,pCloudStruct **);
  gediHDF *hdfData=NULL;
  gediHDF *setUpHDF(control *);
  void writeGEDIwave(control *,waveStruct *,int);
  void setGediFootprint(control *);
  void setGediPulse(control *);
  void setGediGrid(control *);
  void tidySMoothPulse();
  void checkThisFile(lasFile *,control *,int);
  void groundFromDEM(pCloudStruct **,control *,waveStruct *);
  void updateCoord(control *,int,int);
  void checkWaveOverwrite(control *,int);
  void packGEDIhdf(control *,waveStruct *,gediHDF *);

 
  /*read command line*/
  dimage=readCommands(argc,argv);

  /*set up the pulse*/
  setGediPulse(dimage);

  /*set up grid or batch if needed*/
  setGediGrid(dimage);

  /*loop over las files and read*/
  if(!(data=(pCloudStruct **)calloc(dimage->nFiles,sizeof(pCloudStruct *)))){
    fprintf(stderr,"error waveStruct allocation.\n");
    exit(1);
  }
  for(i=0;i<dimage->nFiles;i++){
    /*report progress if rading all data here*/
    if(dimage->doGrid||dimage->readALSonce)fprintf(stdout,"File %d of %d",i+1,dimage->nFiles);
    /*read lasFile*/
    las=readLasHead(dimage->inList[i],dimage->pBuffSize);

    /*read data or write filename if needed*/
    if(dimage->listFiles==0)data[i]=readALSdata(las,dimage);
    else                    checkThisFile(las,dimage,i);
    if(dimage->doGrid||dimage->readALSonce)fprintf(stdout," nPoints %u\n",data[i]->nPoints);

    /*tidy lasFIle*/
    las=tidyLasFile(las);
  }/*file loop*/

  /*set up HDF5 if needed*/
  if(dimage->writeHDF)hdfData=setUpHDF(dimage);

  /*make waveforms*/
  if(dimage->listFiles==0){
    /*loop over waveforms*/
    for(i=0;i<dimage->gNx;i++){
      for(j=0;j<dimage->gNy;j++){
        if(dimage->writeHDF)fprintf(stdout,"Wave %d of %d\n",i*dimage->gNy+j,dimage->gNx*dimage->gNy);
        /*update centre coord*/
        updateCoord(dimage,i,j);

        /*see if that file already exists*/
        checkWaveOverwrite(dimage,i);

        /*if it is not to be overwritten*/
        if(dimage->useFootprint){
          /*set up footprint*/
          setGediFootprint(dimage);

          /*make waveforms*/
          waves=makeGediWaves(dimage,data);
        }

        /*if it is usable*/
        if(dimage->useFootprint){
          /*find the ground if needed*/
          if(dimage->ground&&(dimage->polyGr||dimage->nnGr))groundFromDEM(data,dimage,waves);
  
          /*output results*/
          if(dimage->writeHDF)packGEDIhdf(dimage,waves,hdfData);
          else                writeGEDIwave(dimage,waves,i);
        }

        /*tidy up*/
        TIDY(dimage->nGrid);
        TIDY(dimage->lobe);
        if(waves){
          TTIDY((void **)waves->wave,waves->nWaves);
          TTIDY((void **)waves->canopy,3);
          TTIDY((void **)waves->ground,3);
          TIDY(waves);
        }
      }/*grid y loop*/
    }/*grid x loop*/

    /*write HDF if needed*/
    if(dimage->writeHDF){
      hdfData->nWaves=dimage->hdfCount;  /*account for unusable footprints*/
      writeGEDIhdf(hdfData,dimage->outNamen);
    }
  }/*make and write a waveform if needed*/


  /*tidy up*/
  if(data){
    for(i=0;i<dimage->nFiles;i++)data[i]=tidyPointCloud(data[i]);
    TIDY(data);
  }
  hdfData=tidyGediHDF(hdfData);
  if(dimage){
    if(dimage->pulse){
      TIDY(dimage->pulse->y);
      TIDY(dimage->pulse->x);
      TIDY(dimage->pulse);
    }
    TTIDY((void **)dimage->coords,dimage->gNx);
    TTIDY((void **)dimage->waveIDlist,dimage->gNx);
    TTIDY((void **)dimage->inList,dimage->nFiles);
    TIDY(dimage);
  }
  tidySMoothPulse();
  return(0);
}/*main*/


/*##############################################*/
/*copy waveform into HDF structure*/

void packGEDIhdf(control *dimage,waveStruct *waves,gediHDF *hdfData)
{
  int i=0,j=0,start=0,numb=0;
  int nBins=0,idLength=0;
  float *tot=NULL,*cumul=NULL;
  float *thresh=NULL,buff=0;
  char waveID[100];

  numb=dimage->hdfCount;

  /*trim waveform*/
  buff=30.0;
  /*find energies*/
  tot=falloc(hdfData->nTypeWaves,"tot",0);
  cumul=falloc(hdfData->nTypeWaves,"cumul",0);
  for(j=0;j<hdfData->nTypeWaves;j++){
    tot[j]=cumul[j]=0.0;
    for(i=0;i<waves->nBins;i++)tot[j]+=waves->wave[j][i];
  }
  /*set threshols*/
  thresh=falloc(hdfData->nTypeWaves,"thresh",0);
  for(j=0;j<hdfData->nTypeWaves;j++)thresh[j]=0.01*tot[j];
  TIDY(tot);
  /*find waveform start*/
  start=-1;
  for(i=0;i<waves->nBins;i++){
    for(j=0;j<hdfData->nTypeWaves;j++){
      cumul[j]+=waves->wave[j][i];
      if(cumul[j]>thresh[j]){
        start=i;
        break;
      }
    }
    if(start>=0)break;
  }/*waveform trimming*/
  TIDY(cumul);
  TIDY(thresh);

  start-=buff/dimage->res;
  if(start<0)start=0;

  /*copy data*/
  hdfData->z0[numb]=waves->maxZ-(float)start*dimage->res;
  hdfData->zN[numb]=hdfData->z0[numb]-(float)hdfData->nBins*dimage->res;
  hdfData->lon[numb]=dimage->coord[0];
  hdfData->lat[numb]=dimage->coord[1];
  hdfData->beamDense[numb]=dimage->beamDense;
  hdfData->pointDense[numb]=dimage->pointDense;
  hdfData->zen[numb]=waves->meanScanAng;;


  /*ID*/
  if(dimage->doGrid)sprintf(waveID,"%s.%d.%d",dimage->waveID,(int)dimage->coord[0],(int)dimage->coord[1]);
  else if(dimage->useID)strcpy(waveID,dimage->waveID);
  else                  sprintf(waveID,"%d",numb);
  idLength=(hdfData->idLength<((int)strlen(waveID)+1))?hdfData->idLength:(int)strlen(waveID)+1;
  memcpy(&hdfData->waveID[numb*hdfData->idLength],waveID,idLength);

  /*waveform*/
  nBins=(hdfData->nBins<(waves->nBins-start))?hdfData->nBins:waves->nBins-start;
  for(j=0;j<hdfData->nTypeWaves;j++){
    memcpy(&hdfData->wave[j][numb*hdfData->nBins],&waves->wave[j][start],nBins*sizeof(float));
  }

  /*ground variables if using*/
  if(dimage->ground){
    hdfData->gElev[numb]=waves->gElevSimp;
    hdfData->slope[numb]=waves->gSlopeSimp;
    for(j=0;j<hdfData->nTypeWaves;j++){
      memcpy(&hdfData->ground[j][numb*hdfData->nBins],&waves->ground[j][start],nBins*sizeof(float));
    }
  }

  /*increment counter*/
  dimage->hdfCount++;

  return;
}/*packGEDIhdf*/


/*##############################################*/
/*check whether we are to overwrite or not*/

void checkWaveOverwrite(control *dimage,int numb)
{
  FILE *opoo=NULL;

  /*set output filename*/
  if(dimage->doGrid){
    sprintf(dimage->waveNamen,"%s.%d.%d.wave",dimage->outNamen,(int)dimage->coord[0],(int)dimage->coord[1]);
  }else if(dimage->readALSonce){
    sprintf(dimage->waveNamen,"%s.%s.wave",dimage->outNamen,dimage->waveIDlist[numb]);
  }else{
    strcpy(dimage->waveNamen,dimage->outNamen);
  }

  /*see if the file exists if we are checking for overwriting*/
  if(dimage->overWrite)dimage->useFootprint=1;
  else{
    if((opoo=fopen(dimage->waveNamen,"r"))==NULL){
      dimage->useFootprint=1;
    }else{
      dimage->useFootprint=0;
    }
    if(opoo){
      fclose(opoo);
      opoo=NULL;
    }
  }

  return;
}/*checkWaveOverwrite*/


/*##############################################*/
/*update footprint cordinate*/

void updateCoord(control *dimage,int i,int j)
{

  if(dimage->doGrid){
    dimage->coord[0]=dimage->gMinX+(double)i*dimage->gRes;
    dimage->coord[1]=dimage->gMinY+(double)j*dimage->gRes;
  }else if(dimage->readALSonce){
    dimage->coord[0]=dimage->coords[i][0];
    dimage->coord[1]=dimage->coords[i][1];
  }

  return;
}/*updateCoord*/


/*##############################################*/
/*determine ground properties with e DEM*/

void groundFromDEM(pCloudStruct **data,control *dimage,waveStruct *waves)
{
  int nX=0,nY=0,nBins=0;
  float res=0,rRes=0;
  float *gWave=NULL;
  float *waveFromDEM(double *,int,int,float,double,double,double,double,float,float,double *,int *);
  double minX=0,minY=0,maxX=0,maxY=0;
  double *gDEM=NULL,minZ=0;
  double *findGroundPoly(pCloudStruct **,int,double *,double *,double *,double *,float,int *,int *,double);
  double *findGroundNN(pCloudStruct **,int,double *,double *,float,int *,int *,double);
  void groundProperties(float *,int,double,float,waveStruct *,float);

  res=0.1;
  rRes=0.15;

  /*make DEM*/
  if(dimage->doGrid){
    minX=dimage->minX;
    maxX=dimage->maxX;
    minY=dimage->minY;
    maxY=dimage->maxY;
  }else{
    minX=maxX=100000000000.0;
    maxX=maxY=-100000000000.0;
  }

  if(dimage->polyGr)   gDEM=findGroundPoly(data,dimage->nFiles,&minX,&minY,&maxX,&maxY,res,&nX,&nY,waves->groundBreakElev);
  else if(dimage->nnGr)gDEM=findGroundNN(data,dimage->nFiles,&minX,&minY,res,&nX,&nY,waves->groundBreakElev);

  if(gDEM){
    /*make gap filled ground waveform*/
    gWave=waveFromDEM(gDEM,nX,nY,res,minX,minY,dimage->coord[0],dimage->coord[1],dimage->fSigma,rRes,&minZ,&nBins);
    TIDY(gDEM);

    /*ground properties*/
    groundProperties(gWave,nBins,minZ,rRes,waves,dimage->fSigma);
  }else{
    waves->gElev=waves->gSlope=-10000.0;
  }

  TIDY(gWave);
  return;
}/*groundFromDEM*/


/*####################################*/
/*ground properties*/

void groundProperties(float *gWave,int nBins,double minZ,float rRes,waveStruct *waves,float fSigma)
{
  int i=0;
  float total=0;
  double z=0,gStdev=0;
  char hasGround=0;

  /*CofG*/
  hasGround=1;
  waves->gElev=0.0;
  total=0.0;
  for(i=0;i<nBins;i++){
    z=(double)i*(double)rRes+minZ;
    waves->gElev+=z*(double)gWave[i];
    total+=gWave[i];
  }
  if(total>0.0){
    waves->gElev/=total;
    for(i=0;i<nBins;i++)gWave[i]/=total;
  }else{
    hasGround=0;
  }

  if(hasGround){
    /*stdev*/
    gStdev=total=0.0;
    for(i=0;i<nBins;i++){
      z=(double)i*(double)rRes+minZ;
      gStdev+=pow(z-waves->gElev,2.0)*(double)gWave[i];
      total+=gWave[i];
    }
    gStdev=sqrt(gStdev/total);
    waves->gSlope=atan2(sqrt(gStdev*gStdev),fSigma)*180.0/M_PI;
  }else waves->gSlope=waves->gElev=-10000.0;

  return;
}/*groundProperties*/


/*####################################*/
/*make waveform from DEM*/

float *waveFromDEM(double *gDEM,int nX,int nY,float res,double minX,double minY,double x0,double y0,float fSigma,float rRes,double *minZ,int *nBins)
{
  int i=0,j=0;
  int bin=0,place=0;
  float *gWave=NULL;
  float weight=0;
  double maxZ=0,x=0,y=0;
  double sep=0;

  /*mind min and max*/
  maxZ=-10000000.0;
  *minZ=1000000000.0;


  for(i=nX*nY-1;i>=0;i--){
    if(gDEM[i]<*minZ)*minZ=gDEM[i];
    if(gDEM[i]>maxZ)maxZ=gDEM[i];
  }

  *nBins=(int)((maxZ-*minZ)/rRes+0.5);
  gWave=falloc(*nBins,"ground wave",3);
  for(i=*nBins-1;i>=0;i--)gWave[i]=0.0;

  /*add up DEM*/
  for(i=0;i<nX;i++){
    x=((double)i+0.5)*(double)res+minX;
    for(j=0;j<nY;j++){
      y=((double)j+0.5)*(double)res+minY;
      place=j*nX+i;

      sep=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
      bin=(int)((gDEM[place]-(*minZ))/rRes);
      weight=(float)gaussian(sep,(double)fSigma,0.0);
      if((bin>=0)&&(bin<(*nBins)))gWave[bin]+=weight;
    }
  }

  return(gWave);
}/*waveFromDEM*/


/*####################################*/
/*write name if overlap*/

void checkThisFile(lasFile *las,control *dimage,int i)
{
  if(checkFileBounds(las,dimage->minX,dimage->maxX,dimage->minY,dimage->maxY)){
    fprintf(stdout,"Need %s\n",dimage->inList[i]);
  }
  return;
}/*checkThisFile*/


/*####################################*/
/*read lasFile and save relevant data*/

pCloudStruct *readALSdata(lasFile *las,control *dimage)
{
  int j=0;
  uint32_t i=0;
  uint32_t pUsed=0;    /*number of points used*/
  double x=0,y=0,z=0;
  pCloudStruct *data=NULL;
  char hasWave=0;   /*has waveform data, to save RAM*/
  char useFile=0,usePoint=0;
  char checkMultiFiles(lasFile *,int,double **,double);
  char checkMultiPoints(double,double,double,int,double **,double);

  /*allocate maximum number of points*/
  if(!(data=(pCloudStruct *)calloc(1,sizeof(pCloudStruct)))){
    fprintf(stderr,"error pCloudStruct allocation.\n");
    exit(1);
  }

  /*set nonsense bounds*/
  data->bounds[0]=data->bounds[1]=data->bounds[2]=10000000000.0;
  data->bounds[3]=data->bounds[4]=data->bounds[5]=-10000000000.0;

  /*is file needed?*/
  if(!dimage->readALSonce)useFile=checkFileBounds(las,dimage->globMinX,dimage->globMaxX,dimage->globMinY,dimage->globMaxY);
  else                    useFile=checkMultiFiles(las,dimage->gNx,dimage->coords,dimage->maxSep);

  /*check file is needed*/
  if(useFile){
    data->x=dalloc(las->nPoints,"x",0);
    data->y=dalloc(las->nPoints,"y",0);
    data->z=dalloc(las->nPoints,"z",0);
    data->refl=ialloc(las->nPoints,"refl",0);
    data->class=uchalloc((uint64_t)las->nPoints,"class",0);
    data->nRet=challoc((uint64_t)las->nPoints,"nRet",0);
    data->retNumb=challoc((uint64_t)las->nPoints,"nRet",0);
    data->scanAng=challoc((uint64_t)las->nPoints,"scanAng",0);
    data->packetDes=uchalloc((uint64_t)las->nPoints,"packetDes",0);
    data->grad=fFalloc(las->nPoints,"grad",0);
    for(i=0;i<las->nPoints;i++)data->grad[i]=falloc(3,"grad",i+1);
    data->time=falloc(las->nPoints,"time",0);              /*time in picoseconds of this wave*/
    if(!(data->waveMap=(uint64_t *)calloc(las->nPoints,sizeof(uint64_t)))){
      fprintf(stderr,"error in input filename structure.\n");
      exit(1);
    }
    if(!(data->waveLen=(uint32_t *)calloc(las->nPoints,sizeof(uint32_t)))){
      fprintf(stderr,"error in input filename structure.\n");
      exit(1);
    }

    /*loop over points*/
    pUsed=0;
    hasWave=0;
    for(i=0;i<las->nPoints;i++){
      /*read one point*/
      readLasPoint(las,i);
      setCoords(&x,&y,&z,las);
 
      /*is the point is of use?*/
      if(!dimage->readALSonce){
        if((x>=dimage->globMinX)&&(x<=dimage->globMaxX)&&(y>=dimage->globMinY)&&(y<=dimage->globMaxY)&&\
           (z>-10000.0)&&(z<10000.0)&&(fabs((float)las->scanAng)<=dimage->maxScanAng))usePoint=1;
        else usePoint=0;
      }else usePoint=checkMultiPoints(x,y,z,dimage->gNx,dimage->coords,dimage->maxSep);

      /*if we need to use point*/
      if(usePoint){
        data->x[pUsed]=x;
        data->y[pUsed]=y;
        data->z[pUsed]=z;
        if(las->refl>0)data->refl[pUsed]=(int)las->refl;
        else           data->refl[pUsed]=1;
        data->class[pUsed]=las->classif;
        data->nRet[pUsed]=(char)las->field.nRet;
        data->retNumb[pUsed]=(char)las->field.retNumb;
        data->scanAng[pUsed]=las->scanAng;

        /*determine data bounds*/
        if(x<data->bounds[0])data->bounds[0]=x;
        if(y<data->bounds[1])data->bounds[1]=y;
        if(z<data->bounds[2])data->bounds[2]=z;
        if(x>data->bounds[3])data->bounds[3]=x;
        if(y>data->bounds[4])data->bounds[4]=y;
        if(z>data->bounds[5])data->bounds[5]=z;

        /*record waveform if needed*/
        if(checkOneWave(las)){
          hasWave=1;
          data->packetDes[pUsed]=las->packetDes;
          for(j=0;j<3;j++)data->grad[pUsed][j]=las->grad[j];
          data->time[pUsed]=las->time;
          data->waveMap[pUsed]=las->waveMap;
          data->waveLen[pUsed]=las->waveLen;
        }else{
          data->packetDes[pUsed]=0;
          data->grad[pUsed][0]=data->grad[pUsed][1]=data->grad[pUsed][2]=0.0;
        }

        /*count points here*/
        pUsed++;
      }
    }/*point loop*/

    /*trim data arrays*/
    data->nPoints=pUsed;
    if(pUsed>0){
      if(!(data->x=(double *)realloc(data->x,data->nPoints*sizeof(double)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->y=(double *)realloc(data->y,data->nPoints*sizeof(double)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->z=(double *)realloc(data->z,data->nPoints*sizeof(double)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->refl=(int *)realloc(data->refl,data->nPoints*sizeof(int)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->class=(unsigned char *)realloc(data->class,data->nPoints*sizeof(unsigned char)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->nRet=(char *)realloc(data->nRet,data->nPoints*sizeof(char)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->retNumb=(char *)realloc(data->retNumb,data->nPoints*sizeof(char)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->scanAng=(char *)realloc(data->scanAng,data->nPoints*sizeof(char)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(dimage->useShadow){
        for(i=data->nPoints;i<las->nPoints-data->nPoints;i++)TIDY(data->grad[i]);
        if(!(data->grad=(float **)realloc(data->grad,data->nPoints*sizeof(float *)))){
          fprintf(stderr,"Balls\n");
          exit(1);
        }
      }else if(hasWave==0){
        TTIDY((void **)data->grad,las->nPoints);
        data->grad=NULL;
      }
    }else{
      TIDY(data->x);
      TIDY(data->y);
      TIDY(data->z);
      TIDY(data->refl);
      TIDY(data->class);
      TTIDY((void **)data->grad,las->nPoints);
      data->grad=NULL;
    }
    if(hasWave==1){
      data->waveStart=las->waveStart;
      if(!(data->packetDes=(unsigned char *)realloc(data->packetDes,data->nPoints*sizeof(unsigned char)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->time=(float *)realloc(data->time,data->nPoints*sizeof(float)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->waveMap=(uint64_t *)realloc(data->waveMap,data->nPoints*sizeof(uint64_t)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      if(!(data->waveLen=(uint32_t *)realloc(data->waveLen,data->nPoints*sizeof(uint32_t)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
      for(i=data->nPoints;i<las->nPoints-data->nPoints;i++)TIDY(data->grad[i]);
      if(!(data->grad=(float **)realloc(data->grad,data->nPoints*sizeof(float *)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
    }else{  /*clear out all the waveform bits*/
      TIDY(data->packetDes);
      if(!dimage->useShadow){
        TTIDY((void **)data->grad,las->nPoints);
        data->grad=NULL;
      }
      TIDY(data->time);
      TIDY(data->waveMap);
      TIDY(data->waveLen);
    }
  }else{/*file bounds check*/
    data->nPoints=0;
    data->nRet=NULL;
    data->retNumb=NULL;
    data->packetDes=NULL;
    data->grad=NULL;
    data->time=NULL;
    data->waveMap=NULL;
    data->waveLen=NULL;
    data->x=NULL;
    data->y=NULL;
    data->z=NULL;
    data->refl=NULL;
  }

  data->hasWave=hasWave;
  if(dimage->readWave&&hasWave){  /*only leave files open for waveform*/
    data->ipoo=las->ipoo;
    las->ipoo=NULL;
  }else{
    data->ipoo=NULL;
  }

  return(data);
}/*readALSdata*/


/*##########################################################################*/
/*see if we need to use this file when batch processing ALS data*/

char checkMultiFiles(lasFile *las,int nCoords,double **coords,double maxSep)
{
  int i=0;
  double maxX=0,minX=0;
  double maxY=0,minY=0;
  char useFile=0;

  useFile=0;
  for(i=0;i<nCoords;i++){
    minX=coords[i][0]-maxSep;
    maxX=coords[i][0]+maxSep;
    minY=coords[i][1]-maxSep;
    maxY=coords[i][1]+maxSep;
    if((las->minB[0]<=maxX)&&(las->minB[1]<=maxY)&&(las->maxB[0]>=minX)&&(las->maxB[1]>=minY)){
      useFile=1;
      break;
    }
  }/*footprint loop*/

  return(useFile);
}/*checkMultiFiles*/


/*##########################################################################*/
/*see if we need to use this point when batch processing ALS data*/

char checkMultiPoints(double x,double y,double z,int nCoords,double **coords,double maxSep)
{
  int i=0;
  double maxX=0,minX=0;
  double maxY=0,minY=0;
  char usePoint=0;

  usePoint=0;
  for(i=0;i<nCoords;i++){
    minX=coords[i][0]-maxSep;
    maxX=coords[i][0]+maxSep;
    minY=coords[i][1]-maxSep;
    maxY=coords[i][1]+maxSep;

    if((x>=minX)&&(x<=maxX)&&(y>=minY)&&(y<=maxY)){
      usePoint=1;
      break;
    }
  }/*footprint loop*/

  return(usePoint);
}/*checkMultiPoints*/


/*####################################*/
/*write GEDI waveforms*/

void writeGEDIwave(control *dimage,waveStruct *waves,int numb)
{
  int i=0,j=0;
  float r=0;
  char waveID[200];
  FILE *opoo=NULL;


  /*make waveID if needed*/
  if(dimage->doGrid){
    if(dimage->useID)sprintf(waveID,"%s.%d.%d",dimage->waveID,(int)dimage->coord[0],(int)dimage->coord[1]);
  }else if(dimage->readALSonce){
    strcpy(waveID,dimage->waveIDlist[numb]);
  }else{
    if(dimage->useID)strcpy(waveID,dimage->waveID);
  }


  if((opoo=fopen(dimage->waveNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->waveNamen);
    exit(1);
  }

  /*write header*/
  if(dimage->ground==0)fprintf(opoo,"# 1 elevation, 2 discrete intensity, 3 discrete count, 4 discrete fraction, 5 ALS pulse, 6 ALS and GEDI pulse, 7 ind decon, 8 ind decon GEDI, 9 decon GEDI, 10 ind decon\n");
  else                 fprintf(opoo,"# 1 elevation, 2 discrete intensity, 3 int canopy, 4 int ground, 5 discrete count, 6 count canopy, 7 count ground, 8 discrete fraction, 9 fraction canopy, 10 fraction ground, 11 ALS pulse, 12 ALS and GEDI pulse, 13 ind decon, 14 ind decon GEDI, 15 decon GEDI, 16 ind decon\n");
  fprintf(opoo,"# fSigma %f pSigma %f res %f sideLobes %d\n",dimage->fSigma,dimage->pSigma,dimage->res,dimage->sideLobe);
  fprintf(opoo,"# coord %.2f %.2f\n",dimage->coord[0],dimage->coord[1]);
  fprintf(opoo,"# density point %f beam %f\n",dimage->pointDense,dimage->beamDense);
  fprintf(opoo,"# meanScanAng %f\n",waves->meanScanAng);
  if(dimage->useID)fprintf(opoo,"# waveID %s\n",waveID);
  if(dimage->ground&&(dimage->polyGr||dimage->nnGr))fprintf(opoo,"# ground %f %f\n",waves->gElev,waves->gSlope);
  if(dimage->ground)fprintf(opoo,"# simpleGround %f\n",waves->gElevSimp);

  /*write data*/
  for(i=0;i<waves->nBins;i++){
    r=(float)waves->maxZ-(float)i*dimage->res;

    fprintf(opoo,"%f",r);
    for(j=0;j<waves->nWaves;j++){
      fprintf(opoo," %f",waves->wave[j][i]);
       if(dimage->ground&&(j<3))fprintf(opoo," %f %f",waves->canopy[j][i],waves->ground[j][i]);
    }
    fprintf(opoo,"\n");
  }

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",dimage->waveNamen);
  return;
}/*writeGEDIwave*/


/*####################################*/
/*make GEDI waveforms*/

waveStruct *makeGediWaves(control *dimage,pCloudStruct **data)
{
  int j=0,k=0;
  float tot=0;
  waveStruct *waves=NULL;
  waveStruct *allocateGEDIwaves(control *,pCloudStruct **);
  void gediFromWaveform(pCloudStruct *,uint32_t,float,waveStruct *,control *);
  void processAggragate(control *,waveStruct *);
  void checkFootCovered(control *);
  void cleanOutliers(waveStruct *,control *);
  void waveFromPointCloud(control *,pCloudStruct **,waveStruct *);
  void waveFromShadows(control *,pCloudStruct **,waveStruct *);
  void determineALScoverage(control *,pCloudStruct **);
  denPar *setDeconForGEDI(control *);


  /*determine ALS coverage*/
  determineALScoverage(dimage,data);

  /*check that whole footprint is covered*/
  if(dimage->checkCover)checkFootCovered(dimage);
  else                  dimage->useFootprint=1;

  /*only if it contains data*/
  if(dimage->useFootprint){
    /*allocate*/
    waves=allocateGEDIwaves(dimage,data);

    /*set up denoising if using*/
    if(dimage->doDecon)dimage->decon=setDeconForGEDI(dimage);

    /*make waves*/
    if(dimage->useShadow==0)waveFromPointCloud(dimage,data,waves);
    else                    waveFromShadows(dimage,data,waves);

    /*clean outliers if needed*/
    if(dimage->cleanOut)cleanOutliers(waves,dimage);
    else                waves->groundBreakElev=-100000000.0;

    /*deconvolve aggragated waveform*/
    if(dimage->doDecon)processAggragate(dimage,waves);

    /*normalise integral*/
    for(k=0;k<waves->nWaves;k++){
      tot=0.0;
      for(j=0;j<waves->nBins;j++)tot+=waves->wave[k][j]*dimage->res;
      if(tot>0.0){
        for(j=0;j<waves->nBins;j++)waves->wave[k][j]/=tot;
        if(dimage->ground&&(k<3)){
          for(j=0;j<waves->nBins;j++){
            waves->canopy[k][j]/=tot;
            waves->ground[k][j]/=tot;
          }
        }
      }
    }

    /*tidy arrays*/
    if(dimage->decon){
      TTIDY((void **)dimage->decon->pulse,2);
      dimage->decon->pulse=NULL;
      TIDY(dimage->decon);
    }
  }/*contains data*/

  /*check whether empty*/
  tot=0.0;
  for(j=0;j<waves->nBins;j++)tot+=waves->wave[0][j]*dimage->res;
  if((tot<TOL)||(waves->nBins==0))dimage->useFootprint=0;

  return(waves);
}/*makeGediWaves*/


/*################################################################################*/
/*determine ALS coverage*/

void determineALScoverage(control *dimage,pCloudStruct **data)
{
  int i=0;
  int gX=0,gY=0;
  uint32_t j=0;
  double dx=0,dy=0;
  double sepSq=0;
  float area=0.0;


  for(i=0;i<dimage->nFiles;i++){  /*file loop*/
    for(j=0;j<data[i]->nPoints;j++){ /*point loop*/
      /*check within bounds*/
      if((data[i]->x[j]>=dimage->minX)&&(data[i]->x[j]<=dimage->maxX)&&(data[i]->y[j]>=dimage->minY)&&(data[i]->y[j]<=dimage->maxY)){
        dx=data[i]->x[j]-dimage->coord[0];
        dy=data[i]->y[j]-dimage->coord[1];
        sepSq=dx*dx+dy*dy;

        /*ground grid for ALS coverage*/
        if(dimage->normCover||dimage->checkCover){
          if(data[i]->retNumb[j]==data[i]->nRet[j]){  /*only once per beam*/
            /*mark sampling desnity for normalisation*/
            gX=(int)((data[i]->x[j]-dimage->g0[0])/(double)dimage->gridRes);
            gY=(int)((data[i]->y[j]-dimage->g0[1])/(double)dimage->gridRes);
            if((gX>=0)&&(gX<dimage->gX)&&(gY>=0)&&(gY<dimage->gY)){
              dimage->nGrid[gY*dimage->gX+gX]++;
            }
          }
        }/*ground grid for ALS coverage*/

        /*point and beam density*/
        if(sepSq<=dimage->denseRadSq){
          dimage->pointDense+=1.0;
          if(data[i]->retNumb[j]==data[i]->nRet[j])dimage->beamDense+=1.0;
        }
      }/*bounds check*/
    }/*point loop*/
  }/*file loop*/

  area=M_PI*dimage->denseRadSq;
  dimage->pointDense/=area;
  dimage->beamDense/=area;

  return;
}/*determineALScoverage*/


/*################################################################################*/
/*make waveforms accounting for shadowing*/

void waveFromShadows(control *dimage,pCloudStruct **data,waveStruct *waves)
{
  int i=0;
  float **tempWave=NULL;
  //float iRes=0,grad[3];
  float grad[3];
  void voxelGap(control *,pCloudStruct **,waveStruct *);
  rImageStruct *rImage=NULL;    /*range image, a stack nBins long*/
  lidVoxPar lidPar;


  //iRes=0.02;
  grad[0]=grad[1]=0.0;
  grad[2]=-1.0;

  /*set lidar parameters for a downwards looking ALS*/
  lidPar.minRefl=1.0;             /*minimum refletance value to scale between 0 and 1*/
  lidPar.maxRefl=1.0;             /*maximum refletance value to scale between 0 and 1*/
  lidPar.appRefl=1.0 ;            /*scale between TLS reflectance and size*/
  lidPar.beamTanDiv=0.0;          /*tan of beam divergence*/
  lidPar.beamRad=dimage->beamRad; /*start radius*/
  lidPar.minGap=0.00001;          /*minimum gap fraction correction to apply*/

  /*gap fraction from voxelising data*/
  voxelGap(dimage,data,waves);

  /*create images*/
  //rImage=allocateRangeImage(dimage->nFiles,data,dimage->pRes*4.0,iRes,&(grad[0]),dimage->coord[0],dimage->coord[1],waves->maxZ);
  //rImage=allocateRangeImage(dimage->nFiles,data,NULL,0.15,0.01,&(grad[0]),dimage->coord[0],dimage->coord[1],waves->maxZ,NULL);
  fprintf(stderr,"THis method is no longer operational. Do not use\n");
  exit(1);


  silhouetteImage(dimage->nFiles,data,NULL,rImage,&lidPar,NULL,0,NULL);

  /*convert images to waveform*/
  tempWave=fFalloc(2,"",0);
  for(i=0;i<2;i++)tempWave[i]=falloc(rImage->nBins,"",i+1);
  waveFromImage(rImage,tempWave,1,dimage->fSigma);
  for(i=0;i<rImage->nBins;i++)fprintf(stdout,"%f %f %f\n",waves->maxZ-(double)i*rImage->rRes,tempWave[0][i],tempWave[1][i]);
  TTIDY((void **)tempWave,2);
  tempWave=NULL;


  if(rImage){
    TTIDY((void **)rImage->image,rImage->nBins);
    rImage->image=NULL;
    TIDY(rImage);
  }
  return;
}/*waveFromShadows*/


/*################################################################################*/
/*make a map of voxel gaps*/

void voxelGap(control *dimage,pCloudStruct **data,waveStruct *waves)
{
  int i=0,vInd=0;
  int xBin=0,yBin=0,zBin=0;
  uint32_t j=0;
  double bounds[6];
  voxStruct *vox=NULL;
  voxStruct *tidyVox(voxStruct *);
  voxStruct *voxAllocate(int,float *,double *,char);
  void countVoxGap(double,double,double,float *,voxStruct *,int,int,float,int);


  bounds[0]=dimage->minX;
  bounds[1]=dimage->minY;
  bounds[2]=waves->minZ;
  bounds[3]=dimage->maxX;
  bounds[4]=dimage->maxY;
  bounds[5]=waves->maxZ;


  /*first make a voxel map*/
  vox=voxAllocate(1,&(dimage->vRes[0]),&(bounds[0]),0);

  for(i=0;i<dimage->nFiles;i++){ /*file loop*/
    for(j=0;j<data[i]->nPoints;j++){  /*point loop*/
      countVoxGap(data[i]->x[j],data[i]->y[j],data[i]->z[j],&(data[i]->grad[j][0]),vox,1,1,dimage->beamRad,0);
    }/*point loop*/
  }/*file loop*/

  /*calculate gap fraction for each return*/
  for(i=0;i<dimage->nFiles;i++){ /*file loop*/
    data[i]->gap=falloc(data[i]->nPoints,"point gaps",i+1);
    for(j=0;j<data[i]->nPoints;j++){  /*point loop*/
      xBin=(int)((data[i]->x[j]-vox->bounds[0])/vox->res[0]+0.5);
      yBin=(int)((data[i]->y[j]-vox->bounds[1])/vox->res[1]+0.5);
      zBin=(int)((data[i]->z[j]-vox->bounds[2])/vox->res[2]+0.5);
      vInd=xBin+vox->nX*yBin+vox->nX*vox->nY*zBin;

      if((vox->hits[0][vInd]+vox->miss[0][vInd])>0.0)data[i]->gap[j]=vox->hits[0][vInd]/(vox->hits[0][vInd]+vox->miss[0][vInd]);
      else                                           data[i]->gap[j]=1.0;
    }/*point loop*/
  }/*file loop*/

  vox=tidyVox(vox);
  return;
}/*voxelGap*/


/*################################################################################*/
/*make waveform from point cloud*/

void waveFromPointCloud(control *dimage,pCloudStruct **data,waveStruct *waves)
{
  int numb=0,bin=0,j=0;
  int gX=0,gY=0,n=0;
  uint32_t i=0;
  double sep=0;
  double dX=0,dY=0;
  double totGround=0;     /*contrbution to ground estimate*/
  float refl=0,rScale=0,fracHit=0,totAng=0;
  void gediFromWaveform(pCloudStruct *,uint32_t,float,waveStruct *,control *);
  void processAggragate(control *,waveStruct *);
  void waveFromPointCloud(control *,pCloudStruct **,waveStruct *);
  denPar *setDeconForGEDI(control *);

  /*reset mean scan angle*/
  waves->meanScanAng=totAng=0.0;
  /*ground elevation estimate*/
  if(dimage->ground){
    waves->gElevSimp=0.0;
    totGround=0.0;
  }

  /*make waves*/
  for(n=0;n<dimage->nLobes;n++){
    for(numb=0;numb<dimage->nFiles;numb++){
      for(i=0;i<data[numb]->nPoints;i++){
        dX=data[numb]->x[i]-dimage->lobe[n].coord[0];
        dY=data[numb]->y[i]-dimage->lobe[n].coord[1];
        sep=sqrt(dX*dX+dY*dY);

        if(dimage->topHat==0)rScale=(float)gaussian(sep,(double)dimage->lobe[n].fSigma,0.0);
        else{
          if(sep<=dimage->lobe[n].maxSepSq)rScale=1.0;
          else                             rScale=0.0;
        }

        if(rScale>dimage->iThresh){  /*if bright enough to matter*/
          /*scale by sampling density*/
          if(dimage->normCover){
            gX=(int)((data[numb]->x[i]-dimage->g0[0])/(double)dimage->gridRes);
            gY=(int)((data[numb]->y[i]-dimage->g0[1])/(double)dimage->gridRes);
            if((gX>=0)&&(gX<dimage->gX)&&(gY>=0)&&(gY<dimage->gY)){
              if(dimage->nGrid[gY*dimage->gX+gX]>0)rScale/=(float)dimage->nGrid[gY*dimage->gX+gX];
            }
          }/*scale by sampling density*/


          /*discrete return*/
          refl=(float)data[numb]->refl[i]*rScale;
          if(data[numb]->nRet[i]>0)fracHit=1.0/(float)data[numb]->nRet[i];
          else                     fracHit=1.0;
          for(j=0;j<dimage->pulse->nBins;j++){
            bin=(int)((waves->maxZ-data[numb]->z[i]+(double)dimage->pulse->x[j])/(double)dimage->res);
            if((bin>=0)&&(bin<waves->nBins)){
              waves->wave[0][bin]+=refl*dimage->pulse->y[j];
              waves->wave[1][bin]+=rScale*dimage->pulse->y[j];
              waves->wave[2][bin]+=rScale*fracHit*dimage->pulse->y[j];
              if(dimage->ground){
                if(data[numb]->class[i]==2){
                  waves->ground[0][bin]+=refl*dimage->pulse->y[j];
                  waves->ground[1][bin]+=rScale*dimage->pulse->y[j];
                  waves->ground[2][bin]+=rScale*fracHit*dimage->pulse->y[j];
                }else{
                  waves->canopy[0][bin]+=refl*dimage->pulse->y[j];
                  waves->canopy[1][bin]+=rScale*dimage->pulse->y[j];
                  waves->canopy[2][bin]+=rScale*fracHit*dimage->pulse->y[j];
                }
              }/*ground recording if needed*/
            }/*bin bound check*/
          }/*pulse bin loop*/
          if(dimage->ground){
            if(data[numb]->class[i]==2){
              waves->gElevSimp+=rScale*data[numb]->z[i];
              totGround+=rScale;
            }
          }
          waves->meanScanAng+=rScale*fracHit*(float)abs((int)data[numb]->scanAng[i]);
          totAng+=rScale*fracHit;

          /*full-waveform*/
          if(dimage->readWave&&data[numb]->hasWave){
            if(data[numb]->packetDes[i]){  /*test for waveform*/
              gediFromWaveform(data[numb],i,rScale,waves,dimage);
            }
          }/*waveform test*/
        }
      }/*point loop*/
    }/*file loop*/
  }/*lobe loop*/

  /*normalise mean scan angle*/
  if(totAng>0.0)waves->meanScanAng/=totAng;
  if(totGround>=0.0)waves->gElevSimp/=totGround;
  else              waves->gElevSimp=-9999.0;

  return;
}/*waveFromPointCloud*/


/*####################################*/
/*clean outlier points from waveform*/

void cleanOutliers(waveStruct *waves,control *dimage)
{
  int i=0,j=0,gStart=0;
  char pastGround=0;
  float gGap=0;  /*gap in ground return*/
  float maxGap=0;
  float maxGround=0,gThresh=0;
  float max=0,thresh=0;

  if(!dimage->ground){
    fprintf(stderr,"No need to clean without ground\n");
    exit(1);
  }

  maxGap=10.0;  /*maximum permittable gap in the ground return*/
  gGap=0.0;
  pastGround=0;

  /*determine max ground and max return*/
  maxGround=max=0.0;
  for(i=0;i<waves->nBins;i++){
    if(waves->ground[1][i]>maxGround)maxGround=waves->ground[1][i];
    if(waves->wave[1][i]>max)max=waves->wave[1][i];
  }
  gThresh=maxGround*0.01;
  thresh=max*0.001;

  for(i=0;i<waves->nBins;i++){
    if(waves->ground[1][i]>=gThresh){
      if(pastGround==0)gStart=i;
      pastGround=1;
      gGap=0.0;
    }else{
      if(pastGround)gGap+=dimage->res;
    }
    if(gGap>maxGap){  /*too big a break, delete*/
      waves->groundBreakElev=waves->maxZ-(double)i*(double)dimage->res;
      for(;i<waves->nBins;i++){
        waves->canopy[0][i]=waves->canopy[1][i]=waves->ground[0][i]=waves->ground[1][i]=0.0;
        for(j=0;j<waves->nWaves;j++)waves->wave[j][i]=0.0;
      }
    }/*too big a break, delete*/
  }

  /*look for above canopy outliers*/
  gGap=0.0;
  maxGap=50.0;
  for(i=gStart;i>=0;i--){
    if(waves->wave[1][i]>=thresh)gGap=0.0;
    gGap+=dimage->res;

    if(gGap>=maxGap){
      for(;i>=0;i--){
        waves->canopy[0][i]=waves->canopy[1][i]=waves->ground[0][i]=waves->ground[1][i]=0.0;
        for(j=0;j<waves->nWaves;j++)waves->wave[j][i]=0.0;
      }
    }
  }

  return;
}/*cleanOutliers*/


/*####################################*/
/*check footprint is covered by ALS*/

void checkFootCovered(control *dimage)
{
  int i=0,j=0,nWithin=0;
  int thresh=0,nMissed=0;
  double dX=0,dY=0,sepSq=0;
  double useRad=0,radSq=0;

  useRad=10.0;
  if(useRad>(double)dimage->fSigma)useRad=(double)dimage->fSigma;
  radSq=useRad*useRad;

  for(i=0;i<dimage->gX;i++){
    dX=(double)i*(double)dimage->gridRes-dimage->coord[0];
    for(j=0;j<dimage->gY;j++){
      dY=(double)j*(double)dimage->gridRes-dimage->coord[1];
      sepSq=dX*dX+dY*dY;

      if(sepSq<radSq){
        if(dimage->nGrid[j*dimage->gX+i]==0)nMissed++;
        nWithin++;
      }
    }/*y loop*/
  }/*x loop*/

  thresh=(int)((float)nWithin*2.0/3.0);
  if(nMissed>thresh){
    fprintf(stderr,"Too many missed %d of %d\n",nMissed,nWithin);
    dimage->useFootprint=0;
  }else dimage->useFootprint=1;

  return;
}/*checkFootCovered*/


/*####################################*/
/*deconvolve aggragated wave*/

void processAggragate(control *dimage,waveStruct *waves)
{
  int i=0;
  float *processed=NULL,*smooPro=NULL;
  float *smooth(float,int,float *,float);

  /*Add background noise back*/
  for(i=0;i<waves->nBins;i++)waves->wave[7][i]+=dimage->meanN;

  /*deconvolve and reconvolve*/
  processed=processFloWave(&(waves->wave[7][0]),(int)waves->nBins,dimage->decon,1.0);
  for(i=0;i<waves->nBins;i++)waves->wave[8][i]=processed[i];
  smooPro=smooth(dimage->pSigma,(int)waves->nBins,processed,dimage->res);
  TIDY(processed);

  TIDY(waves->wave[7]);
  waves->wave[7]=smooPro;
  smooPro=NULL;

  return;
}/*processAggragate*/


/*####################################*/
/*GEDI wave from ALS waveforms*/

void gediFromWaveform(pCloudStruct *data,uint32_t i,float rScale,waveStruct *waves,control *dimage)
{
  int j=0,bin=0;
  int buffBins=0;
  uint32_t waveLen=0;
  float grad[3],*smoothed=NULL,*floWave=NULL;
  float *processed=NULL,*smooPro=NULL;
  float *smooth(float,int,float *,float);
  double x=0,y=0,z=0;
  unsigned char *wave=NULL,*temp=NULL;

  for(j=0;j<3;j++)grad[j]=data->grad[j][i];
  wave=readLasWave(data->waveMap[i],data->waveLen[i],data->ipoo,data->waveStart);

  /*buffer to give space for smoothing*/
  buffBins=80;
  waveLen=data->waveLen[i]+(uint32_t)(2*buffBins);
  temp=uchalloc((uint64_t)waveLen,"temp waveform",0);
  for(j=0;j<buffBins;j++)temp[j]=(unsigned char)dimage->meanN;
  for(j=0;j<(int)data->waveLen[i];j++)temp[j+buffBins]=wave[j];
  for(j=(int)data->waveLen[i]+buffBins;j<(int)waveLen;j++)temp[j]=(unsigned char)dimage->meanN;
  TIDY(wave);
  wave=temp;
  temp=NULL;

  /*deconvolve and reconvolve*/
  if(dimage->indDecon){
    processed=processWave(wave,(int)waveLen,dimage->decon,1.0);
    smooPro=smooth(dimage->pSigma,(int)waveLen,processed,dimage->res);
  }

  /*convolve with GEDI pulse*/
  floWave=falloc((int)waveLen,"",0);
  for(j=0;j<(int)waveLen;j++)floWave[j]=(float)wave[j]-dimage->meanN;
  smoothed=smooth(dimage->pSigma,(int)waveLen,floWave,dimage->res);
  TIDY(floWave);

  /*add up*/
  for(j=0;j<(int)waveLen;j++){
    binPosition(&x,&y,&z,j-buffBins,data->x[i],data->y[i],data->z[i],data->time[i],&(grad[0]));
    bin=(int)((waves->maxZ-z)/(double)dimage->res);
    if((bin>=0)&&(bin<waves->nBins)){
      waves->wave[3][bin]+=((float)wave[j]-dimage->meanN)*rScale;  /*with ALS pulse*/
      waves->wave[4][bin]+=smoothed[j]*rScale;
      if(dimage->doDecon){
        if(dimage->indDecon){
          waves->wave[5][bin]+=smooPro[j]*rScale;
          waves->wave[6][bin]+=processed[j]*rScale;
        }
        waves->wave[7][bin]+=((float)wave[j]-dimage->meanN)*rScale;
      }
    }
  }/*wave bin loop*/

  TIDY(wave);
  TIDY(smooPro);
  TIDY(smoothed);
  TIDY(processed);
  return;
}/*gediFromWaveform*/


/*####################################*/
/*set deconvolution parameters for GEDI*/

denPar *setDeconForGEDI(control *dimage)
{
  void setDenoiseDefault(denPar *);
  void readPulse(denPar *);
  denPar *decon=NULL;

  /*set defaults*/
  if(!(decon=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error decon structure allocation.\n");
    exit(1);
  }
  setDenoiseDefault(decon);

  /*particular values for here*/
  decon->deChang=pow(10.0,-8.0);  /*change between decon iterations to stop*/
  decon->thresh=17.0;
  decon->meanN=13.0;
  strcpy(decon->pNamen,"/Users/stevenhancock/data/bess/leica_shape/leicaPulse.dat");  /*pulse filename*/
  decon->deconMeth=0;     /*Gold's method*/
  decon->pScale=1.2;
  decon->noiseTrack=0;

  /*read system pulse*/
  readPulse(decon);

  return(decon);
}/*setDeconForGEDI*/


/*####################################*/
/*allocate wave structure*/

waveStruct *allocateGEDIwaves(control *dimage,pCloudStruct **data)
{
  int j=0,numb=0,k=0;
  uint32_t i=0;
  double maxZ=0,minZ=0;
  double buff=0;
  waveStruct *waves=NULL;
  char hasPoints=0;

  if(!(waves=(waveStruct *)calloc(1,sizeof(waveStruct)))){
    fprintf(stderr,"error waveStruct allocation.\n");
    exit(1);
  }

  /*determine wave bounds*/
  buff=35.0;
  minZ=100000000000.0;
  maxZ=-100000000000.0;
  hasPoints=0;
  for(numb=0;numb<dimage->nFiles;numb++){
    if(data[numb]->nPoints>0)hasPoints=1;
    for(i=0;i<data[numb]->nPoints;i++){
      if(data[numb]->z[i]>maxZ)maxZ=data[numb]->z[i];
      if(data[numb]->z[i]<minZ)minZ=data[numb]->z[i];
    }
  }/*bound finding*/

  if(hasPoints==0){
    fprintf(stderr,"No points included\n");
    exit(1);
  }

  waves->minZ=minZ-buff;
  waves->maxZ=maxZ+buff;

  waves->nBins=(int)((waves->maxZ-waves->minZ)/(double)dimage->res);
  waves->nWaves=9;
  waves->wave=fFalloc(waves->nWaves,"result waveform",0);
  for(j=0;j<waves->nWaves;j++){
    waves->wave[j]=falloc(waves->nBins,"result waveform",j+1);
    for(k=0;k<waves->nBins;k++)waves->wave[j][k]=0.0;
  }
  if(dimage->ground){
    waves->canopy=fFalloc(3,"canopy",0);
    waves->ground=fFalloc(3,"ground",0);
    for(j=0;j<3;j++){
      waves->canopy[j]=falloc(waves->nBins,"canopy waveform",j+1);
      waves->ground[j]=falloc(waves->nBins,"ground waveform",j+1);
      for(k=0;k<waves->nBins;k++)waves->canopy[j][k]=waves->ground[j][k]=0.0;
    }
  }

  return(waves);
}/*allocateGEDIwaves*/


/*####################################*/
/*set GEDI footprint*/

void setGediFootprint(control *dimage)
{
  int i=0;
  float totE=0;
  float az=0;
  double tX=0,tY=0;


  /*footprint width*/
  az=dimage->lobeAng*M_PI/180.0;  /*convert anlge to radians*/

  /*number of lobes and allocate*/
  if(dimage->sideLobe==0)dimage->nLobes=1;
  else                   dimage->nLobes=7;
  if(!(dimage->lobe=(lobeStruct *)calloc(dimage->nLobes,sizeof(lobeStruct)))){
    fprintf(stderr,"error lobeStruct allocation.\n");
    exit(1);
  }

  /*central footprint*/
  i=0;
  dimage->lobe[i].coord[0]=dimage->coord[0];
  dimage->lobe[i].coord[1]=dimage->coord[1];
  dimage->lobe[i].fSigma=dimage->fSigma;
  dimage->lobe[i].maxSepSq=(double)(dimage->maxSep*dimage->maxSep);
  if(dimage->sideLobe==0)dimage->lobe[i].E=1.0;
  else{  /*include side lobes*/
    totE=1.0+0.0599+0.0731+0.0317+0.0319+0.0167+0.0163;
    i=0;
    dimage->lobe[i].E=1.0/totE;

    /*first southern lobe*/
    i=1;
    dimage->lobe[i].E=0.0731/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]-20.0*sin(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]-20.0*cos(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    dimage->lobe[i].maxSepSq=(double)(dimage->maxSep*dimage->maxSep);

    /*first nothern lobe*/
    i=2;
    dimage->lobe[i].E=0.0599/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]+20.0*sin(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]+20.0*cos(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    dimage->lobe[i].maxSepSq=(double)(dimage->maxSep*dimage->maxSep);

    /*western lobe*/
    i=3;
    dimage->lobe[i].E=0.0319/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]-20.0*cos(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]-20.0*sin(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    dimage->lobe[i].maxSepSq=(double)(dimage->maxSep*dimage->maxSep);

    /*eastern lobe*/
    i=4;
    dimage->lobe[i].E=0.0317/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]+20.0*cos(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]+20.0*sin(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    dimage->lobe[i].maxSepSq=(double)(dimage->maxSep*dimage->maxSep);

    /*second southern lobe*/
    i=5;
    dimage->lobe[i].E=0.0167/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]-30.0*sin(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]-30.0*cos(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    dimage->lobe[i].maxSepSq=(double)(dimage->maxSep*dimage->maxSep);

    /*second northern lobe*/
    i=6;
    dimage->lobe[i].E=0.0163/totE;
    dimage->lobe[i].coord[0]=dimage->coord[0]+30.0*cos(az);
    dimage->lobe[i].coord[1]=dimage->coord[1]+30.0*sin(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    dimage->lobe[i].maxSepSq=(double)(dimage->maxSep*dimage->maxSep);
  }/*side lobe test*/


  /*determine min and max bounds*/
  dimage->minX=dimage->minY=100000000000.0;
  dimage->maxX=dimage->maxY=-1000000000000.0;
  for(i=0;i<dimage->nLobes;i++){
    tX=dimage->lobe[i].coord[0]-sqrt(dimage->lobe[i].maxSepSq);
    if(tX<dimage->minX)dimage->minX=tX;
    tX=dimage->lobe[i].coord[0]+sqrt(dimage->lobe[i].maxSepSq);
    if(tX>dimage->maxX)dimage->maxX=tX;
    tY=dimage->lobe[i].coord[1]-sqrt(dimage->lobe[i].maxSepSq);
    if(tY<dimage->minY)dimage->minY=tY;
    tY=dimage->lobe[i].coord[1]+sqrt(dimage->lobe[i].maxSepSq);
    if(tY>dimage->maxY)dimage->maxY=tY;
  }

  /*grid for normalising sampling*/
  if(dimage->normCover||dimage->checkCover){
    dimage->gridRes=1.5;
    dimage->gX=(int)((float)(dimage->maxX-dimage->minX)/dimage->gridRes)+2;
    dimage->gY=(int)((float)(dimage->maxY-dimage->minY)/dimage->gridRes)+2;
    dimage->g0[0]=dimage->minX+(double)dimage->gridRes;
    dimage->g0[1]=dimage->minY+(double)dimage->gridRes;
    dimage->nGrid=ialloc(dimage->gX*dimage->gY,"nGrid",0);
    for(i=dimage->gX*dimage->gY-1;i>=0;i--)dimage->nGrid[i]=0;
  }

  /*radius to calculate density within*/
  if(dimage->topHat==0)dimage->denseRadSq=dimage->fSigma*dimage->fSigma*4.0;
  else                 dimage->denseRadSq=dimage->fSigma;
  dimage->pointDense=dimage->beamDense=0.0;

  return;
}/*setGediFootprint*/


/*####################################*/
/*set GEDI grid or batch*/

void setGediGrid(control *dimage)
{
  void readFeetList(control *);

  /*footprint width*/
  if(dimage->fSigma<0.0)dimage->fSigma=dimage->fWidth;
  if(dimage->topHat==0)dimage->maxSep=determineGaussSep(dimage->fSigma,dimage->iThresh*1.0);
  else                 dimage->maxSep=dimage->fSigma;

  if(dimage->doGrid){  /*it is a grid*/
    /*number of footprints*/
    dimage->gNx=(int)((dimage->gMaxX-dimage->gMinX)/dimage->gRes+1);
    dimage->gNy=(int)((dimage->gMaxY-dimage->gMinY)/dimage->gRes+1);

    /*global bounds*/
    dimage->globMinX=dimage->gMinX-dimage->maxSep;
    dimage->globMaxX=dimage->gMaxX+dimage->maxSep;
    dimage->globMinY=dimage->gMinY-dimage->maxSep;
    dimage->globMaxY=dimage->gMaxY+dimage->maxSep;
  }else if(dimage->readALSonce){ /*it is a batch*/
    /*read list of coords*/
    readFeetList(dimage);
  }else{   /*single footprint*/
    dimage->gNx=dimage->gNy=1;
    dimage->globMinX=dimage->coord[0]-dimage->maxSep;
    dimage->globMaxX=dimage->coord[0]+dimage->maxSep;
    dimage->globMinY=dimage->coord[1]-dimage->maxSep;
    dimage->globMaxY=dimage->coord[1]+dimage->maxSep;
  }

  if((dimage->gNx*dimage->gNy)>dimage->gediIO.nMessages)dimage->gediIO.nMessages=(int)(dimage->gNx*dimage->gNy/dimage->gediIO.nMessages);
  else                                                  dimage->gediIO.nMessages=1;

  return;
}/*setGediGrid*/


/*####################################*/
/*read list of coordinates*/

void readFeetList(control *dimage)
{
  int i=0;
  char line[200],temp1[50];
  char temp2[50],temp3[100];;
  FILE *ipoo=NULL;

  /*open file*/
  if((ipoo=fopen(dimage->coordList,"r"))==NULL){
    fprintf(stderr,"Error opening input file list \"%s\"\n",dimage->coordList);
    exit(1);
  }


  /*count number of lines*/
  i=0;
  while(fgets(line,200,ipoo)!=NULL)if(strncasecmp(line,"#",1))i++;

  /*allocate space*/
  dimage->gNx=i;
  dimage->gNy=1;
  dimage->coords=dDalloc(dimage->gNx,"coord list",0);
  dimage->waveIDlist=chChalloc(dimage->gNx,"wave ID list",0);

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){ 
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read coordinate list*/
  i=0;
  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      dimage->coords[i]=dalloc(2,"coord list",i+1);
      if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){ /*read coord and waveID*/
        dimage->coords[i][0]=atof(temp1);
        dimage->coords[i][1]=atof(temp2);
        dimage->waveIDlist[i]=challoc((int)strlen(temp3)+1,"wave ID list",i+1);
        strcpy(dimage->waveIDlist[i],temp3);
      }else if(sscanf(line,"%s %s",temp1,temp2)==2){
        dimage->coords[i][0]=atof(temp1);
        dimage->coords[i][1]=atof(temp2);
        sprintf(temp3,"%f.%f",dimage->coords[i][0],dimage->coords[i][1]);
        dimage->waveIDlist[i]=challoc((int)strlen(temp3)+1,"wave ID list",i+1);
        strcpy(dimage->waveIDlist[i],temp3);
      }else{
        fprintf(stderr,"coord list reading error \"%s\"\n",line);
        exit(1);
      }
      i++;
    }
  }

  /*close file*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return;
}/*readFeetList*/


/*####################################*/
/*set GEDI pulse*/

void setGediPulse(control *dimage)
{
  int i=0;
  float fwhm=0;   /*FWHM in metres*/
  float x=0,y=0;
  float max=0,tot=0;
  void readSimPulse(control *);


  if(!(dimage->pulse=(pulseStruct *)calloc(1,sizeof(pulseStruct)))){
    fprintf(stderr,"error pulseStruct allocation.\n");
    exit(1);
  }


  if(dimage->readPulse==0){  /*Gaussian pulse*/
    /*pulse length*/
    /*calculate sigma from FWHM*/
    if(dimage->pSigma<0.0){  /*GEDI unless specificed*/
      fwhm=dimage->pFWHM*0.2998/2.0;  /*time for two way*/
      dimage->pSigma=fwhm/2.35482;  /* =2*sqrt(2*ln2) */
    }

    if(dimage->pSigma>0.0){  /*if we are using a pulse width*/
      /*determine number of bins*/
      dimage->pulse->nBins=0;
      x=0.0;
      do{
        y=(float)gaussian((double)x,(double)dimage->pSigma,0.0);
        x+=dimage->pRes;
        dimage->pulse->nBins+=2;  /*both sides of peak*/
      }while(y>=dimage->iThresh);
  
      dimage->pulse->x=falloc(dimage->pulse->nBins,"pulse x",0);
      dimage->pulse->y=falloc(dimage->pulse->nBins,"pulse y",0);
      dimage->pulse->centBin=(int)(dimage->pulse->nBins/2);

      max=-100.0;
      tot=0.0;
      x=-1.0*(float)dimage->pulse->centBin*dimage->pRes;
      for(i=0;i<dimage->pulse->nBins;i++){
        dimage->pulse->x[i]=x;
        dimage->pulse->y[i]=(float)gaussian((double)x,(float)dimage->pSigma,0.0);
        if(dimage->pulse->y[i]>max){
          max=dimage->pulse->y[i];
          dimage->pulse->centBin=i;
        }
        tot+=dimage->pulse->y[i];
        x+=dimage->pRes;
      }
      /*normalise to cope with rounding*/
      for(i=0;i<dimage->pulse->nBins;i++){
        dimage->pulse->y[i]/=tot;
      }
    }else{  /*dirac-delta*/
      dimage->pulse->nBins=1;
      dimage->pulse->x=falloc(dimage->pulse->nBins,"pulse x",0);
      dimage->pulse->y=falloc(dimage->pulse->nBins,"pulse y",0);
      dimage->pulse->centBin=0;

      dimage->pulse->x[0]=0.0;
      dimage->pulse->y[0]=1.0;
    }
  }else{  /*read the pulse from a file*/
    readSimPulse(dimage);
  }

  return;
}/*setGediPulse*/


/*####################################*/
/*read pulse to use for simulator*/

void readSimPulse(control *dimage)
{
  int i=0,maxBin=0;
  float CofG=0,tot=0,centre=0;
  float minSep=0,max=0;
  float wThresh=0;
  float p0=0,p1=0;
  char line[400];
  char temp1[100],temp2[100];
  FILE *ipoo=NULL;

  if((ipoo=fopen(dimage->pulseFile,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",dimage->pulseFile);
    exit(1);
  }


  /*count number of bins*/
  dimage->pulse->nBins=0;
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))dimage->pulse->nBins++;

  dimage->pulse->x=falloc(dimage->pulse->nBins,"pulse x",0);
  dimage->pulse->y=falloc(dimage->pulse->nBins,"pulse y",0);

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  i=0;
  while(fgets(line,400,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s",temp1,temp2)==2){
        dimage->pulse->x[i]=atof(temp1);
        dimage->pulse->y[i]=atof(temp2);
        i++;
      }
    }
  }
  dimage->pRes=fabs(dimage->pulse->x[dimage->pulse->nBins-1]-dimage->pulse->x[0])/(float)(dimage->pulse->nBins-1);

  /*determine maximum to centre and total to normalise*/
  tot=0.0;
  CofG=0.0;
  max=-1000.0;
  for(i=0;i<dimage->pulse->nBins;i++){
    CofG+=dimage->pulse->x[i]*dimage->pulse->y[i];
    if(dimage->pulse->y[i]>max){
      max=dimage->pulse->y[i];
      centre=dimage->pulse->x[i];
      maxBin=i;
    }
    tot+=dimage->pulse->y[i];
  }
  CofG/=tot;

  /*align pulse*/
  minSep=1000.0;
  dimage->pSigma=0.0;
  for(i=0;i<dimage->pulse->nBins;i++){
    dimage->pulse->x[i]-=centre;

    if(fabs(dimage->pulse->x[i])<minSep){
      minSep=fabs(dimage->pulse->x[i]);
      dimage->pulse->centBin=i;
    }
  }

  /*pulse width*/
  wThresh=max*exp(-0.5);
  minSep=1000000.0;
  for(i=0;i<dimage->pulse->nBins;i++){
    if(i<maxBin){
      if(fabs(dimage->pulse->y[i]-wThresh)<minSep){
        minSep=fabs(dimage->pulse->y[i]-wThresh);
        p0=dimage->pulse->x[i];
      }
    }else if(i==maxBin){
      minSep=1000000.0;
    }else{
      if(fabs(dimage->pulse->y[i]-wThresh)<minSep){
        minSep=fabs(dimage->pulse->y[i]-wThresh);
        p1=dimage->pulse->x[i];
      }
    }
  }
  dimage->pSigma=fabs(p1-p0)/2.0;

  /*now normalise*/
  for(i=0;i<dimage->pulse->nBins;i++)dimage->pulse->y[i]/=tot;

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return;
}/*readSimPulse*/


/*####################################*/
/*read ASCII data*/

pCloudStruct *readAsciiData(char *inNamen)
{
  uint64_t i=0;
  pCloudStruct *data=NULL;
  char line[400],temp1[100],temp2[100];
  char temp3[100],temp4[100],temp5[100];
  char temp6[100],temp7[100],temp8[100];
  FILE *ipoo=NULL;


  if(!(data=(pCloudStruct *)calloc(1,sizeof(pCloudStruct)))){
    fprintf(stderr,"error pCloudStruct allocation.\n");
    exit(1);
  }


  if((ipoo=fopen(inNamen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",inNamen);
    exit(1);
  }

  /*count number of points*/
  data->nPoints=0;
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))data->nPoints++;

  /*allocate space*/
  data->x=dalloc(data->nPoints,"x",0);
  data->y=dalloc(data->nPoints,"y",0);
  data->z=dalloc(data->nPoints,"z",0);
  data->refl=ialloc(data->nPoints,"refl",0);


  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){ 
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  i=0;
  while(fgets(line,400,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8)==8){
        data->x[i]=atof(temp1);
        data->y[i]=atof(temp2);
        data->z[i]=atof(temp3);
        data->refl[i]=atoi(temp8);
        i++;
      }
    }
  }

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(data);
}/*readAsciiData*/


/*##############################################*/
/*set up HDF structure and write header*/

gediHDF *setUpHDF(control *dimage)
{
  int i=0;
  gediHDF *hdfData=NULL;

  /*set counter to zero*/
  dimage->hdfCount=0;

  /*allocate space*/
  if(!(hdfData=(gediHDF *)calloc(1,sizeof(gediHDF)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*header*/
  hdfData->nWaves=dimage->gNx*dimage->gNy;
  hdfData->nBins=dimage->maxBins;
  hdfData->nTypeWaves=3;
  hdfData->pSigma=dimage->pSigma;
  hdfData->fSigma=dimage->fSigma;

  /*max id label length*/
  if(dimage->useID){
    if(dimage->doGrid){
      hdfData->idLength=(int)strlen(dimage->waveID)+1+20;
    }else if(dimage->waveIDlist){
      hdfData->idLength=-1;
      for(i=0;i<hdfData->nWaves;i++){
        if(((int)strlen(dimage->waveIDlist[i])+1)>hdfData->idLength)hdfData->idLength=(int)strlen(dimage->waveIDlist[i])+1;
      }
    }else hdfData->idLength=(int)strlen(dimage->waveID)+1;
  }else hdfData->idLength=7;

  if(dimage->readPulse){
    hdfData->pRes=dimage->pRes;
    hdfData->nBins=dimage->pulse->nBins;
    hdfData->pulse=falloc(hdfData->nBins,"hdf pulse",0);
    memcpy(hdfData->pulse,dimage->pulse->y,sizeof(float)*hdfData->nBins);
  }else{
    hdfData->pulse=NULL;
    hdfData->nPbins=0;
  }

  /*allocate arrays*/
  hdfData->wave=fFalloc(hdfData->nTypeWaves,"hdf waveforms",0);
  for(i=0;i<hdfData->nTypeWaves;i++)hdfData->wave[i]=falloc(hdfData->nWaves*hdfData->nBins,"hdf waveforms",i+1);
  if(dimage->ground){
    hdfData->ground=fFalloc(hdfData->nTypeWaves,"hdf ground waveforms",0);
    for(i=0;i<hdfData->nTypeWaves;i++)hdfData->ground[i]=falloc(hdfData->nWaves*hdfData->nBins,"hdf ground waveforms",i+1);
  }
  hdfData->z0=falloc(hdfData->nWaves,"hdf z0",0);
  hdfData->zN=falloc(hdfData->nWaves,"hdf zN",0);
  hdfData->lon=dalloc(hdfData->nWaves,"hdf lon",0);
  hdfData->lat=dalloc(hdfData->nWaves,"hdf lat",0);
  hdfData->slope=falloc(hdfData->nWaves,"hdf slope",0);
  hdfData->gElev=falloc(hdfData->nWaves,"hdf gElev",0);
  hdfData->demElev=falloc(hdfData->nWaves,"hdf demElev",0);
  hdfData->beamDense=falloc(hdfData->nWaves,"hdf beamDense",0);
  hdfData->pointDense=falloc(hdfData->nWaves,"hdf pointDense",0);
  hdfData->zen=falloc(hdfData->nWaves,"hdf zen",0);
  hdfData->waveID=challoc(hdfData->nWaves*hdfData->idLength,"hdf waveID",0);

  return(hdfData);
}/*setUpHDF*/


/*##############################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  dimage->nFiles=1;
  dimage->inList=chChalloc(dimage->nFiles,"inList",0);
  dimage->inList[0]=challoc(200,"inList",0);
  strcpy(&(dimage->inList[0][0]),"/Users/dill/data/teast/maryland_play/sc_79_112_1.las");
  strcpy(dimage->outNamen,"teast.wave");
  dimage->pFWHM=15.0;   /*12 ns FWHM*/
  dimage->fWidth=5.5;  /*86% of energy within a diameter of 20-25m*/
  dimage->res=0.15;
  dimage->pRes=0.01;
  dimage->coord[0]=624366.0;
  dimage->coord[1]=3.69810*pow(10.0,6.0);

  /*switches*/
  dimage->readWave=0;
  dimage->ground=0;
  dimage->sideLobe=0;   /*no side lobes*/
  dimage->lobeAng=0.0;
  dimage->listFiles=0;
  dimage->pBuffSize=(uint64_t)200000000;
  dimage->checkCover=0;
  dimage->normCover=1;
  dimage->cleanOut=0;
  dimage->topHat=0;
  dimage->useID=0;
  dimage->readPulse=0;
  dimage->useShadow=0;
  dimage->vRes[0]=dimage->vRes[1]=dimage->vRes[2]=1.0;
  dimage->beamRad=0.165;    /*33 cm*/
  dimage->maxScanAng=1000000.0;   /*maximum scan angle*/
  dimage->polyGr=0;     /*don't fit a polynomial through the ground*/
  dimage->nnGr=0;       /*don't make a DEM from nearest neighbour*/
  dimage->overWrite=1;  /*over write any files with the same name if they exist*/
  dimage->readALSonce=0;/*read each footprint separately*/
  dimage->writeHDF=0;   /*write output as ascii*/

  /*gridding options*/
  dimage->doGrid=0;           /*gridded switch*/
  dimage->gRes=30.0;          /*grid resolution*/
  dimage->gNx=dimage->gNy=0;
  /*batch*/
  dimage->coords=NULL;       /*list of coordinates*/
  dimage->waveIDlist=NULL;   /*list of waveform IDs*/
  dimage->maxBins=1024;      /*to match LVIS*/
  dimage->gediIO.nMessages=200;

  dimage->iThresh=0.0006;
  dimage->meanN=12.0;
  dimage->doDecon=0;
  dimage->indDecon=0;
  dimage->pSigma=-1.0;  /*leave blannk for default GEDI*/
  dimage->fSigma=-1.0;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->nFiles=1;
        dimage->inList=chChalloc(dimage->nFiles,"input name list",0);
        dimage->inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=readInList(&dimage->nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-coord",6)){
        checkArguments(2,i,argc,"-coord");
        dimage->coord[0]=atof(argv[++i]);
        dimage->coord[1]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-decon",6)){
        dimage->doDecon=1;
      }else if(!strncasecmp(argv[i],"-indDecon",9)){
        dimage->indDecon=1;
        dimage->doDecon=1;
      }else if(!strncasecmp(argv[i],"-LVIS",5)){
        dimage->pSigma=0.6893;  /*two way trip*/
        dimage->fSigma=6.25;
      }else if(!strncasecmp(argv[i],"-pSigma",7)){
        checkArguments(1,i,argc,"-pSigma");
        dimage->pSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pFWHM",6)){
        checkArguments(1,i,argc,"-pFWHM");
        dimage->pFWHM=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fSigma",7)){
        checkArguments(1,i,argc,"-fSigma");
        dimage->fSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-readWave",9)){
        dimage->readWave=1;
      }else if(!strncasecmp(argv[i],"-ground",8)){
        dimage->ground=1;
        dimage->cleanOut=1;
      }else if(!strncasecmp(argv[i],"-sideLobe",9)){
        dimage->sideLobe=1;
      }else if(!strncasecmp(argv[i],"-lobeAng",8)){
        checkArguments(1,i,argc,"-lobeAng");
        dimage->lobeAng=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-listFiles",10)){
        dimage->listFiles=1;
      }else if(!strncasecmp(argv[i],"-pBuff",6)){
        checkArguments(1,i,argc,"-pBuff");
        dimage->pBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0);
      }else if(!strncasecmp(argv[i],"-noNorm",7)){
        dimage->normCover=0;
      }else if(!strncasecmp(argv[i],"-checkCover",11)){
        dimage->checkCover=1;
      }else if(!strncasecmp(argv[i],"-topHat",7)){
        dimage->topHat=1;
      }else if(!strncasecmp(argv[i],"-waveID",7)){
        checkArguments(1,i,argc,"-waveID");
        dimage->useID=1;
        strcpy(dimage->waveID,argv[++i]);
      }else if(!strncasecmp(argv[i],"-readPulse",10)){
        checkArguments(1,i,argc,"-readPulse");
        dimage->readPulse=1;
        strcpy(dimage->pulseFile,argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxScanAng",11)){
        checkArguments(1,i,argc,"-maxScanAng");
        dimage->maxScanAng=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-useShadow",10)){
        dimage->useShadow=1;
      }else if(!strncasecmp(argv[i],"-polyGround",11)){
        dimage->polyGr=1;
      }else if(!strncasecmp(argv[i],"-nnGround",99)){
        dimage->nnGr=1;
      }else if(!strncasecmp(argv[i],"-res",4)){
        checkArguments(1,i,argc,"-res");
        dimage->res=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gridBound",10)){
        checkArguments(4,i,argc,"-gridBound");
        dimage->doGrid=1;
        dimage->useID=1;
        dimage->gMinX=atof(argv[++i]);
        dimage->gMaxX=atof(argv[++i]);
        dimage->gMinY=atof(argv[++i]);
        dimage->gMaxY=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-gridStep",9)){
        checkArguments(1,i,argc,"-gridStep");
        dimage->gRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-keepOld",8)){
        dimage->overWrite=0;
      }else if(!strncasecmp(argv[i],"-listCoord",10)){
        checkArguments(1,i,argc,"-listCoord");
        dimage->readALSonce=1;
        dimage->useID=1;
        strcpy(dimage->coordList,argv[++i]);
      }else if(!strncasecmp(argv[i],"-hdf",4)){
        dimage->writeHDF=1;
      }else if(!strncasecmp(argv[i],"-maxBins",8)){
        checkArguments(1,i,argc,"-maxBins");
        dimage->maxBins=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to create GEDI waveforms from ALS las files\n#####\n\n-input name;     lasfile input filename\n-output name;    output filename\n-inList list;    input file list for multiple files\n-coord lon lat;  footprint coordinate in same system as lasfile\n-listCoord name; list of coordinates\n-gridBound minX maxX minY maxY;    make a grid of waveforms in this box\n-gridStep res;   grid step size\n-waveID id;      supply a waveID to pass to the output\n-hdf;            write output as HDF5. Best with gridded or list of coords\n-maxBins;        for HDF5, limit number of bins to save trimming\n-readPulse file; read pulse shape from a file\n-decon;          deconvolve\n-indDecon;       deconvolve individual beams\n-LVIS;           use LVIS pulse length, sigma=6.25m\n-pSigma sig;     set pulse width\n-pFWHM fhwm;     set pulse width in ns\n-fSigma sig;     set footprint width\n-res res;        range resolution to output in metres\n-readWave;       read full-waveform where available\n-ground;         split ground and canopy  points\n-sideLobe;       use side lobes\n-lobeAng ang;    lobe axis azimuth\n-topHat;         use a top hat wavefront\n-listFiles;      list files. Do not read them\n-pBuff s;        point reading buffer size in Gbytes\n-noNorm;         don't normalise for ALS density\n-checkCover;     check that the footprint is covered by ALS data. Exit if not\n-keepOld;        do not overwrite old files, if they exist\n-maxScanAng ang; maximum scan angle, degrees\n-useShadow;      account for shadowing in discrete return data through voxelisation\n-polyGround;     find mean ground elevation and slope through fitting a polynomial\n-nnGround;       find mean ground elevation and slope through nearest neighbour\n\nQuestions to svenhancock@gmail.com\n\n");
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
/*####################################*/

