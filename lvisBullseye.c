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
#include "ogr_srs_api.h"

/*#############################*/
/*# Uses simulated waveforms #*/
/*# to correlte to LVIS to   #*/
/*# coregister     2017      #*/
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
  /*input/output*/
  gediIOstruct simIO;    /*input/output structure*/
  gediIOstruct lvisIO;   /*input/output structure*/
  gediRatStruct gediRat; /*simulator options*/
  char outNamen[200];
  int nLvis;           /*number of LVIS*/

  /*options*/
  char useLvisHDF;     /*input data format switch*/
  char useLvisLGW;     /*input data format switch*/
  char useGediHDF;     /*input data format switch*/
  float offset;        /*vertical datum offset*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/
  char filtOutli;      /*filter outliers to avoid falling trees*/
  float maxZen;        /*maximum LVIS zenith angle to use*/

  /*bullseye settings*/
  float maxShift;      /*maximum distance to shift*/
  float shiftStep;     /*distance to shift steps*/

  /*bounds for subsets*/
  char useBounds;
  double minX;
  double maxX;
  double minY;
  double maxY;
  int aEPSG;           /*ALS epsg*/
  int lEPSG;           /*LVIS epsg*/
}control;


/*########################################################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct **lvis=NULL;
  dataStruct **readMultiLVIS(control *,float *);
  pCloudStruct **als=NULL;
  pCloudStruct **readMultiALS(control *,dataStruct **);
  void bullseyeCorrel(dataStruct **,pCloudStruct **,control *);

  /*read command Line*/
  dimage=readCommands(argc,argv);

  /*read LVIS data*/
  lvis=readMultiLVIS(dimage,&dimage->simIO.res);

  /*read ALS data*/
  als=readMultiALS(dimage,lvis);

  /*loop over, simulating*/
  bullseyeCorrel(lvis,als,dimage);

  /*tidy up*/
  if(als){
    for(i=0;i<dimage->simIO.nFiles;i++)als[i]=tidyPointCloud(als[i]);
    TIDY(als);
  }
  lvis=tidyAsciiStruct(lvis,dimage->nLvis);
  if(dimage){
    TTIDY((void **)dimage->lvisIO.inList,dimage->lvisIO.nFiles);
    TTIDY((void **)dimage->simIO.inList,dimage->simIO.nFiles);
    if(dimage->lvisIO.den){
      TTIDY((void **)dimage->lvisIO.den->pulse,2);
      TIDY(dimage->lvisIO.den->matchPulse);
      TIDY(dimage->lvisIO.den->hardPulse);
      TIDY(dimage->lvisIO.den);
    }
    if(dimage->lvisIO.gFit){
      TTIDY((void **)dimage->lvisIO.gFit->pulse,2);
      TIDY(dimage->lvisIO.gFit);
    }
    if(dimage->simIO.den){
      TTIDY((void **)dimage->simIO.den->pulse,2);
      TIDY(dimage->simIO.den->matchPulse);
      TIDY(dimage->simIO.den->hardPulse);
      TIDY(dimage->simIO.den);
    }
    if(dimage->simIO.gFit){
      TTIDY((void **)dimage->simIO.gFit->pulse,2);
      TIDY(dimage->simIO.gFit);
    }
    TTIDY((void **)dimage->gediRat.coords,dimage->gediRat.gNx);
    TTIDY((void **)dimage->gediRat.waveIDlist,dimage->gediRat.gNx);
    TIDY(dimage->gediRat.nGrid);
    dimage->gediRat.octree=tidyOctree(dimage->gediRat.octree);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################################*/
/*simulate and calculate correlation*/

void bullseyeCorrel(dataStruct **lvis,pCloudStruct **als,control *dimage)
{
  int i=0,j=0,k=0;
  int nX=0,contN=0;
  int nTypeWaves=0;
  double xOff=0,yOff=0;
  double **coords=NULL;
  double **shiftPrints(double **,double,double,int);
  float **denoiseAllLvis(dataStruct **,control *);
  float **denoised=NULL;
  float **correl=NULL;
  float *waveCorrel(waveStruct *,float *,dataStruct *,gediIOstruct *);
  void writeCorrelStats(float **,int,int,FILE *,double,double,control *);
  waveStruct *waves=NULL;
  FILE *opoo=NULL;


  /*how mamy types of simuation methods*/
  nTypeWaves=(int)(dimage->simIO.useCount+dimage->simIO.useFrac+dimage->simIO.useInt);

  /*open ourput file*/
  if((opoo=fopen(dimage->outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->outNamen);
    exit(1);
  }
  fprintf(opoo,"# 1 xOff, 2 yOff");
  if(dimage->simIO.useInt)fprintf(opoo,", 3 correlInt, 4 stdev, 5 deltaCofGint, 6 numb");
  if(dimage->simIO.useCount)fprintf(opoo,", %d correlCount, %d stdev, %d deltaCofGcount, %d numb",3+4*dimage->simIO.useInt,3+4*dimage->simIO.useInt+1,3+4*dimage->simIO.useInt+2,3+4*dimage->simIO.useInt+3);
  if(dimage->simIO.useFrac)fprintf(opoo,", %d correlFrac, %d stdev, %d deltaCofGfrac, %d numb",3+4*(dimage->simIO.useInt+dimage->simIO.useCount),3+4*(dimage->simIO.useInt+dimage->simIO.useCount)+1\
                                                                                              ,3+4*(dimage->simIO.useInt+dimage->simIO.useCount)+2,3+4*(dimage->simIO.useInt+dimage->simIO.useCount)+3);
  fprintf(opoo,"\n");

  /*set up pulse*/
  setGediPulse(&dimage->simIO,&dimage->gediRat);

  /*number of steps*/
  nX=(int)(2.0*dimage->maxShift/dimage->shiftStep+1);

  /*save original coordinates*/
  coords=dimage->gediRat.coords;
  dimage->gediRat.coords=NULL;

  /*denoise LVIS*/
  denoised=denoiseAllLvis(lvis,dimage);

  /*loop over shifts*/
  for(i=0;i<nX;i++){
    xOff=(double)(i-nX/2)*(double)dimage->shiftStep;
    for(j=0;j<nX;j++){
      yOff=(double)(j-nX/2)*(double)dimage->shiftStep;
      fprintf(stdout,"Testing x %f y %f\n",xOff,yOff);
      /*shift prints*/
      dimage->gediRat.coords=shiftPrints(coords,xOff,yOff,dimage->gediRat.gNx);
      correl=fFalloc(dimage->gediRat.gNx,"correl",0);

      /*loop over footprints*/
      contN=0;
      for(k=0;k<dimage->gediRat.gNx;k++){
        if(lvis[k]->zen>dimage->maxZen)continue;
        /*set one coordinate*/
        updateGediCoord(&dimage->gediRat,k,0);

        /*simulate waves*/
        setGediFootprint(&dimage->gediRat,&dimage->simIO);
        waves=makeGediWaves(&dimage->gediRat,&dimage->simIO,als);

        /*calculate correlation*/
        if(dimage->gediRat.useFootprint){
          correl[contN]=waveCorrel(waves,denoised[k],lvis[k],&dimage->simIO);
        }else correl[contN]=NULL;
     
        /*tidy up*/
        if(waves){
          TTIDY((void **)waves->wave,waves->nWaves);
          TTIDY((void **)waves->canopy,nTypeWaves);
          TTIDY((void **)waves->ground,nTypeWaves);
          TIDY(waves);
        }
        TIDY(dimage->gediRat.lobe);
        TIDY(dimage->gediRat.nGrid);
        contN++;
      }/*footprint loop*/

      /*output results*/
      writeCorrelStats(correl,contN,nTypeWaves,opoo,xOff,yOff,dimage);

      /*tidy up*/
      TTIDY((void **)dimage->gediRat.coords,dimage->gediRat.gNx);
      TTIDY((void **)correl,dimage->gediRat.gNx);
    }/*y loop*/
  }/*x loop*/


  /*tidy up*/
  dimage->gediRat.coords=coords;
  coords=NULL;
  TTIDY((void **)denoised,dimage->nLvis);
  if(dimage->simIO.pulse){
    TIDY(dimage->simIO.pulse->y);
    TIDY(dimage->simIO.pulse->x);
    TIDY(dimage->simIO.pulse);
  }
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",dimage->outNamen);
  return;
}/*bullseyeCorrel*/


/*####################################################*/
/*correlation stats and write*/

void writeCorrelStats(float **correl,int numb,int nTypes,FILE *opoo,double xOff,double yOff,control *dimage)
{
  int i=0,k=0,nUsed=0;
  int usedNew=0;
  float mean=0,stdev=0;
  float newMean=0,meanCofG=0;
  float thresh=0;

  fprintf(opoo,"%f %f",xOff,yOff);

  /*loop over types*/
  for(k=0;k<nTypes;k++){
    mean=stdev=0.0;
    nUsed=0;
    for(i=0;i<numb;i++){
      if(correl[i]){
        mean+=correl[i][2*k];
        nUsed++;
      }
    }
    mean/=(float)nUsed;
    meanCofG/=(float)nUsed;
    for(i=0;i<numb;i++){
      if(correl[i]){
        stdev+=pow(correl[i][2*k]-mean,2.0);
      }
    }
    stdev=sqrt(stdev/(float)nUsed);

    /*check for outliers*/
    if(dimage->filtOutli)thresh=2.5*stdev;
    else                 thresh=1000000000.0;

    usedNew=0;
    newMean=meanCofG=0.0;
    for(i=0;i<numb;i++){
      if(correl[i]){
        if((mean-correl[i][2*k])<thresh){
          newMean+=correl[i][2*k];
          meanCofG+=correl[i][2*k+1];
          usedNew++;
        }
      }
    }
    newMean/=(float)usedNew;
    meanCofG/=(float)usedNew;
    stdev=0.0;
    for(i=0;i<numb;i++){
      if(correl[i]){
        if((mean-=correl[i][2*k])<thresh){
          stdev+=pow(correl[i][2*k]-newMean,2.0);
        }
      }
    }
    stdev=sqrt(stdev/(float)usedNew);
    fprintf(opoo," %f %f %f %d",newMean,stdev,meanCofG,usedNew);
  }
  fprintf(opoo," %d\n",nUsed);

  return;
}/*writeCorrelStats*/


/*####################################################*/
/*calculate correlation*/

float *waveCorrel(waveStruct *sim,float *truth,dataStruct *lvis,gediIOstruct *simIO)
{
  int i=0,j=0,k=0,bin=0;
  int numb=0;
  float *correl=NULL;
  float totS=0,totL=0;
  float CofGl=0,CofGs=0;
  float z=0,thresh=0;
  float sLx=0,eLx=0;
  float sSx=0,eSx=0;
  float startX=0,endX=0;
  float sumProd=0,minSepSq=0;
  float sepSq=0;
  float *smooSim=NULL;
  float meanL=0,meanS=0;
  float stdevL=0,stdevS=0;

  /*allocate space fgor correlation and CofG shift*/
  correl=falloc(2*sim->nWaves,"correlation",0);

  /*total energies*/
  totL=CofGl=0.0;
  for(i=0;i<lvis->nBins;i++){
    totL+=truth[i];
    CofGl+=truth[i]*(float)lvis->z[i];
  }
  CofGl/=totL;

  /*find lvis bounds*/
  thresh=0.0001*totL;
  for(i=0;i<lvis->nBins;i++){
    if(truth[i]>thresh){
      eLx=(float)lvis->z[i];
      break;
    }
  }
  for(i=lvis->nBins-1;i>=0;i--){
    if(truth[i]>thresh){
      sLx=(float)lvis->z[i];
      break;
    }
  }

  /*loop over three wave types*/
  for(k=0;k<sim->nWaves;k++){
    smooSim=processFloWave(sim->wave[k],sim->nBins,simIO->den,1.0);


    totS=CofGs=0.0;
    for(i=0;i<sim->nBins;i++){
      totS+=smooSim[i];
      z=(float)sim->maxZ-(float)i*simIO->res;
      CofGs+=smooSim[i]*z;
    }

    thresh=0.0001*totS;
    CofGs/=totS;
    correl[2*k+1]=CofGl-CofGs;

    /*find sim bounds*/
    for(i=0;i<sim->nBins;i++){
      if(smooSim[i]>thresh){
        z=(float)sim->maxZ-(float)i*simIO->res;
        eSx=z;
        break;
      } 
    }
    for(i=sim->nBins-1;i>=0;i--){
      if(smooSim[i]>thresh){
        z=(float)sim->maxZ-(float)i*simIO->res;
        sSx=z;
        break;
      }
    }

    startX=(sSx<sLx)?sSx:sLx;
    endX=(eSx>eLx)?eSx:eLx;
    numb=(int)(fabs(endX-startX)/simIO->res);

    /*means*/
    meanS=meanL=0.0;
    for(i=0;i<lvis->nBins;i++)if((lvis->z[i]>=startX)&&(lvis->z[i]<=endX))meanL+=truth[i];
    for(i=0;i<sim->nBins;i++){
      z=(float)sim->maxZ-(float)i*simIO->res;
      if((z>=startX)&&(z<=endX))meanS+=smooSim[i];
    }
    meanL/=(float)numb;
    meanS/=(float)numb;
    /*stdev*/
    stdevS=stdevL=0.0;
    for(i=0;i<lvis->nBins;i++)if((lvis->z[i]>=startX)&&(lvis->z[i]<=endX))stdevL+=(truth[i]-meanL)*(truth[i]-meanL);
    for(i=0;i<sim->nBins;i++){
      z=(float)sim->maxZ-(float)i*simIO->res;
      if((z>=startX)&&(z<=endX))stdevS+=(smooSim[i]-meanS)*(smooSim[i]-meanS);
    }
    stdevL=sqrt(stdevL/(float)numb);
    stdevS=sqrt(stdevS/(float)numb);

    /*shared variance*/
    sumProd=0.0;
    for(i=lvis->nBins-1;i>=0;i--){
      if(((float)lvis->z[i]<startX)||((float)lvis->z[i]>endX))continue;
      minSepSq=100000.0;
      for(j=0;j<sim->nBins;j++){
        z=(float)sim->maxZ-(float)j*simIO->res;
        if((z<startX)||(z>endX))continue;
        sepSq=pow((z-(float)lvis->z[i]),2.0);
        if(sepSq<minSepSq){
          minSepSq=sepSq;
          bin=j;
        }
      }
      sumProd+=(truth[i]-meanL)*(smooSim[bin]-meanS);
    }
    correl[2*k]=(sumProd/(float)numb)/(stdevL*stdevS);
    TIDY(smooSim);
  }/*wave type loop*/

  return(correl);
}/*waveCorrel*/


/*####################################################*/
/*denoiseall LVIS data*/

float **denoiseAllLvis(dataStruct **lvis,control *dimage)
{
  int i=0,j=0;
  float **denoised=NULL,tot=0;

  /*allocate space*/
  denoised=fFalloc(dimage->nLvis,"denoised waveforms",0);

  /*loop over LVIS footprints*/
  for(i=0;i<dimage->nLvis;i++){
    denoised[i]=processFloWave(lvis[i]->wave[0],lvis[i]->nBins,dimage->lvisIO.den,1.0);
    /*rescale*/
    tot=0.0;
    for(j=0;j<lvis[i]->nBins;j++)tot+=denoised[i][j];
    for(j=0;j<lvis[i]->nBins;j++)denoised[i][j]/=tot;
  }/*LVIS footprint loop*/

  return(denoised);
}/*denoiseAllLvis*/


/*####################################################*/
/*shift coordinate centres*/

double **shiftPrints(double **coords,double xOff,double yOff,int numb)
{
  int i=0;
  double **newCoords=NULL;

  /*allocate space*/
  newCoords=dDalloc(numb,"newCoords",0);

  /*loop overall points*/
  for(i=0;i<numb;i++){
    newCoords[i]=dalloc(2,"newCoords",i+1);
    newCoords[i][0]=coords[i][0]+xOff;
    newCoords[i][1]=coords[i][1]+yOff;
  }

  return(newCoords);
}/*shiftPrints*/


/*####################################################*/
/*read multiple ALS files and save relevant data*/

pCloudStruct **readMultiALS(control *dimage,dataStruct **lvis)
{
  int i=0;
  pCloudStruct **als=NULL;
  lasFile *las=NULL;
  void copyLvisCoords(gediRatStruct *,dataStruct **,int,int,int);


  /*determine bounds*/
  copyLvisCoords(&dimage->gediRat,lvis,dimage->nLvis,dimage->aEPSG,dimage->lEPSG);
  setGediGrid(&dimage->simIO,&dimage->gediRat);
 /*account for jittering*/
  dimage->gediRat.globMinX-=(double)dimage->maxShift;
  dimage->gediRat.globMaxX+=(double)dimage->maxShift;
  dimage->gediRat.globMinY-=(double)dimage->maxShift;
  dimage->gediRat.globMaxY+=(double)dimage->maxShift;

  /*allocate space*/
  if(!(als=(pCloudStruct **)calloc(dimage->simIO.nFiles,sizeof(pCloudStruct *)))){
    fprintf(stderr,"error in lvis data allocation.\n");
    exit(1);
  } 

  /*loop over ALS files*/
  for(i=0;i<dimage->simIO.nFiles;i++){
    /*read header*/
    las=readLasHead(dimage->simIO.inList[i],dimage->pBuffSize);

    /*read data*/
    als[i]=readALSdata(las,&dimage->gediRat,i);
    if(als[i]->nPoints>0)fprintf(stdout,"%d ALS\n",als[i]->nPoints);

    /*tidy up*/
    las=tidyLasFile(las);
  }/*ALS file loop*/

  return(als);
}/*readMultiALS*/


/*####################################################*/
/*copy LVIS coords into ALS structure*/

void copyLvisCoords(gediRatStruct *gediRat,dataStruct **lvis,int nLvis,int aEPSG,int lEPSG)
{
  int i=0;
  double *x=NULL,*y=NULL,*z=NULL;
  OGRCoordinateTransformationH hTransform;
  OGRSpatialReferenceH hSourceSRS,hTargetSRS;
  OGRErr err;

  gediRat->gNx=nLvis;
  gediRat->gNy=1;
  gediRat->coords=dDalloc(gediRat->gNx,"coord list",0);
  gediRat->waveIDlist=NULL;


  if(aEPSG!=lEPSG){
    x=dalloc(nLvis,"x",0);
    y=dalloc(nLvis,"y",0);
    z=dalloc(nLvis,"z",0);
    for(i=0;i<nLvis;i++){
      x[i]=lvis[i]->lon;
      y[i]=lvis[i]->lat;
      z[i]=0.0;
    }

    hSourceSRS=OSRNewSpatialReference(NULL);
    hTargetSRS=OSRNewSpatialReference(NULL);
    err=OSRImportFromEPSG(hTargetSRS,aEPSG);
    err=OSRImportFromEPSG(hSourceSRS,lEPSG);
    hTransform=OCTNewCoordinateTransformation(hSourceSRS,hTargetSRS);
    OCTTransform(hTransform,nLvis,x,y,z);
    OCTDestroyCoordinateTransformation(hTransform);
    OSRDestroySpatialReference(hSourceSRS);
    OSRDestroySpatialReference(hTargetSRS);

    for(i=0;i<nLvis;i++){
      gediRat->coords[i]=dalloc(2,"coords",i+1);
      gediRat->coords[i][0]=x[i];
      gediRat->coords[i][1]=y[i];
    }
  }else{
    for(i=0;i<nLvis;i++){
      gediRat->coords[i]=dalloc(2,"coords",i+1);
      gediRat->coords[i][0]=lvis[i]->lon;
      gediRat->coords[i][1]=lvis[i]->lat;
    }
  }

  TIDY(x);
  TIDY(y);
  TIDY(z);
  return;
}/*copyLvisCoords*/


/*####################################################*/
/*read multiple LVIS files and save relevant data*/

dataStruct **readMultiLVIS(control *dimage,float *res)
{
  int i=0,j=0;
  dataStruct **lvis=NULL;
  dataStruct **copyLVIShdf(lvisHDF *,dataStruct **,control *,double *);
  dataStruct **copyLVISlgw(char *,dataStruct **,control *,double *);
  dataStruct **copyGEDIhdf(gediHDF *,dataStruct **,control *,double *);
  double bounds[4],offset=0;
  lvisHDF *hdf=NULL;        /*LVIS HDF5 structure*/
  gediHDF *simHDF=NULL;     /*GEDI HDF5 structure*/
  void reprojectBounds(control *,double *);


  dimage->nLvis=0;

  /*reproject bounds*/
  reprojectBounds(dimage,bounds);

  /*loop over lvus files*/
  for(i=0;i<dimage->lvisIO.nFiles;i++){
    if(dimage->useLvisHDF){
      /*read HDF5*/
      hdf=readLVIShdf(dimage->lvisIO.inList[i]);
      /*unpack*/
      lvis=copyLVIShdf(hdf,lvis,dimage,bounds);
      /*tidy up*/
      hdf=tidyLVISstruct(hdf);
    }else if(dimage->useLvisLGW){
      /*unpack*/
      lvis=copyLVISlgw(dimage->lvisIO.inList[i],lvis,dimage,bounds);
    }else if(dimage->useGediHDF){
      /*read HDF*/
      simHDF=readGediHDF(dimage->lvisIO.inList[i],&dimage->lvisIO);
      /*unpack simulated HDF data*/
      lvis=copyGEDIhdf(simHDF,lvis,dimage,bounds);
      /*tidy up*/
      simHDF=tidyGediHDF(simHDF);
    }
  }/*file loop*/

  /*offset vertical datum*/
  if(fabs(dimage->offset)>0.001){
    offset=(double)dimage->offset;
    for(i=0;i<dimage->nLvis;i++){
      for(j=0;j<lvis[i]->nBins;j++)lvis[i]->z[j]+=offset;
    }
  }

  /*res for simulations*/
  *res=0.0;
  for(i=0;i<dimage->nLvis;i++)*res+=lvis[i]->res;
  *res/=(float)dimage->nLvis;


  fprintf(stdout,"Found %d LVIS\n",dimage->nLvis);
  return(lvis);
}/*readMultiLVIS*/


/*####################################################*/
/*reproject the bounds*/

void reprojectBounds(control *dimage,double *bounds)
{
  double *x=NULL,*y=NULL,*z=NULL;
  OGRCoordinateTransformationH hTransform;
  OGRSpatialReferenceH hSourceSRS,hTargetSRS;
  OGRErr err;

  if(dimage->aEPSG!=dimage->lEPSG){
    x=dalloc(2,"x trans",9);
    y=dalloc(2,"y trans",9);
    z=dalloc(2,"z trans",9);

    x[0]=dimage->minX;
    x[1]=dimage->maxX;
    y[0]=dimage->minY;
    y[1]=dimage->maxY;
    z[0]=z[1]=0.0;

    hSourceSRS=OSRNewSpatialReference(NULL);
    hTargetSRS=OSRNewSpatialReference(NULL);
    err=OSRImportFromEPSG(hTargetSRS,dimage->lEPSG);
    err=OSRImportFromEPSG(hSourceSRS,dimage->aEPSG);
    hTransform=OCTNewCoordinateTransformation(hSourceSRS,hTargetSRS);
    OCTTransform(hTransform,2,x,y,z);
    OCTDestroyCoordinateTransformation(hTransform);
    OSRDestroySpatialReference(hSourceSRS);
    OSRDestroySpatialReference(hTargetSRS);

    bounds[0]=x[0];
    bounds[1]=y[0];
    bounds[2]=x[1];
    bounds[3]=y[1];
  }else{
    bounds[0]=dimage->minX;
    bounds[1]=dimage->minY;
    bounds[2]=dimage->maxX;
    bounds[3]=dimage->maxY;
  }

  TIDY(x);
  TIDY(y);
  TIDY(z);
  return;
}/*reprojectBounds*/


/*####################################################*/
/*unpack relevant data from LGW*/

dataStruct **copyLVISlgw(char *namen,dataStruct **lvis,control *dimage,double *bounds)
{
  int i=0,nNew=0;
  int oldN=0;
  double x=0,y=0;
  lvisLGWstruct lgw;
  dataStruct *tempData=NULL;

  /*as this is overwritten later*/
  oldN=dimage->lvisIO.nFiles;

  /*count number of new files*/
  nNew=0;
  lgw.data=NULL;
  lgw.ipoo=NULL;
  tempData=readBinaryLVIS(namen,&lgw,0,&dimage->lvisIO);
  for(i=0;i<lgw.nWaves;i++){
    x=(lgw.data[i].lon0+lgw.data[i].lon431)/2.0;
    y=(lgw.data[i].lat0+lgw.data[i].lat431)/2.0;
    if((dimage->lEPSG==4326)&&(x>180.0))x-=360.0;

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3]))nNew++;
  }
  if(tempData){
    TTIDY((void **)tempData->wave,tempData->nWaveTypes);
    TTIDY((void **)tempData->ground,tempData->nWaveTypes);
    TIDY(tempData->noised);
    TIDY(tempData->totE);
    TIDY(tempData->z);
    TIDY(tempData);
  }


  /*allocate space*/
  if(lvis==NULL){
    if(!(lvis=(dataStruct **)calloc(nNew,sizeof(dataStruct *)))){
      fprintf(stderr,"error in lvis data allocation.\n");
      exit(1);
    }
  }else{
    if(!(lvis=(dataStruct **)realloc(lvis,(dimage->nLvis+nNew)*sizeof(dataStruct *)))){
      fprintf(stderr,"Balls\n");
      exit(1);
    }
  }

  /*copy data*/
  nNew=0;
  for(i=0;i<lgw.nWaves;i++){
    x=(lgw.data[i].lon0+lgw.data[i].lon431)/2.0;
    y=(lgw.data[i].lat0+lgw.data[i].lat431)/2.0;
    if((dimage->lEPSG==4326)&&(x>180.0))x-=360.0;

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3])){
      lvis[nNew+dimage->nLvis]=readBinaryLVIS(namen,&lgw,i,&dimage->lvisIO);
      nNew++;
    }
  }

  /*keep count*/
  dimage->nLvis+=nNew;
  dimage->lvisIO.nFiles=oldN;
  TIDY(lgw.data);
  if(lgw.ipoo){
    fclose(lgw.ipoo);
    lgw.ipoo=NULL;
  }
  return(lvis);
}/*copyLVISlgw*/


/*####################################################*/
/*copy data from simulated HDF5*/

dataStruct **copyGEDIhdf(gediHDF *hdf,dataStruct **lvis,control *dimage,double *bounds)
{
  int i=0,nNew=0;
  double x=0,y=0;


  /*count number within*/
  nNew=0;
  for(i=0;i<hdf->nWaves;i++){
    x=hdf->lon[i];
    y=hdf->lat[i];
    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3]))nNew++;
  }/*point loop*/

  /*allocate space*/
  if(lvis==NULL){
    if(!(lvis=(dataStruct **)calloc(nNew,sizeof(dataStruct *)))){
      fprintf(stderr,"error in lvis data allocation.\n");
      exit(1);
    }
  }else{
    if(!(lvis=(dataStruct **)realloc(lvis,(dimage->nLvis+nNew)*sizeof(dataStruct *)))){
      fprintf(stderr,"Balls\n");
      exit(1);
    }
  }

  /*copy data*/
  nNew=0;
  for(i=0;i<hdf->nWaves;i++){
    x=hdf->lon[i];
    y=hdf->lat[i];

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3])){
      lvis[nNew+dimage->nLvis]=unpackHDFgedi(NULL,&dimage->lvisIO,&hdf,i);
      nNew++;
    }
  }

  /*keep count*/
  dimage->nLvis+=nNew;

  return(lvis);
}/*copyGEDIhdf*/


/*####################################################*/
/*unpack relevant data from HDF*/

dataStruct **copyLVIShdf(lvisHDF *hdf,dataStruct **lvis,control *dimage,double *bounds)
{
  int i=0,nNew=0;
  double x=0,y=0;

  /*count number within*/
  nNew=0;
  for(i=0;i<hdf->nWaves;i++){
    x=(hdf->lon0[i]+hdf->lon1023[i])/2.0;
    y=(hdf->lat0[i]+hdf->lat1023[i])/2.0;

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3]))nNew++;
  }/*wave loop*/
  if(!nNew)return(lvis);

  /*allocate space*/
  if(lvis==NULL){
    if(!(lvis=(dataStruct **)calloc(nNew,sizeof(dataStruct *)))){
      fprintf(stderr,"error in lvis data allocation.\n");
      exit(1);
    }
  }else{
    if(!(lvis=(dataStruct **)realloc(lvis,(dimage->nLvis+nNew)*sizeof(dataStruct *)))){
      fprintf(stderr,"Balls\n");
      exit(1);
    }
  }

  /*copy data*/
  nNew=0;
  for(i=0;i<hdf->nWaves;i++){
    x=(hdf->lon0[i]+hdf->lon1023[i])/2.0;
    y=(hdf->lat0[i]+hdf->lat1023[i])/2.0;

    if((x>=bounds[0])&&(y>=bounds[1])&&(x<=bounds[2])&&(y<=bounds[3])){
      lvis[nNew+dimage->nLvis]=unpackHDFlvis(NULL,&hdf,&dimage->lvisIO,i);
      nNew++;
    }
  }

  /*keep count*/
  dimage->nLvis+=nNew;

  return(lvis);
}/*copyLVIShdf*/


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
  if(!(dimage->simIO.den=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  dimage->simIO.gFit=NULL;
  if(!(dimage->lvisIO.den=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  dimage->lvisIO.gFit=NULL;


  /*defaults*/
  strcpy(dimage->outNamen,"teast.correl");
  dimage->useLvisHDF=1;
  dimage->useLvisLGW=0;
  dimage->useGediHDF=0;
  dimage->aEPSG=32632;
  dimage->lEPSG=4326;
  dimage->offset=0.0;
  dimage->pBuffSize=(uint64_t)200000000;
  dimage->filtOutli=1;    /*filter outliers*/
  dimage->maxZen=100000.0;

  /*octree*/
  dimage->gediRat.useOctree=1;
  dimage->gediRat.octLevels=0;
  dimage->gediRat.nOctTop=30;
  dimage->gediRat.octree=NULL;

  /*LVIS params for sim*/
  dimage->simIO.fSigma=4.31;
  dimage->simIO.pSigma=0.6893;
  dimage->gediRat.iThresh=0.002;
  dimage->simIO.useCount=1;
  dimage->simIO.useInt=0;
  dimage->simIO.useFrac=0;

  /*waveform reading settings*/
  dimage->lvisIO.useInt=dimage->lvisIO.useFrac=0;
  dimage->lvisIO.useCount=1;
  dimage->lvisIO.ground=0;
  dimage->lvisIO.nMessages=200;

  /*LVIS denoising*/
  setDenoiseDefault(dimage->lvisIO.den);
  dimage->lvisIO.den->varNoise=1;
  dimage->lvisIO.den->statsLen=12.0;
  dimage->lvisIO.den->noiseTrack=1;
  dimage->lvisIO.den->threshScale=4.0;

  /*bullseyte params*/
  dimage->maxShift=10.0;      /*maximum distance to shift*/
  dimage->shiftStep=1.0;     /*distance to shift steps*/

  /*simulation settings*/
  dimage->gediRat.readALSonce=1;    /*read all ALS data once*/
  dimage->gediRat.readWave=0;       /*do not read waveform switch*/
  dimage->gediRat.useShadow=0;      /*do not account for shadowing through voxelisation*/
  dimage->gediRat.maxScanAng=1000000000.0;    /*maximum scan angle*/
  dimage->gediRat.coords=NULL;     /*list of coordinates*/
  dimage->gediRat.readWave=0;
  dimage->simIO.ground=0;
  dimage->gediRat.sideLobe=0;   /*no side lobes*/
  dimage->gediRat.lobeAng=0.0;
  dimage->gediRat.checkCover=1;
  dimage->gediRat.normCover=1;
  dimage->gediRat.cleanOut=0;
  dimage->gediRat.topHat=0;
  dimage->simIO.readPulse=0;
  dimage->gediRat.useShadow=0;
  dimage->gediRat.maxScanAng=1000000.0;   /*maximum scan angle*/
  dimage->gediRat.coords=NULL;       /*list of coordinates*/
  dimage->gediRat.waveIDlist=NULL;   /*list of waveform IDs*/
  dimage->simIO.nMessages=200;
  dimage->gediRat.doDecon=0;
  dimage->gediRat.indDecon=0;
  dimage->simIO.pRes=0.01;
  setDenoiseDefault(dimage->simIO.den);
  dimage->simIO.den->meanN=0.0;
  dimage->simIO.den->thresh=0.0;
  dimage->simIO.den->minWidth=0.0;


  /*read the command line*/
  for(i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-lvis",5)){
        checkArguments(1,i,argc,"-lvis");
        TTIDY((void **)dimage->lvisIO.inList,dimage->lvisIO.nFiles);
        dimage->lvisIO.inList=NULL;
        dimage->lvisIO.nFiles=1;
        dimage->lvisIO.inList=chChalloc(dimage->lvisIO.nFiles,"input name list",0);
        dimage->lvisIO.inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->lvisIO.inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-listLvis",9)){
        checkArguments(1,i,argc,"-listLvis");
        TTIDY((void **)dimage->lvisIO.inList,dimage->lvisIO.nFiles);
        dimage->lvisIO.inList=readInList(&dimage->lvisIO.nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-als",4)){
        checkArguments(1,i,argc,"-als");
        TTIDY((void **)dimage->simIO.inList,dimage->simIO.nFiles);
        dimage->simIO.inList=NULL;
        dimage->simIO.nFiles=1;
        dimage->simIO.inList=chChalloc(dimage->simIO.nFiles,"input name list",0);
        dimage->simIO.inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->simIO.inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-listAls",8)){
        checkArguments(1,i,argc,"-listAls");
        TTIDY((void **)dimage->simIO.inList,dimage->simIO.nFiles);
        dimage->simIO.inList=readInList(&dimage->simIO.nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(4,i,argc,"-bounds");
        dimage->useBounds=1;
        dimage->minX=atof(argv[++i]);
        dimage->minY=atof(argv[++i]);
        dimage->maxX=atof(argv[++i]);
        dimage->maxY=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-lgw",4)){
        dimage->useLvisHDF=0;
        dimage->useLvisLGW=1;
        dimage->useGediHDF=0;
      }else if(!strncasecmp(argv[i],"-readHDFgedi",11)){
        dimage->useLvisHDF=0;
        dimage->useLvisLGW=0;
        dimage->useGediHDF=1;
      }else if(!strncasecmp(argv[i],"-aEPSG",6)){
        checkArguments(1,i,argc,"-aEPSG");
        dimage->aEPSG=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-lEPSG",6)){
        checkArguments(1,i,argc,"-lEPSG");
        dimage->lEPSG=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pBuff",6)){
        checkArguments(1,i,argc,"-pBuff");
        dimage->pBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0);
      }else if(!strncasecmp(argv[i],"-readPulse",10)){
        checkArguments(1,i,argc,"-readPulse");
        dimage->simIO.readPulse=1;
        strcpy(dimage->simIO.pulseFile,argv[++i]);
      }else if(!strncasecmp(argv[i],"-pSigma",7)){
        checkArguments(1,i,argc,"-pSigma");
        dimage->simIO.pSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fSigma",7)){
        checkArguments(1,i,argc,"-fSigma");
        dimage->simIO.fSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-smooth",7)){
        checkArguments(1,i,argc,"-smooth");
        dimage->lvisIO.den->sWidth=dimage->simIO.den->sWidth=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxShift",9)){
        checkArguments(1,i,argc,"-maxShift");
        dimage->maxShift=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-step",5)){
        checkArguments(1,i,argc,"-step");
        dimage->shiftStep=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-offset",8)){
        checkArguments(1,i,argc,"-offset");
        dimage->offset=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noNorm",7)){
        dimage->gediRat.normCover=0;
      }else if(!strncasecmp(argv[i],"-noFilt",7)){
        dimage->filtOutli=0;
      }else if(!strncasecmp(argv[i],"-noOctree",9)){
        dimage->gediRat.useOctree=0;
      }else if(!strncasecmp(argv[i],"-octLevels",9)){
        checkArguments(1,i,argc,"-octLevels");
        dimage->gediRat.octLevels=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-nOctPix",8)){
        checkArguments(1,i,argc,"-nOctPix");
        dimage->gediRat.nOctTop=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxZen",7)){
        checkArguments(1,i,argc,"-maxZen");
        dimage->maxZen=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-allSimMeth",11)){
        dimage->simIO.useCount=dimage->simIO.useInt=dimage->simIO.useFrac=1;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to calculate GEDI waveform metrics\n#####\n\n-output name;     output filename\n-listAls list;    input file list for multiple als files\n-als file;        input als file\n-lvis file;       single input LVIS file\n-listLvis file;   list of multiple LVIS files\n-lgw;             LVIS is in lgw (default is LVIS hdf5)\n-readHDFgedi;     read GEDI HDF5 input (default is LVIS hdf5)\n-lEPSG epsg;      LVIS projection\n-aEPSG epsg;      ALS projection\n-pSigma x;        pulse length, sigma in metres\n-fSigma x;        footprint width, sigma in metres\n-readPulse file;  pulse shape\n-smooth sig;      smooth both waves before comparing\n-maxShift x;      distance to search over\n-step x;          steps to take\n-offset x;        vertical datum offset\n-bounds minX minY maxX maxY;    bounds to use, in ALS projection\n-noNorm;          don't correct sims for ALS densiy variations\n-noFilt;          don't filter outliers from correlation\n-allSimMeth;      use all simulation methods\n\n# Octree\n-noOctree;      do not use an octree\n-octLevels n;   number of octree levels to use\n-nOctPix n;     number of octree pixels along a side for the top level\n-maxZen zen;     maximum zenith angle to use, degrees\n\n");
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
/*########################################################################*/

