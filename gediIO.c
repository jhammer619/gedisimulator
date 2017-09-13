#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "hdf5.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "gediIO.h"


/*tolerances*/
#define TOL 0.00001



/*####################################################*/
/*tidy data structure*/

dataStruct **tidyAsciiStruct(dataStruct **data,int nFiles)
{
  int i=0;

  if(data){
    for(i=0;i<nFiles;i++){
      TIDY(data[i]->wave);
      TIDY(data[i]->ground);
      TIDY(data[i]->noised);
      TIDY(data[i]->z);
      TIDY(data[i]);
    }
    TIDY(data);
  }

  return(data);
}/*tidyAsciiStruct*/


/*####################################################*/
/*read ASCII data*/

dataStruct *readASCIIdata(char *namen,gediIOstruct *gediIO)
{
  int i=0;
  dataStruct *data=NULL;
  char line[400],temp1[100],temp2[100];
  char temp3[100],temp4[100],temp5[100];
  char temp6[100],temp7[100],temp8[100];
  char temp9[100],temp10[100];
  FILE *ipoo=NULL;

  /*open input*/
  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  data->pSigma=-1.0;    /*nonesense pulse length*/
  data->fSigma=-1.0;    /*nonesense footprint width*/
  data->usable=1;
  data->zen=0.0;        /*straight down*/

  /*count number of wavebins*/
  data->nBins=0;
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))data->nBins++;

  data->wave=falloc(data->nBins,"waveform",0);
  data->z=dalloc(data->nBins,"z",0);
  if(gediIO->ground)data->ground=falloc(data->nBins,"ground",0);
  data->useID=0;
  data->demGround=0;

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  i=0;
  while(fgets(line,400,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(gediIO->useInt){  /*read intensity*/
        if(gediIO->ground==0){   /*don't read ground*/
          if(sscanf(line,"%s %s",temp1,temp2)==2){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp2);
            i++;
          }
        }else{                   /*read ground*/
          if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp2);
            data->ground[i]=atof(temp4);
            i++;
          }
        }/*ground switch*/
      }else if(gediIO->useFrac){ /*read fraction*/
        if(gediIO->ground==0){   /*don't read ground*/
         if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp3);
            i++;
          }
        }else{                   /*read ground*/
          if(sscanf(line,"%s %s %s %s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10)==10){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp8);
            data->ground[i]=atof(temp10);
            i++;
          }
        }/*ground switch*/
      }else{ /*read count*/
        if(gediIO->ground==0){   /*don't read ground*/
          if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp3);
            i++;
          }
        }else{                   /*read ground*/
          if(sscanf(line,"%s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7)==7){
            data->z[i]=atof(temp1);
            data->wave[i]=atof(temp5);
            data->ground[i]=atof(temp7);
            i++;
          }
        }/*ground switch*/
      }/*intensity or count switch*/
    }else{  /*read the header*/
      if(sscanf(line,"%s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7)==7){
        if(!strncasecmp(temp2,"fSigma",6)){
          data->fSigma=atof(temp3);
          data->pSigma=atof(temp5);
        }
      }
      if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
        if(!strncasecmp(temp2,"waveID",6)){
          data->useID=1;
          strcpy(&(data->waveID[0]),temp3);
        }else if(!strncasecmp(temp2,"meanScanAng",11)){
          data->zen=atof(temp3);
        }
      }
      if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
        if(!strncasecmp(temp2,"coord",5)){
          data->lon=atof(temp3);
          data->lat=atof(temp4);
        }else if(!strncasecmp(temp2,"ground",6)){
          if(!gediIO->dontTrustGround){
            data->gElev=atof(temp3);
            data->slope=atof(temp4);
            data->demGround=1;
          }else data->demGround=0;
        }
      }
      if(sscanf(line,"%s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6)==6){
        if(!strncasecmp(temp2,"density",7)){
          data->pointDense=atof(temp4);
          data->beamDense=atof(temp6);
        }else if(!strncasecmp(temp2,"lvis",4)){
          data->res=atof(temp4);
          data->zen=atof(temp6);
        }
      }
    }
  }/*line loop*/

  if(data->res<=0.0)data->res=gediIO->res;

  /*add up energy*/
  data->totE=0.0;
  for(i=0;i<data->nBins;i++)data->totE+=data->wave[i];


  gediIO->res=gediIO->den->res=gediIO->gFit->res=fabs(data->z[1]-data->z[0]);
  if(gediIO->den->res<TOL)data->usable=0;
  if(data->totE<=0.0)data->usable=0;
  if(gediIO->ground==0){   /*set to blank*/
    data->cov=-1.0;
    data->gLap=-1.0;
    data->gMinimum=-1.0;
    data->gInfl=-1.0;
  }
  if(data->fSigma<0.0)data->fSigma=gediIO->fSigma;
  if(data->pSigma<0.0)data->pSigma=gediIO->pSigma;

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  /*set up number of messages*/
  if((gediIO->nFiles>gediIO->nMessages)&&(gediIO->nMessages>0))gediIO->nMessages=(int)(gediIO->nFiles/gediIO->nMessages);
  else                                gediIO->nMessages=1;

  return(data);
}/*readASCIIdata*/


/*####################################################*/
/*turn data into HDF5 structure*/

gediHDF *arrangeGEDIhdf(dataStruct **data,gediIOstruct *gediIO)
{
  gediHDF *hdfData=NULL;
  void trimDataLength(dataStruct **,gediHDF *);

  /*allocate space for all*/
  if(!(hdfData=(gediHDF *)calloc(1,sizeof(gediHDF)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  hdfData->nWaves=gediIO->nFiles;
  hdfData->pSigma=data[0]->pSigma;
  hdfData->fSigma=data[0]->fSigma;
  /*beams*/
  hdfData->z0=falloc(hdfData->nWaves,"top elevations",0);       /*wave elevations*/
  hdfData->zN=falloc(hdfData->nWaves,"bottom elevations",0);       /*wave elevations*/
  hdfData->lon=dalloc(hdfData->nWaves,"lon",0);     /*longitudes*/
  hdfData->lat=dalloc(hdfData->nWaves,"lat",0);    /*latitudes*/
  hdfData->slope=falloc(hdfData->nWaves,"slope",0);    /*ground slope*/
  hdfData->cov=falloc(hdfData->nWaves,"cov",0);      /*canopy cover*/
  hdfData->gElev=falloc(hdfData->nWaves,"ground elevation, CofG",0);    /*ground elevation, CofG*/
  hdfData->demElev=falloc(hdfData->nWaves,"ground elevation, DEM",0);  /*ground elevation, DEM*/
  hdfData->beamDense=falloc(hdfData->nWaves,"beamDense",0);/*beam density*/
  hdfData->pointDense=falloc(hdfData->nWaves,"pointDense",0);/*point density*/
  hdfData->zen=falloc(hdfData->nWaves,"zen",0);      /*scan angles, or mean angles*/

  /*trim and copy data*/
  trimDataLength(data,hdfData);

  return(hdfData);
}/*arrangeGEDIhdf*/


/*####################################################*/
/*trim all arrays to be the same length*/

void trimDataLength(dataStruct **data,gediHDF *hdfData)
{
  int i=0,j=0,maxBins=0,maxID=0;
  uint64_t place=0;
  int *start=NULL,*end=NULL;
  float tot=0,thresh=0,res=0;
  float buffer=0,cumul=0;

  /*buffer length in metres*/
  buffer=30.0;

  /*allocate usable bins*/
  start=ialloc(hdfData->nWaves,"wave starts",0);
  end=ialloc(hdfData->nWaves,"end starts",0);


  /*determine start and end of all waves*/
  maxBins=maxID=-1;
  for(i=0;i<hdfData->nWaves;i++){
    /*total energy for a threshod*/
    tot=0.0;
    for(j=0;j<data[i]->nBins;j++)tot+=data[i]->wave[j];
    thresh=0.005*tot;

    /*determine used bounds*/
    cumul=0.0;
    for(j=0;j<data[i]->nBins;j++){
      cumul+=data[i]->wave[j];
      if(cumul>=thresh){
        start[i]=j;
        break;
      }
    }
    cumul=0.0;
    for(j=data[i]->nBins-1;j>=0;j--){
      cumul+=data[i]->wave[j];
      if(cumul>=thresh){
        end[i]=j;
        break;
      }
    }

    /*add a buffer for later smoothing*/
    res=fabs(data[i]->z[0]-data[i]->z[data[i]->nBins-1])/(float)(data[i]->nBins-1);
    start[i]-=(int)(buffer/res);
    end[i]+=(int)(buffer/res);
    if(start[i]<0)start[i]=0;
    if(end[i]>=data[i]->nBins)end[i]=data[i]->nBins-1;

    /*determine max*/
    if((end[i]-start[i])>maxBins)maxBins=end[i]-start[i];
    if(data[i]->useID){
      if(((int)strlen(data[i]->waveID)+1)>maxID)maxID=(int)strlen(data[i]->waveID)+1;
    }
  }
  hdfData->nBins=maxBins;
  if(maxID>0)hdfData->idLength=maxID;
  else       hdfData->idLength=7;

  /*allocate funny arrays needed by HDF5 library*/
  hdfData->wave=falloc(hdfData->nWaves*hdfData->nBins,"waveforms",0);
  hdfData->ground=falloc(hdfData->nWaves*hdfData->nBins,"ground waves",0);
  hdfData->waveID=challoc(hdfData->nWaves*hdfData->idLength,"wave IDs",0);

  /*copy arrays*/
  for(i=0;i<hdfData->nWaves;i++){
    /*range and resolution*/
    res=fabs(data[i]->z[0]-data[i]->z[data[i]->nBins-1])/(float)(data[i]->nBins-1);
    hdfData->z0[i]=data[i]->z[start[i]];
    hdfData->zN[i]=hdfData->z0[i]-(float)hdfData->nBins*res;
    /*copy data*/
    for(j=start[i];j<end[i];j++){
      place=(uint64_t)i*(uint64_t)hdfData->nBins+(uint64_t)j-(uint64_t)start[i];
      hdfData->wave[place]=data[i]->wave[j];
      hdfData->ground[place]=data[i]->ground[j];
    }
    /*pad end if not long enough*/
    for(j=end[i];j<(maxBins+start[i]);j++){
      place=(uint64_t)i*(uint64_t)hdfData->nBins+(uint64_t)j-(uint64_t)start[i];
      hdfData->wave[place]=0.0;
      hdfData->ground[place]=0.0;
    }

    hdfData->lon[i]=data[i]->lon;
    hdfData->lat[i]=data[i]->lat;
    hdfData->slope[i]=data[i]->slope;
    hdfData->cov[i]=data[i]->cov;
    hdfData->gElev[i]=data[i]->gElev;
    hdfData->demElev[i]=data[i]->demGround;
    if(maxID>0)strcpy(&(hdfData->waveID[i*hdfData->idLength]),data[i]->waveID);
    else       sprintf(&(hdfData->waveID[i*hdfData->idLength]),"%d",i);
    hdfData->beamDense[i]=data[i]->beamDense;
    hdfData->pointDense[i]=data[i]->pointDense;
    hdfData->zen[i]=data[i]->zen;
  }/*wave loop*/

  TIDY(start);
  TIDY(end);
  return;
}/*trimDataLength*/

/*the end*/
/*####################################################*/

