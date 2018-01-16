#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "hdf5.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "libLasRead.h"
#include "libLidVoxel.h"
#include "libOctree.h"
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
      TTIDY((void **)data[i]->wave,data[i]->nWaveTypes);
      TTIDY((void **)data[i]->ground,data[i]->nWaveTypes);
      TIDY(data[i]->noised);
      TIDY(data[i]->totE);
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
  int i=0,ind=0,*numb=NULL;
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
  data->nWaveTypes=(int)(gediIO->useInt+gediIO->useCount+gediIO->useFrac);
  if(data->nWaveTypes==1)data->useType=0;

  /*count number of wavebins*/
  data->nBins=0;
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))data->nBins++;

  /*is there usable data?*/
  if(data->nBins==0){
    data->usable=0;
  }else{
    data->wave=fFalloc(data->nWaveTypes,"waves",0);
    if(gediIO->ground)data->ground=fFalloc(data->nWaveTypes,"ground",0);
    data->z=dalloc(data->nBins,"z",0);
    for(i=0;i<data->nWaveTypes;i++){
    data->wave[i]=falloc(data->nBins,"waveform",0);
      if(gediIO->ground)data->ground[i]=falloc(data->nBins,"ground",0);
    }
    data->useID=0;
    data->demGround=0;

    /*rewind to start of file*/
    if(fseek(ipoo,(long)0,SEEK_SET)){
      fprintf(stderr,"fseek error\n");
      exit(1);
    }

    /*read data*/
    numb=ialloc(data->nWaveTypes,"number",0);
    for(i=0;i<data->nWaveTypes;i++)numb[i]=0;
    while(fgets(line,400,ipoo)!=NULL){
      if(strncasecmp(line,"#",1)){
        if(gediIO->useInt){  /*read intensity*/
          ind=0;
          if(gediIO->ground==0){   /*don't read ground*/
            if(sscanf(line,"%s %s",temp1,temp2)==2){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp2);
              numb[ind]++;
            }
          }else{                   /*read ground*/
            if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp2);
              data->ground[ind][numb[ind]]=atof(temp4);
              numb[ind]++;
            }
          }/*ground switch*/
        }
        if(gediIO->useFrac){ /*read fraction*/
          ind=(int)(gediIO->useInt+gediIO->useCount);
          if(gediIO->ground==0){   /*don't read ground*/
           if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp3);
              numb[ind]++;
            }
          }else{                   /*read ground*/
            if(sscanf(line,"%s %s %s %s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10)==10){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp8);
              data->ground[ind][numb[ind]]=atof(temp10);
              numb[ind]++;
            }
          }/*ground switch*/
        }
        if(gediIO->useCount){ /*read count*/
          ind=(int)gediIO->useInt;
          if(gediIO->ground==0){   /*don't read ground*/
            if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp3);
              numb[ind]++;
            }
          }else{                   /*read ground*/
            if(sscanf(line,"%s %s %s %s %s %s %s",temp1,temp2,temp3,temp4,temp5,temp6,temp7)==7){
              data->z[numb[ind]]=atof(temp1);
              data->wave[ind][numb[ind]]=atof(temp5);
              data->ground[ind][numb[ind]]=atof(temp7);
              numb[ind]++;
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
    TIDY(numb);

    if(data->res<=0.0)data->res=gediIO->res;

    /*add up energy*/
    data->totE=falloc(data->nWaveTypes,"",0);
    for(ind=0;ind<data->nWaveTypes;ind++){
      data->totE[ind]=0.0;
      for(i=0;i<data->nBins;i++)data->totE[ind]+=data->wave[ind][i];
    }

    gediIO->res=gediIO->den->res=gediIO->gFit->res=fabs(data->z[1]-data->z[0]);
    if(gediIO->den->res<TOL)data->usable=0;
    if(data->totE[data->useType]<=0.0)data->usable=0;
    if(gediIO->ground==0){   /*set to blank*/
      data->cov=-1.0;
      data->gLap=-1.0;
      data->gMinimum=-1.0;
      data->gInfl=-1.0;
    }
    if(data->fSigma<0.0)data->fSigma=gediIO->fSigma;
    if(data->pSigma<0.0)data->pSigma=gediIO->pSigma;
  }/*check there is dara*/

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
  int i=0;
  gediHDF *hdfData=NULL;
  void trimDataLength(dataStruct **,gediHDF *);

  /*allocate space for all*/
  if(!(hdfData=(gediHDF *)calloc(1,sizeof(gediHDF)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*count nuber of usable waves*/
  hdfData->nWaves=0;
  for(i=0;i<gediIO->nFiles;i++)if(data[i]->usable)hdfData->nWaves++;

  /*per scan*/
  hdfData->pSigma=data[0]->pSigma;
  hdfData->fSigma=data[0]->fSigma;
  hdfData->nPbins=0;
  hdfData->pulse=NULL;
  /*per beam*/
  hdfData->z0=falloc(hdfData->nWaves,"top elevations",0);       /*wave elevations*/
  hdfData->zN=falloc(hdfData->nWaves,"bottom elevations",0);       /*wave elevations*/
  hdfData->lon=dalloc(hdfData->nWaves,"lon",0);     /*longitudes*/
  hdfData->lat=dalloc(hdfData->nWaves,"lat",0);    /*latitudes*/
  hdfData->slope=falloc(hdfData->nWaves,"slope",0);    /*ground slope*/
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
  int ind=0,numb=0;;
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
    if(!data[i]->usable)continue;
    for(ind=0;ind<data[i]->nWaveTypes;ind++){
      /*total energy for a threshod*/
      tot=0.0;
      for(j=0;j<data[i]->nBins;j++)tot+=data[i]->wave[ind][j];
      thresh=0.005*tot;

      /*determine used bounds*/
      cumul=0.0;
      for(j=0;j<data[i]->nBins;j++){
        cumul+=data[i]->wave[ind][j];
        if(cumul>=thresh){
          start[i]=j;
          break;
        }
      }
      cumul=0.0;
      for(j=data[i]->nBins-1;j>=0;j--){
        cumul+=data[i]->wave[ind][j];
        if(cumul>=thresh){
          end[i]=j;
          break;
        }
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
  hdfData->nTypeWaves=data[0]->nWaveTypes;
  hdfData->wave=fFalloc(hdfData->nTypeWaves,"waveforms",0);
  hdfData->ground=fFalloc(hdfData->nTypeWaves,"ground waves",0);
  for(ind=0;ind<hdfData->nTypeWaves;ind++){
    hdfData->wave[ind]=falloc(hdfData->nWaves*hdfData->nBins,"waveforms",0);
    hdfData->ground[ind]=falloc(hdfData->nWaves*hdfData->nBins,"ground waves",0);
  }
  hdfData->waveID=challoc(hdfData->nWaves*hdfData->idLength,"wave IDs",0);

  /*copy arrays*/
  numb=0;
  for(i=0;i<hdfData->nWaves;i++){
    if(!data[i]->usable)continue;
    /*range and resolution*/
    res=fabs(data[i]->z[0]-data[i]->z[data[i]->nBins-1])/(float)(data[i]->nBins-1);
    hdfData->z0[numb]=data[i]->z[start[i]];
    hdfData->zN[numb]=hdfData->z0[numb]-(float)hdfData->nBins*res;
    /*copy data*/
    for(ind=0;ind<hdfData->nTypeWaves;ind++){
      for(j=start[i];j<end[i];j++){
        place=(uint64_t)numb*(uint64_t)hdfData->nBins+(uint64_t)j-(uint64_t)start[i];
        hdfData->wave[ind][place]=data[i]->wave[ind][j];
        hdfData->ground[ind][place]=data[i]->ground[ind][j];
      }
      /*pad end if not long enough*/
      for(j=end[i];j<(maxBins+start[i]);j++){
        place=(uint64_t)numb*(uint64_t)hdfData->nBins+(uint64_t)j-(uint64_t)start[i];
        hdfData->wave[ind][place]=0.0;
        hdfData->ground[ind][place]=0.0;
      }
    }

    hdfData->lon[numb]=data[i]->lon;
    hdfData->lat[numb]=data[i]->lat;
    hdfData->slope[numb]=data[i]->slope;
    hdfData->gElev[numb]=data[i]->gElev;
    hdfData->demElev[numb]=data[i]->demGround;
    if(maxID>0)strcpy(&(hdfData->waveID[numb*hdfData->idLength]),data[i]->waveID);
    else       sprintf(&(hdfData->waveID[numb*hdfData->idLength]),"%d",i);
    hdfData->beamDense[numb]=data[i]->beamDense;
    hdfData->pointDense[numb]=data[i]->pointDense;
    hdfData->zen[numb]=data[i]->zen;
    numb++;
  }/*wave loop*/

  TIDY(start);
  TIDY(end);
  return;
}/*trimDataLength*/


/*####################################################*/
/*write data to HDF5*/

void writeGEDIhdf(gediHDF *hdfData,char *namen)
{
  hid_t file;         /* Handles */

  /*open new file*/
  file=H5Fcreate(namen,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

 /*write header*/
  write1dIntHDF5(file,"NWAVES",&hdfData->nWaves,1);
  write1dIntHDF5(file,"NBINS",&hdfData->nBins,1);
  write1dIntHDF5(file,"NTYPEWAVES",&hdfData->nTypeWaves,1);
  write1dIntHDF5(file,"IDLENGTH",&hdfData->idLength,1);
  write1dFloatHDF5(file,"PSIGMA",&hdfData->pSigma,1);
  write1dFloatHDF5(file,"FSIGMA",&hdfData->fSigma,1);
  write1dIntHDF5(file,"NPBINS",&hdfData->nPbins,1);
  /*write datasets*/
  write1dDoubleHDF5(file,"LON0",hdfData->lon,hdfData->nWaves);
  write1dDoubleHDF5(file,"LAT0",hdfData->lat,hdfData->nWaves);
  if(hdfData->ground){
    write1dFloatHDF5(file,"SLOPE",hdfData->slope,hdfData->nWaves);
    write1dFloatHDF5(file,"ZG",hdfData->gElev,hdfData->nWaves);
    if(hdfData->demElev)write1dFloatHDF5(file,"ZGDEM",hdfData->demElev,hdfData->nWaves);
  }
  write1dFloatHDF5(file,"BEAMDENSE",hdfData->beamDense,hdfData->nWaves);
  write1dFloatHDF5(file,"POINTDENSE",hdfData->pointDense,hdfData->nWaves);
  write1dFloatHDF5(file,"INCIDENTANGLE",hdfData->zen,hdfData->nWaves);
  if(hdfData->nTypeWaves==3){
    write2dFloatHDF5(file,"RXWAVEINT",hdfData->wave[0],hdfData->nWaves,hdfData->nBins);
    if(hdfData->ground)write2dFloatHDF5(file,"GRWAVEINT",hdfData->ground[0],hdfData->nWaves,hdfData->nBins);
    write2dFloatHDF5(file,"RXWAVECOUNT",hdfData->wave[1],hdfData->nWaves,hdfData->nBins);
    if(hdfData->ground)write2dFloatHDF5(file,"GRWAVECOUNT",hdfData->ground[1],hdfData->nWaves,hdfData->nBins);
    write2dFloatHDF5(file,"RXWAVEFRAC",hdfData->wave[2],hdfData->nWaves,hdfData->nBins);
    if(hdfData->ground)write2dFloatHDF5(file,"GRWAVEFRAC",hdfData->ground[2],hdfData->nWaves,hdfData->nBins);
  }else if(hdfData->nTypeWaves==1){
    fprintf(stderr,"We are not set up for wrinting only one HDF5 wave yet");
    exit(1);
    write2dFloatHDF5(file,"RXWAVE",hdfData->wave[0],hdfData->nWaves,hdfData->nBins);
    if(hdfData->ground)write2dFloatHDF5(file,"GRWAVE",hdfData->ground[0],hdfData->nWaves,hdfData->nBins);
  }else{
    fprintf(stderr,"Can't handle this number of waveforms");
    exit(1);
  }
  write1dFloatHDF5(file,"Z0",hdfData->z0,hdfData->nWaves);
  write1dFloatHDF5(file,"ZN",hdfData->zN,hdfData->nWaves);
  write2dCharHDF5(file,"WAVEID",hdfData->waveID,hdfData->nWaves,hdfData->idLength);

  if(hdfData->nPbins>0){
    write1dFloatHDF5(file,"PRES",&hdfData->pRes,1);
    write1dFloatHDF5(file,"PULSE",hdfData->pulse,hdfData->nPbins);
  }

  /*close file*/
  if(H5Fclose(file)){
    fprintf(stderr,"Issue closing file\n");
    exit(1);
  }
  fprintf(stdout,"Waveforms written to %s\n",namen);
  return;
}/*writeGEDIhdf*/


/*####################################################*/
/*read HDF5 GEDI data into structure*/

gediHDF *readGediHDF(char *namen,gediIOstruct *gediIO)
{
  int nWaves=0,nBins=0;
  int *tempI=NULL,ind=0;
  float *tempF=NULL;
  gediHDF *hdfData=NULL;
  hid_t file;         /* Handles */
  void checkNwavesDF(int,int);

  /*allocate space for all*/
  if(!(hdfData=(gediHDF *)calloc(1,sizeof(gediHDF)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*open file*/
  fprintf(stdout,"Reading %s\n",namen);
  file=H5Fopen(namen,H5F_ACC_RDONLY,H5P_DEFAULT);

  /*read the header*/
  tempF=read1dFloatHDF5(file,"PSIGMA",&nWaves);
  hdfData->pSigma=*tempF;
  TIDY(tempF);
  tempF=read1dFloatHDF5(file,"FSIGMA",&nWaves);
  hdfData->fSigma=*tempF;
  TIDY(tempF);
  tempI=read1dIntHDF5(file,"NWAVES",&nWaves);
  hdfData->nWaves=*tempI;
  TIDY(tempI);
  tempI=read1dIntHDF5(file,"NBINS",&nWaves);
  hdfData->nBins=*tempI;
  TIDY(tempI);
  tempI=read1dIntHDF5(file,"NPBINS",&nWaves);
  hdfData->nPbins=*tempI;
  TIDY(tempI);
  tempI=read1dIntHDF5(file,"NTYPEWAVES",&nWaves);
  hdfData->nTypeWaves=*tempI;
  TIDY(tempI);
  tempI=read1dIntHDF5(file,"IDLENGTH",&nWaves);
  hdfData->idLength=*tempI;
  TIDY(tempI);

  /*read ancillary data*/
  hdfData->z0=read1dFloatHDF5(file,"Z0",&nWaves);
  checkNwavesDF(nWaves,hdfData->nWaves);
  hdfData->zN=read1dFloatHDF5(file,"ZN",&nWaves);
  checkNwavesDF(nWaves,hdfData->nWaves);
  if(gediIO->ground){
    hdfData->slope=read1dFloatHDF5(file,"SLOPE",&nWaves);
    checkNwavesDF(nWaves,hdfData->nWaves);
    hdfData->gElev=read1dFloatHDF5(file,"ZG",&nWaves);
    checkNwavesDF(nWaves,hdfData->nWaves);
  }
  hdfData->beamDense=read1dFloatHDF5(file,"BEAMDENSE",&nWaves);
  checkNwavesDF(nWaves,hdfData->nWaves);
  hdfData->pointDense=read1dFloatHDF5(file,"POINTDENSE",&nWaves);
  checkNwavesDF(nWaves,hdfData->nWaves);
  hdfData->zen=read1dFloatHDF5(file,"INCIDENTANGLE",&nWaves);
  checkNwavesDF(nWaves,hdfData->nWaves);
  hdfData->lon=read1dDoubleHDF5(file,"LON0",&nWaves);
  checkNwavesDF(nWaves,hdfData->nWaves);
  hdfData->lat=read1dDoubleHDF5(file,"LAT0",&nWaves);
  checkNwavesDF(nWaves,hdfData->nWaves);
  hdfData->waveID=read15dCharHDF5(file,"WAVEID",&nWaves,&nBins);
  checkNwavesDF(nWaves,hdfData->nWaves);
  if(hdfData->nPbins>0){
    hdfData->pulse=read1dFloatHDF5(file,"PULSE",&nBins);
    checkNwavesDF(nBins,hdfData->nPbins);
    tempF=read1dFloatHDF5(file,"PRES",&nWaves);
    hdfData->pRes=*tempF;
    TIDY(tempF);
  }else hdfData->pulse=NULL;

  //float *demElev;  /*ground elevation, DEM*/


  /*determine how many waveforms we want to read*/
  if((int)(gediIO->useInt+gediIO->useCount+gediIO->useFrac)>hdfData->nTypeWaves){
    fprintf(stderr,"Not enough waveform types for that option set. %d %d",hdfData->nTypeWaves,gediIO->useInt+gediIO->useCount+gediIO->useFrac);
    exit(1);
  }else hdfData->nTypeWaves=(int)(gediIO->useInt+gediIO->useCount+gediIO->useFrac);

  /*allocate waveform space*/
  hdfData->wave=fFalloc(hdfData->nTypeWaves,"hdf waveforms",0);
  if(gediIO->ground)hdfData->ground=fFalloc(hdfData->nTypeWaves,"hdf waveforms",0);

  /*read data*/
  if(gediIO->useInt){
    ind=0;
    hdfData->wave[ind]=read15dFloatHDF5(file,"RXWAVEINT",&nWaves,&nBins);
    checkNwavesDF(nWaves,hdfData->nWaves);
    checkNwavesDF(nBins,hdfData->nBins);
    if(gediIO->ground){
      hdfData->ground[ind]=read15dFloatHDF5(file,"GRWAVEINT",&nWaves,&nBins);
      checkNwavesDF(nWaves,hdfData->nWaves);
      checkNwavesDF(nBins,hdfData->nBins);
    }
  }
  if(gediIO->useCount){
    ind=(int)gediIO->useInt;
    hdfData->wave[ind]=read15dFloatHDF5(file,"RXWAVECOUNT",&nWaves,&nBins);
    checkNwavesDF(nWaves,hdfData->nWaves);
    checkNwavesDF(nBins,hdfData->nBins);
    if(gediIO->ground){
      hdfData->ground[ind]=read15dFloatHDF5(file,"GRWAVECOUNT",&nWaves,&nBins);
      checkNwavesDF(nWaves,hdfData->nWaves);
      checkNwavesDF(nBins,hdfData->nBins);
    }
  }
  if(gediIO->useFrac){
    ind=(int)(gediIO->useInt+gediIO->useCount);
    hdfData->wave[ind]=read15dFloatHDF5(file,"RXWAVEFRAC",&nWaves,&nBins);
    checkNwavesDF(nWaves,hdfData->nWaves);
    checkNwavesDF(nBins,hdfData->nBins);
    if(gediIO->ground){
      hdfData->ground[ind]=read15dFloatHDF5(file,"GRWAVEFRAC",&nWaves,&nBins);
      checkNwavesDF(nWaves,hdfData->nWaves);
      checkNwavesDF(nBins,hdfData->nBins);
    }
  }

  /*close file*/
  if(H5Fclose(file)){
    fprintf(stderr,"Issue closing file\n");
    exit(1);
  }
  return(hdfData);
}/*readGediHDF*/


/*####################################################*/
/*check that number of waves match*/

void checkNwavesDF(int nRead,int nWaves)
{
  if(nRead!=nWaves){
    fprintf(stderr,"number of waves mismatch: read %d, expecting %d\n",nRead,nWaves);
    exit(1);
  }

  return;
}/*checkNwavesDF*/


/*####################################################*/
/*tidy GEDI HDF data structire*/

gediHDF *tidyGediHDF(gediHDF *hdfData)
{

  if(hdfData){
    TTIDY((void **)hdfData->wave,hdfData->nTypeWaves);
    TTIDY((void **)hdfData->ground,hdfData->nTypeWaves);
    TIDY(hdfData->waveID);
    TIDY(hdfData->pulse);
    TIDY(hdfData->z0);       /*wave top elevations*/
    TIDY(hdfData->zN);       /*wave bottom elevations*/
    TIDY(hdfData->lon);     /*longitudes*/
    TIDY(hdfData->lat);     /*latitudes*/
    TIDY(hdfData->slope);    /*ground slope*/
    TIDY(hdfData->gElev);    /*ground elevation, CofG*/
    TIDY(hdfData->demElev);  /*ground elevation, DEM*/
    TIDY(hdfData->beamDense);/*beam density*/
    TIDY(hdfData->pointDense);/*point density*/
    TIDY(hdfData->zen);      /*scan angles, or mean angles*/
    TIDY(hdfData);
  }

  return(hdfData);
}/*tidyHDFdata*/


/*####################################################*/
/*read LVIS HDF file*/

dataStruct *unpackHDFlvis(char *namen,lvisHDF **hdfLvis,gediIOstruct *gediIO,int numb)
{
  int i=0,botBin=0;
  int findLvisBottom(float *wave,int nBins);
  dataStruct *data=NULL;
  float *tempPulse=NULL;
  float pulseLenFromTX(float *,int);
  double dx=0,dy=0,scale=0;  /*for padded LVIS files*/

  /*read data if needed*/
  if(*hdfLvis==NULL){
    *hdfLvis=readLVIShdf(namen);
    gediIO->nFiles=hdfLvis[0]->nWaves;
    gediIO->ground=0;
  }

  /*allocate space*/
  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  data->useID=1;
  data->nBins=hdfLvis[0]->nBins;
  data->nWaveTypes=1;
  data->useType=0;
  data->wave=fFalloc(data->nWaveTypes,"waveform",0);
  data->wave[0]=falloc(data->nBins,"waveform",0);
  data->totE=falloc(data->nWaveTypes,"totE",0);
  data->z=dalloc(data->nBins,"z",0);
  data->ground=NULL;
  data->demGround=0;
  data->pSigma=-1.0;    /*nonesense pulse length*/
  data->fSigma=-1.0;    /*nonesense footprint width*/
  data->usable=1;

  /*copy data to structure*/
  data->zen=hdfLvis[0]->zen[numb];
  data->res=fabs(hdfLvis[0]->z0[numb]-hdfLvis[0]->z1023[numb])/(float)hdfLvis[0]->nBins;
  if(gediIO->den)gediIO->den->res=data->res;
  if(gediIO->gFit)gediIO->gFit->res=data->res;
  data->totE[data->useType]=0.0;
  for(i=0;i<hdfLvis[0]->nBins;i++){
    data->wave[data->useType][i]=(float)hdfLvis[0]->wave[numb][i];
    data->z[i]=(double)(hdfLvis[0]->z0[numb]-(float)i*data->res);
    data->totE[data->useType]+=data->wave[data->useType][i];
  }
  if(gediIO->den->res<TOL)data->usable=0;
  if(data->totE[data->useType]<=0.0)data->usable=0;
  if(gediIO->ground==0){   /*set to blank*/
    data->cov=-1.0;
    data->gLap=-1.0;
    data->gMinimum=-1.0;
    data->gInfl=-1.0;
  }

  /*for setting coordinate, use range of lowest ground return*/
  botBin=findLvisBottom(data->wave[data->useType],data->nBins);

  dx=hdfLvis[0]->lon1023[numb]-hdfLvis[0]->lon0[numb];
  dy=hdfLvis[0]->lat1023[numb]-hdfLvis[0]->lat0[numb];
  scale=(double)botBin/1024.0;
  data->lon=hdfLvis[0]->lon0[numb]+scale*dx;
  data->lat=hdfLvis[0]->lat0[numb]+scale*dy;
  data->lfid=hdfLvis[0]->lfid[numb];
  data->shotN=hdfLvis[0]->shotN[numb];
  sprintf(data->waveID,"%d.%d",hdfLvis[0]->lfid[numb],hdfLvis[0]->shotN[numb]);

  /*analyse pulse*/
  if(gediIO->readPsigma){
    tempPulse=falloc(hdfLvis[0]->pBins,"temp pulse",0);
    for(i=0;i<hdfLvis[0]->pBins;i++)tempPulse[i]=(float)hdfLvis[0]->pulse[numb][i];
    data->pSigma=pulseLenFromTX(tempPulse,hdfLvis[0]->pBins);
    TIDY(tempPulse);
  }else data->pSigma=gediIO->pSigma;

  if(data->fSigma<0.0)data->fSigma=gediIO->fSigma;
  if(data->pSigma<0.0)data->pSigma=gediIO->pSigma;

  /*set up number of messages*/
  if((gediIO->nMessages>0)&&(hdfLvis[0]->nWaves>gediIO->nMessages))gediIO->nMessages=(int)(hdfLvis[0]->nWaves/gediIO->nMessages);
  else                                 gediIO->nMessages=1;

  return(data);
}/*unpackHDFlvis*/


/*####################################################*/
/*find bottom bin for determing coordinate*/

int findLvisBottom(float *wave,int nBins)
{
  int i=0;
  float *tempWave=NULL;
  denPar den;
  char found=0;
  void setDenoiseDefault(denPar *);

  /*set denoising parameters*/
  setDenoiseDefault(&den);
  den.varNoise=1;
  den.statsLen=10.0;
  den.noiseTrack=0;
  den.threshScale=4.0;

  /*find point to calcuklate coordinate from*/
  tempWave=processFloWave(wave,nBins,&den,1.0);
  found=0;
  for(i=nBins-1;i>=0;i--){
    if(tempWave[i]>TOL){
      found=1;
      break;
    }
  }
  TIDY(tempWave);
  if(found==0)i=nBins/2;

  return(i);
}/*findLvisBottom*/


/*####################################################*/
/*determine pulse width from TXwave*/

float pulseLenFromTX(float *pulse,int nBins)
{
  int i=0;
  float pSigma=0;
  float *denoised=NULL;
  float tot=0,CofG=0;
  denPar den;
  void setDenoiseDefault(denPar *);

  /*denoise*/
  setDenoiseDefault(&den);
  den.varNoise=1;
  den.threshScale=5.0;
  den.noiseTrack=1;
  den.minWidth=5;
  den.statsLen=3.0;
  den.res=0.15;
  denoised=processFloWave(pulse,nBins,&den,1.0);

  /*CofG*/
  CofG=tot=0.0;
  for(i=0;i<nBins;i++){
    CofG+=(float)i*den.res*denoised[i];
    tot+=denoised[i];
  }
  if(tot>0.0){
    CofG/=tot;

    pSigma=0.0;
    for(i=0;i<nBins;i++)pSigma+=((float)i*den.res*-CofG)*((float)i*den.res-CofG)*denoised[i];
    pSigma=sqrt(pSigma/tot);

  }else pSigma=-1.0;

  TIDY(denoised);
  return(pSigma);
}/*pulseLenFromTX*/


/*####################################################*/
/*read LVIS binary data*/

dataStruct *readBinaryLVIS(char *namen,lvisLGWstruct *lvis,int numb,gediIOstruct *gediIO)
{
  int i=0,botBin=0;
  int findLvisBottom(float *,int);
  dataStruct *data=NULL;
  float *tempPulse=NULL;
  float pulseLenFromTX(float *,int);
  double dx=0,dy=0,scale=0;

  /*do we need to read all the data*/
  if(lvis->data==NULL){
    /*read data*/
    lvis->data=readLVISlgw(namen,lvis);
    gediIO->nFiles=lvis->nWaves;
  }


  /*allocate space*/
  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  data->useID=1;
  data->nBins=lvis->nBins;
  data->useType=0;
  data->nWaveTypes=1;
  data->wave=fFalloc(data->nWaveTypes,"waveform",0);
  data->wave[data->useType]=falloc(data->nBins,"waveform",0);
  data->z=dalloc(data->nBins,"z",0);
  data->ground=NULL;
  data->demGround=0;
  data->pSigma=-1.0;    /*nonesense pulse length*/
  data->fSigma=-1.0;    /*nonesense footprint width*/
  data->usable=1;
  data->totE=falloc(data->nWaveTypes,"totE",0);

  /*copy data to structure*/
  data->zen=lvis->data[numb].zen;
  data->res=(lvis->data[numb].z0-lvis->data[numb].z431)/(float)lvis->nBins;
  if(gediIO->den)gediIO->den->res=data->res;
  if(gediIO->gFit)gediIO->gFit->res=data->res;
  data->totE[data->useType]=0.0;
  for(i=0;i<lvis->nBins;i++){
    data->wave[data->useType][i]=(float)lvis->data[numb].rxwave[i];
    data->z[i]=(double)(lvis->data[numb].z0-(float)i*data->res);
    data->totE[data->useType]+=data->wave[data->useType][i];
  }
  if(gediIO->den->res<TOL)data->usable=0;
  if(data->totE[data->useType]<=0.0)data->usable=0;
  if(gediIO->ground==0){   /*set to blank*/
    data->cov=-1.0;
    data->gLap=-1.0;
    data->gMinimum=-1.0;
    data->gInfl=-1.0;
  }

  /*find point to calculate coordinate from*/
  botBin=findLvisBottom(data->wave[data->useType],data->nBins);
  dx=lvis->data[numb].lon431-lvis->data[numb].lon0;
  dy=lvis->data[numb].lat431-lvis->data[numb].lat0;
  scale=(double)botBin/432.0;
  data->lon=lvis->data[numb].lon0+dx*scale;
  data->lat=lvis->data[numb].lat0+dy*scale;
  sprintf(data->waveID,"%d.%d",lvis->data[numb].lfid,lvis->data[numb].shotN);

  /*analyse pulse*/
  if(gediIO->readPsigma){
    tempPulse=falloc(80,"temp pulse",0);
    for(i=0;i<80;i++)tempPulse[i]=(float)lvis->data[numb].rxwave[i];
    data->pSigma=pulseLenFromTX(tempPulse,80);
    TIDY(tempPulse);
  }else data->pSigma=gediIO->pSigma;

  if(data->fSigma<0.0)data->fSigma=gediIO->fSigma;
  if(data->pSigma<0.0)data->pSigma=gediIO->pSigma;


  /*set up number of messages*/
  if((gediIO->nMessages>0)&&(lvis->nWaves>gediIO->nMessages))gediIO->nMessages=(int)(lvis->nWaves/gediIO->nMessages);
  else                              gediIO->nMessages=1;

  return(data);
}/*readBinaryLVIS*/


/*####################################*/
/*read lasFile and save relevant data*/

pCloudStruct *readALSdata(lasFile *las,gediRatStruct *gediRat,int nFile)
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
  if(!gediRat->readALSonce)useFile=checkFileBounds(las,gediRat->globMinX,gediRat->globMaxX,gediRat->globMinY,gediRat->globMaxY);
  else                     useFile=checkMultiFiles(las,gediRat->gNx,gediRat->coords,gediRat->maxSep);

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
      if(!gediRat->readALSonce){
        if((x>=gediRat->globMinX)&&(x<=gediRat->globMaxX)&&(y>=gediRat->globMinY)&&(y<=gediRat->globMaxY)&&\
           (z>-10000.0)&&(z<10000.0)&&(fabs((float)las->scanAng)<=gediRat->maxScanAng))usePoint=1;
        else usePoint=0;
      }else usePoint=checkMultiPoints(x,y,z,gediRat->gNx,gediRat->coords,gediRat->maxSep);

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

        /*map to octree if needed*/
        if(gediRat->useOctree)fillOctree(x,y,z,nFile,pUsed,gediRat->octree);

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
      if(gediRat->useShadow){
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
      if(!gediRat->useShadow){
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
  if(gediRat->readWave&&hasWave){  /*only leave files open for waveform*/
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
/*set GEDI grid or batch*/

void setGediGrid(gediIOstruct *gediIO,gediRatStruct *gediRat)
{
  void readFeetList(gediRatStruct *);
  void setRatBounds(gediRatStruct *);
  void readWavefront(gediRatStruct *,gediIOstruct *);

  /*footprint width*/
  if(gediRat->defWfront==0){   /*regular footprint*/
    if(gediRat->topHat==0)gediRat->maxSep=determineGaussSep(gediIO->fSigma,gediRat->iThresh);
    else                  gediRat->maxSep=gediIO->fSigma;
  }else{    /*read assymetric footprint*/
    readWavefront(gediRat,gediIO);
  }/*footprint width setting*/

  if(gediRat->doGrid){  /*it is a grid*/
    /*number of footprints*/
    gediRat->gNx=(int)((gediRat->gMaxX-gediRat->gMinX)/gediRat->gRes+1);
    gediRat->gNy=(int)((gediRat->gMaxY-gediRat->gMinY)/gediRat->gRes+1);

    /*global bounds*/
    gediRat->globMinX=gediRat->gMinX-gediRat->maxSep;
    gediRat->globMaxX=gediRat->gMaxX+gediRat->maxSep;
    gediRat->globMinY=gediRat->gMinY-gediRat->maxSep;
    gediRat->globMaxY=gediRat->gMaxY+gediRat->maxSep;
  }else if(gediRat->readALSonce){ /*it is a batch*/
    /*read list of coords*/
    if(gediRat->coords==NULL)readFeetList(gediRat);
    setRatBounds(gediRat);
  }else{   /*single footprint*/
    gediRat->gNx=gediRat->gNy=1;
    gediRat->globMinX=gediRat->coord[0]-gediRat->maxSep;
    gediRat->globMaxX=gediRat->coord[0]+gediRat->maxSep;
    gediRat->globMinY=gediRat->coord[1]-gediRat->maxSep;
    gediRat->globMaxY=gediRat->coord[1]+gediRat->maxSep;
  }

  if((gediIO->nMessages>1)&&(gediRat->gNx*gediRat->gNy)>gediIO->nMessages)gediIO->nMessages=(int)(gediRat->gNx*gediRat->gNy/gediIO->nMessages);
  else                                             gediIO->nMessages=1;

  /*allocate octree if needed*/
  if(gediRat->useOctree){
    gediRat->octree=allocateOctree(gediRat->octLevels,gediRat->nOctTop,\
            gediRat->globMinX,gediRat->globMaxX,gediRat->globMinY,gediRat->globMaxY);
  }else gediRat->octree=NULL;

  return;
}/*setGediGrid*/


/*###################################################*/

void readWavefront(gediRatStruct *gediRat,gediIOstruct *gediIO)
{
  int i=0,j=0,maxI=0;
  float len=0; /*total=0*/
  char line[20000];
  char *token=NULL;
  void setWavefrontRes(wFrontStruct *,float);
  FILE *ipoo=NULL;

  /*open file*/
  if((ipoo=fopen(gediRat->wavefront->frontFile,"r"))==NULL){
    fprintf(stderr,"Error opening wavefront file \"%s\"\n",gediRat->wavefront->frontFile);
    exit(1);
  }

  /*find file size*/
  j=0;
  maxI=0;
  while(fgets(line,20000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      i=0;
      token=strtok(line,",");
      while(token){
        token=strtok(NULL,",");
        i++;
      }
      if(i>maxI)maxI=i;
      j++;
    }
  }/*file size counting*/
  gediRat->wavefront->nX=maxI;
  gediRat->wavefront->nY=j;

  /*allocate space*/
  gediRat->wavefront->front=fFalloc(gediRat->wavefront->nX,"wavefront",0);
  for(i=0;i<gediRat->wavefront->nX;i++){
    gediRat->wavefront->front[i]=falloc(gediRat->wavefront->nY,"wavefront",j);
    for(j=0;j<gediRat->wavefront->nY;j++)gediRat->wavefront->front[i][j]=-1.0;  /*mark as blank*/
  }/*allocation*/

  /*rewind*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  /*total=0.0;*/
  j=gediRat->wavefront->nY-1;
  while(fgets(line,20000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      i=0;
      token=strtok(line,",");
      while(token){
        if(strncasecmp(token,",",1)){  /*is there a result*/
          gediRat->wavefront->front[i][j]=atof(token);
          /*total+=gediRat->wavefront->front[i][j];*/
        }
        token=strtok(NULL,",");
        i++;
      }
      j--;
    }
  }/*data reading*/

  /*tidy up*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  /*normalise*/
  /*for(i=0;i<gediRat->wavefront->nX;i++){
    for(j=0;j<gediRat->wavefront->nY;j++){
      if(gediRat->wavefront->front[i][j]>0.0)gediRat->wavefront->front[i][j]/=total;
    }
  }*//*normalisation*/

  /*determine resolution from footprint width*/
  setWavefrontRes(gediRat->wavefront,gediIO->fSigma);

  /*set bounds for reading*/
  gediRat->maxSep=-100000.0;
  len=(gediRat->wavefront->x0>(gediRat->wavefront->nX/2))?(float)gediRat->wavefront->x0*gediRat->wavefront->res:\
                  (float)(gediRat->wavefront->nX-gediRat->wavefront->x0)*gediRat->wavefront->res;
  gediRat->maxSep=(len>gediRat->maxSep)?len:gediRat->maxSep;
  len=(gediRat->wavefront->y0>(gediRat->wavefront->nY/2))?(float)gediRat->wavefront->y0*gediRat->wavefront->res:\
                  (float)(gediRat->wavefront->nY-gediRat->wavefront->y0)*gediRat->wavefront->res;
  gediRat->maxSep=(len>gediRat->maxSep)?len:gediRat->maxSep;

  return;
}/*readWavefront*/


/*##########################################################*/
/*determine wavefront resolution and peak*/

void setWavefrontRes(wFrontStruct *wavefront,float fSigma)
{
  int i=0,j=0;
  int ii=0,jj=0;
  int window=0,contN=0;
  float max=0,tot=0,mean=0;
  float xStdev=0,yStdev=0;

  /*find centre*/
  window=2;
  max=-1.0;
  for(i=0;i<wavefront->nX;i++){
    for(j=0;j<wavefront->nY;j++){
      mean=0.0;
      contN=0;
      for(ii=i-window;ii<(i+window);ii++){
        if((ii<0)||(ii>=wavefront->nX))continue;
        for(jj=j-window;jj<(j+window);jj++){
          if((jj<0)||(jj>=wavefront->nY))continue;
          if(wavefront->front[i][j]>=0.0){
            mean+=wavefront->front[ii][jj];
            contN++;
          }
        }
      }
      if(contN>0){
        mean/=(float)contN;
        if(mean>max){
          max=mean;
          wavefront->x0=i;
          wavefront->y0=j;
        }
      }
    }
  }

  /*determine width in two axes*/
  tot=xStdev=0.0;
  for(i=0;i<wavefront->nX;i++){
    if(wavefront->front[i][wavefront->y0]>=0.0){
      xStdev+=pow((float)(i-wavefront->x0),2.0)*wavefront->front[i][wavefront->y0];
      tot+=wavefront->front[i][wavefront->y0];
    }
  }
  xStdev=sqrt(xStdev/tot);
  tot=yStdev=0.0;
  for(j=0;j<wavefront->nY;j++){
    if(wavefront->front[wavefront->x0][j]>=0.0){
      yStdev+=pow((float)(j-wavefront->y0),2.0)*wavefront->front[wavefront->x0][j];
      tot+=wavefront->front[wavefront->x0][j];
    }
  }
  yStdev=sqrt(yStdev/tot);
  wavefront->res=2.0*fSigma/(xStdev+yStdev);

  return;
}/*setWavefrontRes*/


/*####################################*/
/*set bounds from list of coords*/

void setRatBounds(gediRatStruct *gediRat)
{
  int i=0;
  double minX=0,maxX=0;
  double minY=0,maxY=0;

  minX=minY=1000000000.0;
  maxX=maxY=-1000000000.0;
  for(i=0;i<gediRat->gNx;i++){
    if(gediRat->coords[i][0]<minX)minX=gediRat->coords[i][0];
    if(gediRat->coords[i][1]<minY)minY=gediRat->coords[i][1];
    if(gediRat->coords[i][0]>maxX)maxX=gediRat->coords[i][0];
    if(gediRat->coords[i][1]>maxY)maxY=gediRat->coords[i][1];
  }

  gediRat->globMinX=minX-gediRat->maxSep;
  gediRat->globMaxX=maxX+gediRat->maxSep;
  gediRat->globMinY=minY-gediRat->maxSep;
  gediRat->globMaxY=maxY+gediRat->maxSep;

  return;
}/*setRatBounds*/


/*####################################*/
/*read list of coordinates*/

void readFeetList(gediRatStruct *gediRat)
{
  int i=0;
  char line[200],temp1[50];
  char temp2[50],temp3[100];;
  FILE *ipoo=NULL;

  /*open file*/
  if((ipoo=fopen(gediRat->coordList,"r"))==NULL){
    fprintf(stderr,"Error opening input file list \"%s\"\n",gediRat->coordList);
    exit(1);
  }


  /*count number of lines*/
  i=0;
  while(fgets(line,200,ipoo)!=NULL)if(strncasecmp(line,"#",1))i++;

  /*allocate space*/
  gediRat->gNx=i;
  gediRat->gNy=1;
  gediRat->coords=dDalloc(gediRat->gNx,"coord list",0);
  gediRat->waveIDlist=chChalloc(gediRat->gNx,"wave ID list",0);

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read coordinate list*/
  i=0;
  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      gediRat->coords[i]=dalloc(2,"coord list",i+1);
      if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){ /*read coord and waveID*/
        gediRat->coords[i][0]=atof(temp1);
        gediRat->coords[i][1]=atof(temp2);
        gediRat->waveIDlist[i]=challoc((int)strlen(temp3)+1,"wave ID list",i+1);
        strcpy(gediRat->waveIDlist[i],temp3);
      }else if(sscanf(line,"%s %s",temp1,temp2)==2){
        gediRat->coords[i][0]=atof(temp1);
        gediRat->coords[i][1]=atof(temp2);
        sprintf(temp3,"%f.%f",gediRat->coords[i][0],gediRat->coords[i][1]);
        gediRat->waveIDlist[i]=challoc((int)strlen(temp3)+1,"wave ID list",i+1);
        strcpy(gediRat->waveIDlist[i],temp3);
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

void setGediPulse(gediIOstruct *gediIO,gediRatStruct *gediRat)
{
  int i=0;
  float fwhm=0;   /*FWHM in metres*/
  float x=0,y=0;
  float max=0,tot=0;
  void readSimPulse(gediIOstruct *,gediRatStruct *);


  if(!(gediIO->pulse=(pulseStruct *)calloc(1,sizeof(pulseStruct)))){
    fprintf(stderr,"error pulseStruct allocation.\n");
    exit(1);
  }


  if(gediIO->readPulse==0){  /*Gaussian pulse*/
    /*pulse length*/
    /*calculate sigma from FWHM*/
    if(gediIO->pSigma<0.0){  /*GEDI unless specificed*/
      fwhm=gediIO->pFWHM*0.2998/2.0;  /*time for two way*/
      gediIO->pSigma=fwhm/2.35482;  /* =2*sqrt(2*ln2) */
    }

    if(gediIO->pSigma>0.0){  /*if we are using a pulse width*/
      /*determine number of bins*/
      gediIO->pulse->nBins=0;
      x=0.0;
      do{
        y=(float)gaussian((double)x,(double)gediIO->pSigma,0.0);
        x+=gediIO->pRes;
        gediIO->pulse->nBins+=2;  /*both sides of peak*/
      }while(y>=gediRat->iThresh);

      gediIO->pulse->x=falloc(gediIO->pulse->nBins,"pulse x",0);
      gediIO->pulse->y=falloc(gediIO->pulse->nBins,"pulse y",0);
      gediIO->pulse->centBin=(int)(gediIO->pulse->nBins/2);

      max=-100.0;
      tot=0.0;
      x=-1.0*(float)gediIO->pulse->centBin*gediIO->pRes;
      for(i=0;i<gediIO->pulse->nBins;i++){
        gediIO->pulse->x[i]=x;
        gediIO->pulse->y[i]=(float)gaussian((double)x,(float)gediIO->pSigma,0.0);
        if(gediIO->pulse->y[i]>max){
          max=gediIO->pulse->y[i];
          gediIO->pulse->centBin=i;
        }
        tot+=gediIO->pulse->y[i];
        x+=gediIO->pRes;
      }
      /*normalise to cope with rounding*/
      for(i=0;i<gediIO->pulse->nBins;i++){
        gediIO->pulse->y[i]/=tot;
      }
    }else{  /*dirac-delta*/
      gediIO->pulse->nBins=1;
      gediIO->pulse->x=falloc(gediIO->pulse->nBins,"pulse x",0);
      gediIO->pulse->y=falloc(gediIO->pulse->nBins,"pulse y",0);
      gediIO->pulse->centBin=0;

      gediIO->pulse->x[0]=0.0;
      gediIO->pulse->y[0]=1.0;
    }
  }else{  /*read the pulse from a file*/
    readSimPulse(gediIO,gediRat);
  }

  return;
}/*setGediPulse*/


/*####################################*/
/*read pulse to use for simulator*/

void readSimPulse(gediIOstruct *gediIO,gediRatStruct *gediRat)
{
  int i=0;
  float CofG=0,tot=0,centre=0;
  float minSep=0,max=0;
  char line[400];
  char temp1[100],temp2[100];
  FILE *ipoo=NULL;

  if((ipoo=fopen(gediIO->pulseFile,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",gediIO->pulseFile);
    exit(1);
  }


  /*count number of bins*/
  gediIO->pulse->nBins=0;
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))gediIO->pulse->nBins++;

  gediIO->pulse->x=falloc(gediIO->pulse->nBins,"pulse x",0);
  gediIO->pulse->y=falloc(gediIO->pulse->nBins,"pulse y",0);

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
        gediIO->pulse->x[i]=atof(temp1);
        gediIO->pulse->y[i]=atof(temp2);
        i++;
      }
    }
  }
  gediIO->pRes=fabs(gediIO->pulse->x[gediIO->pulse->nBins-1]-gediIO->pulse->x[0])/(float)(gediIO->pulse->nBins-1);

  /*determine maximum to centre and total to normalise*/
  tot=0.0;
  CofG=0.0;
  max=-1000.0;
  for(i=0;i<gediIO->pulse->nBins;i++){
    CofG+=gediIO->pulse->x[i]*gediIO->pulse->y[i];
    if(gediIO->pulse->y[i]>max){
      max=gediIO->pulse->y[i];
      centre=gediIO->pulse->x[i];
    }
    tot+=gediIO->pulse->y[i];
  }
  CofG/=tot;
  CofG-=centre;

  /*align pulse*/
  minSep=1000.0;
  gediIO->pSigma=0.0;
  for(i=0;i<gediIO->pulse->nBins;i++){
    gediIO->pulse->x[i]-=centre;

    if(fabs(gediIO->pulse->x[i])<minSep){
      minSep=fabs(gediIO->pulse->x[i]);
      gediIO->pulse->centBin=i;
    }
  }

  /*pulse width*/
  gediIO->pSigma=0.0;
  for(i=0;i<gediIO->pulse->nBins;i++){
    gediIO->pSigma+=(gediIO->pulse->x[i]-CofG)*(gediIO->pulse->x[i]-CofG)*gediIO->pulse->y[i];
  }
  gediIO->pSigma=sqrt(gediIO->pSigma/tot);

  /*now normalise*/
  for(i=0;i<gediIO->pulse->nBins;i++)gediIO->pulse->y[i]/=tot;

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return;
}/*readSimPulse*/


/*####################################*/
/*set GEDI footprint*/

void setGediFootprint(gediRatStruct *gediRat,gediIOstruct *gediIO)
{
  int i=0;
  float totE=0;
  float az=0;
  double tX=0,tY=0;
  void intersectOctree(gediRatStruct *);


  /*footprint width*/
  az=gediRat->lobeAng*M_PI/180.0;  /*convert anlge to radians*/

  /*number of lobes and allocate*/
  if(gediRat->sideLobe==0)gediRat->nLobes=1;
  else                   gediRat->nLobes=7;
  if(!(gediRat->lobe=(lobeStruct *)calloc(gediRat->nLobes,sizeof(lobeStruct)))){
    fprintf(stderr,"error lobeStruct allocation.\n");
    exit(1);
  }

  /*central footprint*/
  i=0;
  gediRat->lobe[i].coord[0]=gediRat->coord[0];
  gediRat->lobe[i].coord[1]=gediRat->coord[1];
  gediRat->lobe[i].fSigma=gediIO->fSigma;
  gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);
  if(gediRat->sideLobe==0)gediRat->lobe[i].E=1.0;
  else{  /*include side lobes*/
    totE=1.0+0.0599+0.0731+0.0317+0.0319+0.0167+0.0163;
    i=0;
    gediRat->lobe[i].E=1.0/totE;

    /*first southern lobe*/
    i=1;
    gediRat->lobe[i].E=0.0731/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]-20.0*sin(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]-20.0*cos(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);

    /*first nothern lobe*/
    i=2;
    gediRat->lobe[i].E=0.0599/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]+20.0*sin(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]+20.0*cos(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);

    /*western lobe*/
    i=3;
    gediRat->lobe[i].E=0.0319/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]-20.0*cos(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]-20.0*sin(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);

    /*eastern lobe*/
    i=4;
    gediRat->lobe[i].E=0.0317/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]+20.0*cos(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]+20.0*sin(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);

    /*second southern lobe*/
    i=5;
    gediRat->lobe[i].E=0.0167/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]-30.0*sin(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]-30.0*cos(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);

    /*second northern lobe*/
    i=6;
    gediRat->lobe[i].E=0.0163/totE;
    gediRat->lobe[i].coord[0]=gediRat->coord[0]+30.0*cos(az);
    gediRat->lobe[i].coord[1]=gediRat->coord[1]+30.0*sin(az);
    gediRat->lobe[i].fSigma=gediIO->fSigma;
    gediRat->lobe[i].maxSepSq=(double)(gediRat->maxSep*gediRat->maxSep);
  }/*side lobe test*/


  /*determine min and max bounds*/
  gediRat->minX=gediRat->minY=100000000000.0;
  gediRat->maxX=gediRat->maxY=-1000000000000.0;
  for(i=0;i<gediRat->nLobes;i++){
    tX=gediRat->lobe[i].coord[0]-sqrt(gediRat->lobe[i].maxSepSq);
    if(tX<gediRat->minX)gediRat->minX=tX;
    tX=gediRat->lobe[i].coord[0]+sqrt(gediRat->lobe[i].maxSepSq);
    if(tX>gediRat->maxX)gediRat->maxX=tX;
    tY=gediRat->lobe[i].coord[1]-sqrt(gediRat->lobe[i].maxSepSq);
    if(tY<gediRat->minY)gediRat->minY=tY;
    tY=gediRat->lobe[i].coord[1]+sqrt(gediRat->lobe[i].maxSepSq);
    if(tY>gediRat->maxY)gediRat->maxY=tY;
  }

  /*grid for normalising sampling*/
  if(gediRat->normCover||gediRat->checkCover){
    gediRat->gridRes=1.5;
    gediRat->gX=(int)((float)(gediRat->maxX-gediRat->minX)/gediRat->gridRes)+2;
    gediRat->gY=(int)((float)(gediRat->maxY-gediRat->minY)/gediRat->gridRes)+2;
    gediRat->g0[0]=gediRat->minX+(double)gediRat->gridRes;
    gediRat->g0[1]=gediRat->minY+(double)gediRat->gridRes;
    gediRat->nGrid=ialloc(gediRat->gX*gediRat->gY,"nGrid",0);
    for(i=gediRat->gX*gediRat->gY-1;i>=0;i--)gediRat->nGrid[i]=0;
  }

  /*radius to calculate density within*/
  if(gediRat->topHat==0)gediRat->denseRadSq=gediIO->fSigma*gediIO->fSigma*4.0;
  else                  gediRat->denseRadSq=gediIO->fSigma;
  gediRat->pointDense=gediRat->beamDense=0.0;

  /*determine which octree cells are intresected*/
  if(gediRat->useOctree)intersectOctree(gediRat);

  return;
}/*setGediFootprint*/


/*#####################################################*/
/*determine which top level octree pixels intersect*/

void intersectOctree(gediRatStruct *gediRat)
{
  int i=0,j=0;
  int minI=0,maxI=0;
  int minJ=0,maxJ=0;
  int *markInt(int,int *,int);

  /*reset counters*/
  gediRat->nOct=0;
  TIDY(gediRat->octList);
  gediRat->octList=NULL;

  /*determine bounds to search*/
  minI=(int)((gediRat->minX-gediRat->octree->minX)/(double)gediRat->octree->res);
  maxI=(int)((gediRat->maxX-gediRat->octree->minX)/(double)gediRat->octree->res+0.5);
  minJ=(int)((gediRat->minY-gediRat->octree->minY)/(double)gediRat->octree->res);
  maxJ=(int)((gediRat->maxY-gediRat->octree->minY)/(double)gediRat->octree->res+0.5);

  /*loop over and test*/
  for(i=minI;i<=maxI;i++){
    if((i<0)||(i>=gediRat->octree->nX))continue;
    for(j=minJ;j<=maxJ;j++){
      if((j<0)||(j>=gediRat->octree->nY))continue;
      gediRat->octList=markInt(gediRat->nOct,gediRat->octList,i+j*gediRat->octree->nX);
      gediRat->nOct++;
    }
  }
  return;
}/*intersectOctree*/


/*##############################################*/
/*update footprint cordinate*/

void updateGediCoord(gediRatStruct *gediRat,int i,int j)
{

  if(gediRat->doGrid){
    gediRat->coord[0]=gediRat->gMinX+(double)i*gediRat->gRes;
    gediRat->coord[1]=gediRat->gMinY+(double)j*gediRat->gRes;
  }else if(gediRat->readALSonce){
    gediRat->coord[0]=gediRat->coords[i][0];
    gediRat->coord[1]=gediRat->coords[i][1];
  }

  return;
}/*updateGediCoord*/


/*####################################*/
/*allocate wave structure*/

waveStruct *allocateGEDIwaves(gediIOstruct *gediIO,gediRatStruct *gediRat,pCloudStruct **data,pointMapStruct *pointmap)
{
  int j=0,numb=0,k=0;
  uint32_t i=0,n=0;
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

  for(n=0;n<pointmap->nPoints;n++){
    numb=pointmap->fList[n];
    i=pointmap->pList[n];
    if(data[numb]->nPoints>0)hasPoints=1;
    if(data[numb]->z[i]>maxZ)maxZ=data[numb]->z[i];
    if(data[numb]->z[i]<minZ)minZ=data[numb]->z[i];
  }/*bound finding*/

  if(hasPoints==0){
    fprintf(stderr,"No points included\n");
    exit(1);
  }

  waves->minZ=minZ-buff;
  waves->maxZ=maxZ+buff;

  waves->nBins=(int)((waves->maxZ-waves->minZ)/(double)gediIO->res);
  if(gediRat->readWave)waves->nWaves=9;
  else                 waves->nWaves=3;
  waves->wave=fFalloc(waves->nWaves,"result waveform",0);
  for(j=0;j<waves->nWaves;j++){
    waves->wave[j]=falloc(waves->nBins,"result waveform",j+1);
    for(k=0;k<waves->nBins;k++)waves->wave[j][k]=0.0;
  }
  if(gediIO->ground){
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


/*################################################################################*/
/*determine ALS coverage*/

void determineALScoverage(gediIOstruct *gediIO,gediRatStruct *gediRat,pCloudStruct **data,pointMapStruct *pointmap)
{
  int i=0;
  int gX=0,gY=0;
  uint32_t j=0,k=0;
  double dx=0,dy=0;
  double sepSq=0;
  float area=0.0;
  char hasPoints=0;


  /*reset counter*/
  for(i=gediRat->gX*gediRat->gX-1;i>=0;i--)gediRat->nGrid[i]=0;

  /*loop over points needed*/
  for(k=0;k<pointmap->nPoints;k++){
    i=pointmap->fList[k];
    j=pointmap->pList[k];
    /*check within bounds*/
    if((data[i]->x[j]>=gediRat->minX)&&(data[i]->x[j]<=gediRat->maxX)&&(data[i]->y[j]>=gediRat->minY)&&(data[i]->y[j]<=gediRat->maxY)){
      dx=data[i]->x[j]-gediRat->coord[0];
      dy=data[i]->y[j]-gediRat->coord[1];
      sepSq=dx*dx+dy*dy;

      /*ground grid for ALS coverage*/
      if(gediRat->normCover||gediRat->checkCover){
        if(data[i]->retNumb[j]==data[i]->nRet[j]){  /*only once per beam*/
          /*mark sampling desnity for normalisation*/
          gX=(int)((data[i]->x[j]-gediRat->g0[0])/(double)gediRat->gridRes);
          gY=(int)((data[i]->y[j]-gediRat->g0[1])/(double)gediRat->gridRes);
          if((gX>=0)&&(gX<gediRat->gX)&&(gY>=0)&&(gY<gediRat->gY)){
            gediRat->nGrid[gY*gediRat->gX+gX]++;
          }
        }
      }/*ground grid for ALS coverage*/

      /*point and beam density*/
      if(sepSq<=gediRat->denseRadSq){
        hasPoints=1;
        gediRat->pointDense+=1.0;
        if(data[i]->retNumb[j]==data[i]->nRet[j])gediRat->beamDense+=1.0;
      }
    }/*bounds check*/
  }/*point loop*/

  area=M_PI*gediRat->denseRadSq;
  gediRat->pointDense/=area;
  gediRat->beamDense/=area;

  if(hasPoints==0){
    for(i=gediRat->gX*gediRat->gX-1;i>=0;i--)gediRat->nGrid[i]=0;
  }

  return;
}/*determineALScoverage*/


/*####################################*/
/*check footprint is covered by ALS*/

void checkFootCovered(gediIOstruct *gediIO,gediRatStruct *gediRat)
{
  int i=0,j=0,nWithin=0;
  int thresh=0,nMissed=0;
  double dX=0,dY=0,sepSq=0;
  double useRad=0,radSq=0;

  useRad=10.0;
  if(useRad>(double)gediIO->fSigma)useRad=(double)gediIO->fSigma;
  radSq=useRad*useRad;

  for(i=0;i<gediRat->gX;i++){
    dX=(double)i*(double)gediRat->gridRes-gediRat->coord[0];
    for(j=0;j<gediRat->gY;j++){
      dY=(double)j*(double)gediRat->gridRes-gediRat->coord[1];
      sepSq=dX*dX+dY*dY;

      if(sepSq<radSq){
        if(gediRat->nGrid[j*gediRat->gX+i]==0)nMissed++;
        nWithin++;
      }
    }/*y loop*/
  }/*x loop*/

  thresh=(int)((float)nWithin*2.0/3.0);
  if(nMissed>thresh){
    fprintf(stderr,"Too many missed %d of %d\n",nMissed,nWithin);
    gediRat->useFootprint=0;
  }else gediRat->useFootprint=1;

  return;
}/*checkFootCovered*/


/*####################################*/
/*set deconvolution parameters for GEDI*/

denPar *setDeconForGEDI(gediRatStruct *gediRat)
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
  decon->pScale=1.0;
  decon->noiseTrack=0;

  /*read system pulse*/
  readPulse(decon);

  return(decon);
}/*setDeconForGEDI*/


/*####################################*/
/*GEDI wave from ALS waveforms*/

void gediFromWaveform(pCloudStruct *data,uint32_t i,float rScale,waveStruct *waves,gediRatStruct *gediRat,gediIOstruct *gediIO)
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
  for(j=0;j<buffBins;j++)temp[j]=(unsigned char)gediRat->meanN;
  for(j=0;j<(int)data->waveLen[i];j++)temp[j+buffBins]=wave[j];
  for(j=(int)data->waveLen[i]+buffBins;j<(int)waveLen;j++)temp[j]=(unsigned char)gediRat->meanN;
  TIDY(wave);
  wave=temp;
  temp=NULL;

  /*deconvolve and reconvolve*/
  if(gediRat->indDecon){
    processed=processWave(wave,(int)waveLen,gediRat->decon,1.0);
    smooPro=smooth(gediIO->pSigma,(int)waveLen,processed,gediIO->res);
  }

  /*convolve with GEDI pulse*/
  floWave=falloc((int)waveLen,"",0);
  for(j=0;j<(int)waveLen;j++)floWave[j]=(float)wave[j]-gediRat->meanN;
  smoothed=smooth(gediIO->pSigma,(int)waveLen,floWave,gediIO->res);
  TIDY(floWave);

  /*add up*/
  for(j=0;j<(int)waveLen;j++){
    binPosition(&x,&y,&z,j-buffBins,data->x[i],data->y[i],data->z[i],data->time[i],&(grad[0]));
    bin=(int)((waves->maxZ-z)/(double)gediIO->res);
    if((bin>=0)&&(bin<waves->nBins)){
      waves->wave[3][bin]+=((float)wave[j]-gediRat->meanN)*rScale;  /*with ALS pulse*/
      waves->wave[4][bin]+=smoothed[j]*rScale;
      if(gediRat->doDecon){
        if(gediRat->indDecon){
          waves->wave[5][bin]+=smooPro[j]*rScale;
          waves->wave[6][bin]+=processed[j]*rScale;
        }
        waves->wave[7][bin]+=((float)wave[j]-gediRat->meanN)*rScale;
      }
    }
  }/*wave bin loop*/

  TIDY(wave);
  TIDY(smooPro);
  TIDY(smoothed);
  TIDY(processed);
  return;
}/*gediFromWaveform*/


/*################################################################################*/
/*make waveform from point cloud*/

void waveFromPointCloud(gediRatStruct *gediRat, gediIOstruct *gediIO,pCloudStruct **data,waveStruct *waves,pointMapStruct *pointmap)
{
  int numb=0,bin=0,j=0;
  int gX=0,gY=0,n=0;
  int xInd=0,yInd=0;
  uint32_t i=0,k=0;
  double sep=0;
  double dX=0,dY=0;
  double totGround=0;     /*contrbution to ground estimate*/
  float refl=0,rScale=0,fracHit=0,totAng=0;
  void gediFromWaveform(pCloudStruct *,uint32_t,float,waveStruct *,gediRatStruct *,gediIOstruct *);

  /*reset mean scan angle*/
  waves->meanScanAng=totAng=0.0;
  /*ground elevation estimate*/
  if(gediIO->ground){
    waves->gElevSimp=0.0;
    totGround=0.0;
  }

  /*make waves*/
  for(n=0;n<gediRat->nLobes;n++){
    for(k=0;k<pointmap->nPoints;k++){
      numb=pointmap->fList[k];
      i=pointmap->pList[k];

      /*determine laser intensity at this point*/
      dX=data[numb]->x[i]-gediRat->lobe[n].coord[0];
      dY=data[numb]->y[i]-gediRat->lobe[n].coord[1];
      if(gediRat->defWfront==0){    /*symmetric wavefront*/
        sep=sqrt(dX*dX+dY*dY);
        if(gediRat->topHat==0)rScale=(float)gaussian(sep,(double)gediRat->lobe[n].fSigma,0.0);
        else{
          if(sep<=gediRat->lobe[n].maxSepSq)rScale=1.0;
          else                             rScale=0.0;
        }
      }else{     /*read assymmetric pulse*/
        xInd=(int)((dX*cos(gediRat->lobeAng)+dY*sin(gediRat->lobeAng))/(double)gediRat->wavefront->res)+gediRat->wavefront->x0;
        yInd=(int)((dY*cos(gediRat->lobeAng)-dX*sin(gediRat->lobeAng))/(double)gediRat->wavefront->res)+gediRat->wavefront->y0;
        if((xInd>=0)&&(xInd<gediRat->wavefront->nX)&&(yInd>=0)&&(yInd<gediRat->wavefront->nY)){
          if(gediRat->wavefront->front[xInd][yInd]>0.0)rScale=gediRat->wavefront->front[xInd][yInd];
          else rScale=0.0;
        }else rScale=0.0;
      }/*determine laser intensity at this point*/

      if(rScale>gediRat->iThresh){  /*if bright enough to matter*/
        /*scale by sampling density*/
        if(gediRat->normCover){
          gX=(int)((data[numb]->x[i]-gediRat->g0[0])/(double)gediRat->gridRes);
          gY=(int)((data[numb]->y[i]-gediRat->g0[1])/(double)gediRat->gridRes);
          if((gX>=0)&&(gX<gediRat->gX)&&(gY>=0)&&(gY<gediRat->gY)){
            if(gediRat->nGrid[gY*gediRat->gX+gX]>0)rScale/=(float)gediRat->nGrid[gY*gediRat->gX+gX];
          }
        }/*scale by sampling density*/


        /*discrete return*/
        refl=(float)data[numb]->refl[i]*rScale;
        if(data[numb]->nRet[i]>0)fracHit=1.0/(float)data[numb]->nRet[i];
        else                     fracHit=1.0;
        for(j=0;j<gediIO->pulse->nBins;j++){
          bin=(int)((waves->maxZ-data[numb]->z[i]+(double)gediIO->pulse->x[j])/(double)gediIO->res);
          if((bin>=0)&&(bin<waves->nBins)){
            waves->wave[0][bin]+=refl*gediIO->pulse->y[j];
            waves->wave[1][bin]+=rScale*gediIO->pulse->y[j];
            waves->wave[2][bin]+=rScale*fracHit*gediIO->pulse->y[j];
            if(gediIO->ground){
              if(data[numb]->class[i]==2){
                waves->ground[0][bin]+=refl*gediIO->pulse->y[j];
                waves->ground[1][bin]+=rScale*gediIO->pulse->y[j];
                waves->ground[2][bin]+=rScale*fracHit*gediIO->pulse->y[j];
              }else{
                waves->canopy[0][bin]+=refl*gediIO->pulse->y[j];
                waves->canopy[1][bin]+=rScale*gediIO->pulse->y[j];
                waves->canopy[2][bin]+=rScale*fracHit*gediIO->pulse->y[j];
              }
            }/*ground recording if needed*/
          }/*bin bound check*/
        }/*pulse bin loop*/
        if(gediIO->ground){
          if(data[numb]->class[i]==2){
            waves->gElevSimp+=rScale*data[numb]->z[i];
            totGround+=rScale;
          }
        }
        waves->meanScanAng+=rScale*fracHit*(float)abs((int)data[numb]->scanAng[i]);
        totAng+=rScale*fracHit;

        /*full-waveform*/
        if(gediRat->readWave&&data[numb]->hasWave){
          if(data[numb]->packetDes[i]){  /*test for waveform*/
            gediFromWaveform(data[numb],i,rScale,waves,gediRat,gediIO);
          }
        }/*waveform test*/
      }
    }/*point loop*/
  }/*lobe loop*/

  /*normalise mean scan angle*/
  if(totAng>0.0)waves->meanScanAng/=totAng;
  if(totGround>=0.0)waves->gElevSimp/=totGround;
  else              waves->gElevSimp=-9999.0;

  return;
}/*waveFromPointCloud*/


/*################################################################################*/
/*make a map of voxel gaps*/

void voxelGap(gediRatStruct *gediRat,gediIOstruct *gediIO,pCloudStruct **data,waveStruct *waves)
{
  int i=0,vInd=0;
  int xBin=0,yBin=0,zBin=0;
  uint32_t j=0;
  double bounds[6];
  voxStruct *vox=NULL;
  voxStruct *tidyVox(voxStruct *);
  voxStruct *voxAllocate(int,float *,double *,char);
  void countVoxGap(double,double,double,float *,voxStruct *,int,int,float,int);


  bounds[0]=gediRat->minX;
  bounds[1]=gediRat->minY;
  bounds[2]=waves->minZ;
  bounds[3]=gediRat->maxX;
  bounds[4]=gediRat->maxY;
  bounds[5]=waves->maxZ;


  /*first make a voxel map*/
  vox=voxAllocate(1,&(gediRat->vRes[0]),&(bounds[0]),0);

  for(i=0;i<gediIO->nFiles;i++){ /*file loop*/
    for(j=0;j<data[i]->nPoints;j++){  /*point loop*/
      countVoxGap(data[i]->x[j],data[i]->y[j],data[i]->z[j],&(data[i]->grad[j][0]),vox,1,1,gediRat->beamRad,0);
    }/*point loop*/
  }/*file loop*/

  /*calculate gap fraction for each return*/
  for(i=0;i<gediIO->nFiles;i++){ /*file loop*/
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
/*make waveforms accounting for shadowing*/

void waveFromShadows(gediRatStruct *gediRat,gediIOstruct *gediIO,pCloudStruct **data,waveStruct *waves,pointMapStruct *pointmap)
{
  int i=0;
  float **tempWave=NULL;
  //float iRes=0,grad[3];
  float grad[3];
  void voxelGap(gediRatStruct *,gediIOstruct *,pCloudStruct **,waveStruct *);
  rImageStruct *rImage=NULL;    /*range image, a stack nBins long*/
  lidVoxPar lidPar;

  fprintf(stderr,"Silouhette images do not currently wqork with octrees\n");
  exit(1);

  /*iRes=0.02;*/
  grad[0]=grad[1]=0.0;
  grad[2]=-1.0;

  /*set lidar parameters for a downwards looking ALS*/
  lidPar.minRefl=1.0;             /*minimum refletance value to scale between 0 and 1*/
  lidPar.maxRefl=1.0;             /*maximum refletance value to scale between 0 and 1*/
  lidPar.appRefl=1.0 ;            /*scale between TLS reflectance and size*/
  lidPar.beamTanDiv=0.0;          /*tan of beam divergence*/
  lidPar.beamRad=gediRat->beamRad; /*start radius*/
  lidPar.minGap=0.00001;          /*minimum gap fraction correction to apply*/

  /*gap fraction from voxelising data*/
  voxelGap(gediRat,gediIO,data,waves);

  /*create images*/
  /*rImage=allocateRangeImage(gediIO->nFiles,data,gediIO->pRes*4.0,iRes,&(grad[0]),gediRat->coord[0],gediRat->coord[1],waves->maxZ);*/
  /*rImage=allocateRangeImage(gediIO->nFiles,data,NULL,0.15,0.01,&(grad[0]),gediRat->coord[0],gediRat->coord[1],waves->maxZ,NULL);*/
  fprintf(stderr,"THis method is no longer operational. Do not use\n");
  exit(1);


  silhouetteImage(gediIO->nFiles,data,NULL,rImage,&lidPar,NULL,0,NULL);


  /*convert images to waveform*/
  tempWave=fFalloc(2,"",0);
  for(i=0;i<2;i++)tempWave[i]=falloc(rImage->nBins,"",i+1);
  waveFromImage(rImage,tempWave,1,gediIO->fSigma);
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


/*####################################*/
/*clean outlier points from waveform*/

void cleanOutliers(waveStruct *waves,gediIOstruct *gediIO)
{
  int i=0,j=0,gStart=0;
  char pastGround=0;
  float gGap=0;  /*gap in ground return*/
  float maxGap=0;
  float maxGround=0,gThresh=0;
  float max=0,thresh=0;

  if(!gediIO->ground){
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
      if(pastGround)gGap+=gediIO->res;
    }
    if(gGap>maxGap){  /*too big a break, delete*/
      waves->groundBreakElev=waves->maxZ-(double)i*(double)gediIO->res;
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
    gGap+=gediIO->res;

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
/*deconvolve aggragated wave*/

void processAggragate(gediRatStruct *gediRat,gediIOstruct *gediIO,waveStruct *waves)
{
  int i=0;
  float *processed=NULL,*smooPro=NULL;
  float *smooth(float,int,float *,float);

  /*Add background noise back*/
  for(i=0;i<waves->nBins;i++)waves->wave[7][i]+=gediRat->meanN;

  /*deconvolve and reconvolve*/
  processed=processFloWave(&(waves->wave[7][0]),(int)waves->nBins,gediRat->decon,1.0);
  for(i=0;i<waves->nBins;i++)waves->wave[8][i]=processed[i];
  smooPro=smooth(gediIO->pSigma,(int)waves->nBins,processed,gediIO->res);
  TIDY(processed);

  TIDY(waves->wave[7]);
  waves->wave[7]=smooPro;
  smooPro=NULL;

  return;
}/*processAggragate*/


/*####################################*/
/*make GEDI waveforms*/

waveStruct *makeGediWaves(gediRatStruct *gediRat,gediIOstruct *gediIO,pCloudStruct **data)
{
  int j=0,k=0;
  float tot=0;
  waveStruct *waves=NULL;
  waveStruct *allocateGEDIwaves(gediIOstruct *,gediRatStruct *,pCloudStruct **,pointMapStruct *);
  void processAggragate(gediRatStruct *,gediIOstruct *,waveStruct *);
  void checkFootCovered(gediIOstruct *,gediRatStruct *);
  void cleanOutliers(waveStruct *,gediIOstruct *);
  void waveFromPointCloud(gediRatStruct *,gediIOstruct *,pCloudStruct **,waveStruct *,pointMapStruct *);
  void waveFromShadows(gediRatStruct *,gediIOstruct *,pCloudStruct **,waveStruct *,pointMapStruct *);
  void determineALScoverage(gediIOstruct *,gediRatStruct *,pCloudStruct **,pointMapStruct *);
  denPar *setDeconForGEDI(gediRatStruct *);
  pointMapStruct *findIntersectingMap(gediRatStruct *,gediIOstruct *,pCloudStruct **);
  pointMapStruct *pointmap=NULL;


  /*determine list of opints to use*/
  pointmap=findIntersectingMap(gediRat,gediIO,data);

  /*determine ALS coverage*/
  determineALScoverage(gediIO,gediRat,data,pointmap);

  /*check that whole footprint is covered*/
  if(pointmap->nPoints==0)gediRat->useFootprint=0;
  else if(gediRat->checkCover)checkFootCovered(gediIO,gediRat);
  else                   gediRat->useFootprint=1;

  /*only if it contains data*/
  if(gediRat->useFootprint){
    /*allocate*/
    waves=allocateGEDIwaves(gediIO,gediRat,data,pointmap);

    /*set up denoising if using*/
    if(gediRat->doDecon)gediRat->decon=setDeconForGEDI(gediRat);

    /*make waves*/
    if(gediRat->useShadow==0)waveFromPointCloud(gediRat,gediIO,data,waves,pointmap);
    else                    waveFromShadows(gediRat,gediIO,data,waves,pointmap);

    /*clean outliers if needed*/
    if(gediRat->cleanOut)cleanOutliers(waves,gediIO);
    else                waves->groundBreakElev=-100000000.0;

    /*deconvolve aggragated waveform*/
    if(gediRat->doDecon)processAggragate(gediRat,gediIO,waves);

    /*normalise integral*/
    for(k=0;k<waves->nWaves;k++){
      tot=0.0;
      for(j=0;j<waves->nBins;j++)tot+=waves->wave[k][j]*gediIO->res;
      if(tot>0.0){
        for(j=0;j<waves->nBins;j++)waves->wave[k][j]/=tot;
        if(gediIO->ground&&(k<3)){
          for(j=0;j<waves->nBins;j++){
            waves->canopy[k][j]/=tot;
            waves->ground[k][j]/=tot;
          }
        }
      }
    }

    /*tidy arrays*/
    if(gediRat->decon){
      TTIDY((void **)gediRat->decon->pulse,2);
      gediRat->decon->pulse=NULL;
      TIDY(gediRat->decon);
    }
    tot=0.0;
    for(j=0;j<waves->nBins;j++)tot+=waves->wave[0][j]*gediIO->res;
  }/*contains data*/

  /*check whether empty*/
  if((tot<TOL)||(waves->nBins==0))gediRat->useFootprint=0;
  if(pointmap){
    TIDY(pointmap->fList);
    TIDY(pointmap->pList);
    TIDY(pointmap);
  }

  return(waves);
}/*makeGediWaves*/


/*####################################################*/
/*map points and file indices to use*/

pointMapStruct *findIntersectingMap(gediRatStruct *gediRat,gediIOstruct *gediIO,pCloudStruct **data)
{
  int i=0;
  uint32_t j=0,ind=0;
  pointMapStruct *pointmap=NULL;

  /*search octree or copy all points*/
  if(gediRat->useOctree){
    pointmap=mapFromOctree(gediRat->octList,gediRat->nOct,gediRat->octree,gediRat->minX,gediRat->maxX,gediRat->minY,gediRat->maxY);
  }else{   /*use all points*/
    /*allocate space*/
    if(!(pointmap=(pointMapStruct *)calloc(1,sizeof(pointMapStruct)))){
      fprintf(stderr,"error pointMapStruct allocation.\n");
      exit(1);
    }
    pointmap->nPoints=0;
    pointmap->fList=NULL;
    pointmap->pList=NULL;

    for(i=0;i<gediIO->nFiles;i++){
      if(data[i]->nPoints==0)continue;
      if(!(pointmap->fList=(int *)realloc(pointmap->fList,(pointmap->nPoints+data[i]->nPoints)*sizeof(int)))){
        fprintf(stderr,"Error allocating memory\n");
        exit(1);
      }
      if(!(pointmap->pList=(uint32_t *)realloc(pointmap->pList,(pointmap->nPoints+data[i]->nPoints)*sizeof(uint32_t)))){
        fprintf(stderr,"Error allocating memory\n");
        exit(1);
      }
      for(j=0;j<+data[i]->nPoints;j++){
        ind=pointmap->nPoints+j;
        pointmap->fList[ind]=i;
        pointmap->pList[ind]=j;
      }
      pointmap->nPoints+=data[i]->nPoints;
    }
  }
  return(pointmap);
}/*findIntersectingMap*/


/*####################################################*/
/*allocate space and copy wavefront filename*/

wFrontStruct *copyFrontFilename(char *namen)
{
  wFrontStruct *wavefront=NULL;

  if(!(wavefront=(wFrontStruct *)calloc(1,sizeof(wFrontStruct)))){
    fprintf(stderr,"error in wavefront allocation.\n");
    exit(1);
  }

  wavefront->front=NULL;
  strcpy(wavefront->frontFile,namen);

  return(wavefront);
}/*copyFrontFilename*/


/*the end*/
/*####################################################*/

