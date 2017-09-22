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

dataStruct *unpackHDFlvis(char *namen,lvisHDF *hdfLvis,gediIOstruct *gediIO,int numb)
{
  int i=0;
  dataStruct *data=NULL;
  float *tempPulse=NULL;
  float pulseLenFromTX(float *,int);

  /*read data if needed*/
  if(hdfLvis==NULL){
    hdfLvis=readLVIShdf(namen);
    gediIO->nFiles=hdfLvis->nWaves;
    gediIO->ground=0;
  }

  /*allocate space*/
  if(!(data=(dataStruct *)calloc(1,sizeof(dataStruct)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  data->useID=1;
  data->nBins=hdfLvis->nBins;
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
  data->zen=hdfLvis->zen[numb];
  data->res=fabs(hdfLvis->z0[numb]-hdfLvis->z1023[numb])/(float)hdfLvis->nBins;
  if(gediIO->den)gediIO->den->res=data->res;
  if(gediIO->gFit)gediIO->gFit->res=data->res;
  data->totE[data->useType]=0.0;
  for(i=0;i<hdfLvis->nBins;i++){
    data->wave[data->useType][i]=(float)hdfLvis->wave[numb][i];
    data->z[i]=(double)(hdfLvis->z0[numb]-(float)i*data->res);
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
  data->lon=(hdfLvis->lon0[numb]+hdfLvis->lon1023[numb])/2.0;
  data->lat=(hdfLvis->lat0[numb]+hdfLvis->lat1023[numb])/2.0;
  data->lfid=hdfLvis->lfid[numb];
  data->shotN=hdfLvis->shotN[numb];
  sprintf(data->waveID,"%d.%d",hdfLvis->lfid[numb],hdfLvis->shotN[numb]);

  /*analyse pulse*/
  if(gediIO->readPsigma){
    tempPulse=falloc(hdfLvis->pBins,"temp pulse",0);
    for(i=0;i<hdfLvis->pBins;i++)tempPulse[i]=(float)hdfLvis->pulse[numb][i];
    data->pSigma=pulseLenFromTX(tempPulse,hdfLvis->pBins);
    TIDY(tempPulse);
  }else data->pSigma=gediIO->pSigma;

  if(data->fSigma<0.0)data->fSigma=gediIO->fSigma;
  if(data->pSigma<0.0)data->pSigma=gediIO->pSigma;

  /*set up number of messages*/
  if((gediIO->nMessages>0)&&(hdfLvis->nWaves>gediIO->nMessages))gediIO->nMessages=(int)(hdfLvis->nWaves/gediIO->nMessages);
  else                                 gediIO->nMessages=1;

  return(data);
}/*unpackHDFlvis*/


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
    for(i=0;i<nBins;i++)pSigma+=((float)i*den.res*denoised[i]-CofG)*((float)i*den.res*denoised[i]-CofG);
    pSigma=sqrt(pSigma/tot);

  }else pSigma=-1.0;

  TIDY(denoised);
  return(pSigma);
}/*pulseLenFromTX*/


/*####################################################*/
/*read LVIS binary data*/

dataStruct *readBinaryLVIS(char *namen,lvisLGWstruct *lvis,int numb,gediIOstruct *gediIO)
{
  int i=0;
  dataStruct *data=NULL;
  float *tempPulse=NULL;
  float pulseLenFromTX(float *,int);

  /*do we need to read all the data*/
  if(lvis->data==NULL){
    lvis->verMaj=1; /*dimage->verMaj;*/
    lvis->verMin=3; /*dimage->verMin;*/
    /*check version number*/
    if((lvis->verMaj!=1)||(lvis->verMin!=3)){
      fprintf(stderr,"Currently only works for version 1.3, not %d.%d\n",lvis->verMaj,lvis->verMin);
      exit(1);
    }
    /*read data*/
    lvis->nBins=432;
    lvis->data=readLVISlgw(namen,&lvis->nWaves);
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
  data->res=gediIO->den->res=gediIO->gFit->res=(lvis->data[numb].z0-lvis->data[numb].z431)/(float)lvis->nBins;
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
  data->lon=(lvis->data[numb].lon0+lvis->data[numb].lon431)/2.0;
  data->lat=(lvis->data[numb].lat0+lvis->data[numb].lat431)/2.0;
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
  if(lvis->nWaves>gediIO->nMessages)gediIO->nMessages=(int)(lvis->nWaves/gediIO->nMessages);
  else                              gediIO->nMessages=1;

  return(data);
}/*readBinaryLVIS*/


/*####################################*/
/*read lasFile and save relevant data*/

pCloudStruct *readALSdata(lasFile *las,gediRatStruct *gediRat)
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
  else                    useFile=checkMultiFiles(las,gediRat->gNx,gediRat->coords,gediRat->maxSep);

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

  /*footprint width*/
fprintf(stdout,"bits %f %f\n",gediIO->fSigma,gediRat->iThresh);
  if(gediRat->topHat==0)gediRat->maxSep=determineGaussSep(gediIO->fSigma,gediRat->iThresh);
  else                  gediRat->maxSep=gediIO->fSigma;

  if(gediRat->doGrid){  /*it is a grid*/
    /*number of footprints*/
    gediRat->gNx=(int)((gediRat->gMaxX-gediRat->gMinX)/gediRat->gRes+1);
    gediRat->gNy=(int)((gediRat->gMaxY-gediRat->gMinY)/gediRat->gRes+1);

    /*global bounds*/
    gediRat->globMinX=gediRat->gMinX-gediRat->maxSep;
    gediRat->globMaxX=gediRat->gMaxX+gediRat->maxSep;
    gediRat->globMinY=gediRat->gMinY-gediRat->maxSep;
    gediRat->globMaxY=gediRat->gMaxY+gediRat->maxSep;
  }else if(gediRat->coords){  /*we have provided a list of coordinates*/
    setRatBounds(gediRat);
  }else if(gediRat->readALSonce){ /*it is a batch*/
    /*read list of coords*/
    readFeetList(gediRat);
  }else{   /*single footprint*/
    gediRat->gNx=gediRat->gNy=1;
    gediRat->globMinX=gediRat->coord[0]-gediRat->maxSep;
    gediRat->globMaxX=gediRat->coord[0]+gediRat->maxSep;
    gediRat->globMinY=gediRat->coord[1]-gediRat->maxSep;
    gediRat->globMaxY=gediRat->coord[1]+gediRat->maxSep;
  }

  if((gediIO->nMessages>1)&&(gediRat->gNx*gediRat->gNy)>gediIO->nMessages)gediIO->nMessages=(int)(gediRat->gNx*gediRat->gNy/gediIO->nMessages);
  else                                             gediIO->nMessages=1;

  return;
}/*setGediGrid*/


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
}/*readFeetList*/


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
  int i=0,maxBin=0;
  float CofG=0,tot=0,centre=0;
  float minSep=0,max=0;
  float wThresh=0;
  float p0=0,p1=0;
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
      maxBin=i;
    }
    tot+=gediIO->pulse->y[i];
  }
  CofG/=tot;

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
  wThresh=max*exp(-0.5);
  minSep=1000000.0;
  for(i=0;i<gediIO->pulse->nBins;i++){
    if(i<maxBin){
      if(fabs(gediIO->pulse->y[i]-wThresh)<minSep){
        minSep=fabs(gediIO->pulse->y[i]-wThresh);
        p0=gediIO->pulse->x[i];
      }
    }else if(i==maxBin){
      minSep=1000000.0;
    }else{
      if(fabs(gediIO->pulse->y[i]-wThresh)<minSep){
        minSep=fabs(gediIO->pulse->y[i]-wThresh);
        p1=gediIO->pulse->x[i];
      }
    }
  }
  gediIO->pSigma=fabs(p1-p0)/2.0;

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

  return;
}/*setGediFootprint*/


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

/*the end*/
/*####################################################*/

