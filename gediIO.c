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
  data->res=gediIO->den->res=gediIO->gFit->res=fabs(hdfLvis->z0[numb]-hdfLvis->z1023[numb])/(float)hdfLvis->nBins;
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
  if(hdfLvis->nWaves>gediIO->nMessages)gediIO->nMessages=(int)(hdfLvis->nWaves/gediIO->nMessages);
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

/*the end*/
/*####################################################*/

