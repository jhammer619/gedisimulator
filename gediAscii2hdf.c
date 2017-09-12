#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "hdf5.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"
#include "libLidarHDF.h"
#include "gediIO.h"



/*tolerances*/
#define TOL 0.00001



/*####################################*/
/*control structure*/

typedef struct{
  /*input/output*/
  int nFiles;   /*number of waveforms*/
  char **inList;
  char listNamen[200];
  char outNamen[200];
  FILE *opooGauss;  /*Gaussian parameter output*/
  FILE *opooMet;    /*waveform metric output*/
  int maxGauss;     /*maximum number of Gaussians for output*/

  /*level2 LVIS for ZG*/
  char l2namen[200]; /*list of level2 filenames*/
  char readL2;      /*switch to read L2 or not*/

  /*switches*/
  char ground;      /*read separateground wave or not*/
  char writeFit;    /*write fitted wave switch*/
  char useInt;      /*use discrete intensity instead of count*/
  char useFrac;     /*use fraction of hits per beam for weighting*/
  float rhRes;      /*rh resolution*/
  char bayesGround; /*Bayseian ground finding*/
  char noRHgauss;   /*do not do Gaussian fitting*/
  char renoiseWave; /*remove noise before adding*/
  char dontTrustGround; /*don't trust ground included with waveforms*/
  char readBinLVIS;  /*read binary LVIS rather than a list of ASCII files*/
  char readHDFlvis;  /*read HDF5 LVIS rather than ASCII*/
  char readPsigma;   /*read psigma from files or not*/
  char coord2dp;     /*round up coords to 2dp when writing*/
  char useBounds;    /*when we will process only a subset of bounds*/

  /*denoising parameters*/
  denPar *den;   /*for denoising*/
  denPar *gFit;  /*for Gaussian fitting*/

  /*noise parameters*/
  float meanN;
  float nSig;
  float bThresh;   /*bounds threshold*/
  float hNoise;    /*hard threshold noise as a fraction of integral*/
  char missGround; /*force to miss ground to get RH errors*/
  float minGap;    /*minimum detectavle gap fraction*/
  char linkNoise;  /*use link noise or not*/
  float linkM;     /*link margin*/
  float linkCov;   /*cover at which link margin is defined*/
  float linkSig;   /*link noise sigma*/
  float linkFsig;  /*footprint sigma used for link margin*/
  float linkPsig;  /*pulse sigma used for link margin*/
  float trueSig;   /*true noise sigma in DN*/
  float deSig;     /*detector sigma*/
  char bitRate;   /*digitiser bit rate*/
  float maxDN;    /*maximum DN we need to digitise*/
  float offset;   /*waveform DN offset*/

  /*pulse parameters*/
  float pSigma;    /*pulse length*/
  float fSigma;    /*footprint width*/
  float newPsig;   /*new pulse sigma*/

  /*LVIS data*/
  int verMaj;           /*major version*/
  int verMin;           /*minor version*/
  //lvisLGWstruct lvis;   /*LVIS lgw structure*/
  lvisHDF *hdf;         /*LVIS HDF5 structure*/
  //lvisL2struct *lvisL2; /*LVIS level2 data*/

  /*bounds for subsets*/
  double minX;
  double maxX;
  double minY;
  double maxY;

  /*others*/
  float rhoRatio; /*ration of canopy to ground reflectance*/
  float res;      /*range resolution*/
  float gTol;     /*toleranve used to label ALS ground finding*/
  float zen;      /*zenith angle*/
  int nMessages;  /*number of progress messages*/
}control;


/*#####################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  dataStruct **data=NULL;
  dataStruct **readMultiData(control *);
  gediHDF *hdfData=NULL;
  gediHDF *arrangeGEDIhdf(dataStruct **,control *);


  /*read command line*/
  dimage=readCommands(argc,argv);

  /*read data*/
  data=readMultiData(dimage);

  /*copy into HDF structure*/
  hdfData=arrangeGEDIhdf(data,dimage);

  /*tidy data array*/
  data=tidyAsciiStruct(data,dimage->nFiles);

  /*write HDF5*/
  writeGEDIhdf(hdfData,dimage->outNamen);

  /*tidy arrays*/
  hdfData=tidyGediHDF(hdfData);
  if(dimage){
    TTIDY((void **)dimage->inList,dimage->nFiles);
    TIDY(dimage->gFit);
    TIDY(dimage->den);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################################*/
/*turn data into HDF5 structure*/

gediHDF *arrangeGEDIhdf(dataStruct **data,control *dimage)
{
  gediHDF *hdfData=NULL;
  void trimDataLength(dataStruct **,gediHDF *);

  /*allocate space for all*/
  if(!(hdfData=(gediHDF *)calloc(1,sizeof(gediHDF)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  hdfData->nWaves=dimage->nFiles;
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


/*####################################################*/
/*read multiple data files*/

dataStruct **readMultiData(control *dimage)
{
  int i=0;
  dataStruct **data=NULL;
  dataStruct *readASCIIdata(char *,control *);

  /*allocate space for all*/
  if(!(data=(dataStruct **)calloc(dimage->nFiles,sizeof(dataStruct *)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*read data*/
  for(i=0;i<dimage->nFiles;i++)data[i]=readASCIIdata(dimage->inList[i],dimage);

  return(data);
}/*readMultiData*/


/*####################################################*/
/*read ASCII data*/

dataStruct *readASCIIdata(char *namen,control *dimage)
{
  int i=0;
  dataStruct *data=NULL;
  char line[1000],temp1[100],temp2[100];
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
  while(fgets(&(line[0]),1000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      data->nBins++;
    }
  }

  data->wave=falloc(data->nBins,"waveform",0);
  data->z=dalloc(data->nBins,"z",0);
  if(dimage->ground)data->ground=falloc(data->nBins,"ground",0);
  data->useID=0;
  data->demGround=0;

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read data*/
  i=0;
  while(fgets(line,1000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(dimage->useInt){  /*read intensity*/
        if(dimage->ground==0){   /*don't read ground*/
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
      }else if(dimage->useFrac){ /*read fraction*/
        if(dimage->ground==0){   /*don't read ground*/
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
        if(dimage->ground==0){   /*don't read ground*/
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
          if(!dimage->dontTrustGround){
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

  if(data->res<=0.0)data->res=dimage->res;

  /*add up energy*/
  data->totE=0.0;
  for(i=0;i<data->nBins;i++)data->totE+=data->wave[i];


  dimage->res=dimage->den->res=dimage->gFit->res=fabs(data->z[1]-data->z[0]);
  if(dimage->den->res<TOL)data->usable=0;
  if(data->totE<=0.0)data->usable=0;
  if(dimage->ground==0){   /*set to blank*/
    data->cov=-1.0;
    data->gLap=-1.0;
    data->gMinimum=-1.0;
    data->gInfl=-1.0;
  }
  if(data->fSigma<0.0)data->fSigma=dimage->fSigma;
  if(data->pSigma<0.0)data->pSigma=dimage->pSigma;

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  /*set up number of messages*/
  //if(dimage->nFiles>dimage->nMessages)dimage->nMessages=(int)(dimage->nFiles/dimage->nMessages);
  //else                                dimage->nMessages=1;

  return(data);
}/*readASCIIdata*/


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
  if(!(dimage->den=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  if(!(dimage->gFit=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }


  /*defaults*/
  /*input/output*/
  dimage->nFiles=1;
  dimage->inList=chChalloc(dimage->nFiles,"inList",0);
  dimage->inList[0]=challoc(200,"inList",0);
  strcpy(&(dimage->inList[0][0]),"/Users/stevenhancock/data/gedi/simulations/USDA_CO/gedi.USDA_CO.4112.wave");
  strcpy(&(dimage->outNamen[0]),"teast.h5");
  dimage->ground=1;


  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=NULL;
        dimage->nFiles=1;
        dimage->inList=chChalloc(dimage->nFiles,"input name list",0);
        dimage->inList[0]=challoc((uint64_t)strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=readInList(&dimage->nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to convert ASCII GEDI waveforms to HDF5\n#####\n\n-input name;     single input filename\n-output name;    output filename\n-inList list;    input file list for multiple files\n\n");
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
/*#####################################*/

