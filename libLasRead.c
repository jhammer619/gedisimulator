#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libLasRead.h" 


/*######################*/
/*# A library for      #*/
/*# handling las files #*/
/*# S Hancock, 2015    #*/
/*######################*/

int nFilesOpen;    /*record the number of files open*/
#define MAX_OPEN 250


/*##################################################################*/
/*read the jheader*/

lasFile *readLasHead(char *namen,uint64_t pBuffSize)
{
  int i=0;
  int offset=0;         /*to step around header byte arrays*/
  int tempLen=0;
  char *pubHead=NULL;   /*public header*/
  lasFile *las=NULL;

  if(!(las=(lasFile *)calloc(1,sizeof(lasFile)))){
    fprintf(stderr,"error lasFile allocation.\n");
    exit(1);
  }

  /*open file*/
  if((las->ipoo=fopen(namen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }
  nFilesOpen++;
  strcpy(las->namen,namen);


  /*determine type*/
  tempLen=96;
  pubHead=challoc(tempLen,"pubHead",0);
  if(fread(&(pubHead[0]),sizeof(char),tempLen,las->ipoo)!=tempLen){
    fprintf(stderr,"error reading data from %s\n",namen);
    exit(1);
  }
  if(strncasecmp(pubHead,"LASF",4)){  /*check is a las file*/
    fprintf(stderr,"Incorrect filetype for %s\n",namen);
    exit(1);
  }/*check it is a LASF file*/

  offset=24;  /*version major*/
  memcpy(&las->vMajor,&pubHead[offset],1);
  offset=25;  /*version minor*/
  memcpy(&las->vMinor,&pubHead[offset],1);
  /*fprintf(stdout,"Version %d.%d\n",las->vMajor,las->vMinor);*/
  offset=94;  /*header size*/
  memcpy(&las->headSize,&pubHead[offset],2);
  TIDY(pubHead);

  if((las->vMajor!=1)||(las->vMinor>3)){
    fprintf(stderr,"Version too new for this program\n");
    fprintf(stderr,"Version %d.%d\n",las->vMajor,las->vMinor);
    exit(1);
  }

  /*rewind*/
  if(fseek(las->ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*set header length depending on version*/
  pubHead=challoc(las->headSize,"pubHead",0);

  if(fread(&(pubHead[0]),sizeof(char),las->headSize,las->ipoo)!=las->headSize){
    fprintf(stderr,"error reading data from %s\n",namen);
    exit(1);
  }

  /*copy memory bits over*/
  offset=90;   /*day of year*/
  memcpy(&las->doy,&pubHead[offset],2);
  offset=92;  /*year*/
  memcpy(&las->year,&pubHead[offset],2);
  offset=94;  /*header size*/
  memcpy(&las->headSize,&pubHead[offset],2);
  offset=96;  /*offset to point data*/
  memcpy(&las->offsetToP,&pubHead[offset],4);
  offset=100;  /*number of variable records*/
  memcpy(&las->nVarRec,&pubHead[offset],4);
  offset=104;  /*point data format ID*/
  memcpy(&las->pointFormat,&pubHead[offset],1);
  offset=105;  /*point data format ID*/
  memcpy(&las->pRecLen,&pubHead[offset],1);
  offset=107;  /*number of point records*/
  memcpy(&las->nPoints,&pubHead[offset],4);
  offset=111;  /*number of points by return*/
  for(i=0;i<7;i++)memcpy(&las->nPbyRet[i],&pubHead[offset+4*i],4);

  /*we have lost 8 bytes somewhere here*/
  offset=139-8;
  for(i=0;i<3;i++){  /*point scaling factors*/
    memcpy(&las->posScale[i],&pubHead[offset],8);
    offset+=8;
  }/*read scaling factors*/
  offset=163-8;
  for(i=0;i<3;i++){  /*point offsets*/
    memcpy(&las->posOffset[i],&pubHead[offset+i*8],8);
  }/*read offsets*/

  offset=187-8;
  for(i=0;i<3;i++){  /*point bounds*/
    memcpy(&las->maxB[i],&pubHead[offset+i*2*8],8);
    memcpy(&las->minB[i],&pubHead[offset+i*2*8+8],8);
    //fprintf(stdout,"Bounds %d %f %f\n",i,las->minB[i],las->maxB[i]);
  }/*read bounds*/

  if((las->vMajor==1)&&(las->vMinor==3)){
    offset=235-8;   /*Waveform packet start*/
    memcpy(&(las->waveStart),&(pubHead[offset]),sizeof(uint64_t));
  }

  TIDY(pubHead);

  /*make sure we don't open too many files*/
  if(nFilesOpen>=MAX_OPEN){
    if(las->ipoo){
      fclose(las->ipoo);
      las->ipoo=NULL;
    }
    nFilesOpen--;
  }

  /*set up point reading buffer size*/
  if(pBuffSize>0)las->maxBuffLen=(uint32_t)(pBuffSize/(uint64_t)las->pRecLen);
  else           las->maxBuffLen=1;
  las->maxBuffSize=(uint64_t)las->maxBuffLen*(uint64_t)las->pRecLen;
  las->buffStart=0;    /*buffer start in number of points*/
  las->buffLen=0;      /*buffer length in number of points*/
  las->buffByte=0;     /*buffer length in number of bytes*/

  return(las);
}/*readLasHead*/


/*##############################################*/
/*read geolocation*/

void readLasGeo(lasFile *las)
{
  int i=0,offset=0;
  int headLen=0,varLen=0;
  uint16_t recID=0,reserve=0;
  char **varHead=NULL,*varData=NULL;
  char namen[16];

  /*open file*/
  if((las->ipoo=fopen(las->namen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",las->namen);
    exit(1);
  }
  nFilesOpen++;

  if(fseeko(las->ipoo,(long)las->headSize,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*read the variable header - needed for waveform information*/
  headLen=54;
  las->epsg=0;
  varHead=chChalloc(las->nVarRec,"variable headers",0);
  for(i=0;i<las->nVarRec;i++){
    varHead[i]=challoc(headLen,"variable headers",i+1);
    if(fread(&(varHead[i][0]),sizeof(char),headLen,las->ipoo)!=headLen){
      fprintf(stderr,"Error reading variable header\n");
      exit(1);
    }
    offset=0;
    memcpy(&reserve,&varHead[i][offset],2);
    offset=2;
    memcpy(&(namen[0]),&varHead[i][offset],16);
    offset=18;
    memcpy(&recID,&varHead[i][offset],2);
    offset=20;
    memcpy(&varLen,&varHead[i][offset],2);

    /*read variable part*/
    varData=challoc((int)varLen,"variable data",0);
    if(fread(&(varData[0]),sizeof(char),varLen,las->ipoo)!=varLen){
      fprintf(stderr,"Error reading variable header\n");
      exit(1);
    }
    if(!strncasecmp(namen,"LASF_Projection",16)){  /*geo projection*/
      if(recID==34735){
fprintf(stdout,"vrLen %d\n",(int)varLen);
        offset=102;
        memcpy(&las->epsg,&varData[offset],2);
        fprintf(stdout,"EPSG %d\n",las->epsg);
      }
    }else{/*geo projection*/
      fprintf(stdout,"%s\n",namen);
    }
    TIDY(varData);
  }/*variable header loop*/

  if(las->epsg==0){
    fprintf(stderr,"No geolocation information. Setting default\n");
    las->epsg=32619;
  }
  TTIDY((void **)varHead,las->nVarRec);

  /*close of too many files*/
  if(nFilesOpen>=MAX_OPEN){
    if(las->ipoo){
      fclose(las->ipoo);
      las->ipoo=NULL;
    }
    nFilesOpen--;
  }
  return;
}/*readLasGeo*/


/*##############################################*/
/*read a single point*/

void readLasPoint(lasFile *las,uint32_t j)
{
  uint64_t offset=0;
  uint64_t offTo=0;


  /*if not already, open file*/
  if(las->ipoo==NULL){
    if((las->ipoo=fopen(las->namen,"rb"))==NULL){
      fprintf(stderr,"Error opening input file %s\n",las->namen);
      exit(1);
    }
    nFilesOpen++;
  }


  /*do we need to read new buffer*/
  if((las->pointBuff==NULL)||(j>=(las->buffStart+las->buffLen))||(j<las->buffStart)){
    TIDY(las->pointBuff);
    las->buffStart=j;
    if((j+(uint32_t)las->maxBuffLen)<las->nPoints)las->buffByte=las->maxBuffSize;
    else                                          las->buffByte=(uint64_t)(las->nPoints%las->maxBuffLen)*(uint64_t)las->pRecLen;
    las->buffLen=(uint32_t)(las->buffByte/(uint64_t)las->pRecLen);

    if(!(las->pointBuff=(char *)calloc(las->buffByte,sizeof(char)))){
      fprintf(stderr,"error point buffer allocation.\n");
      exit(1);
    } 

    if(fseeko(las->ipoo,(long)las->offsetToP+(long)((uint64_t)las->buffStart*(uint64_t)las->pRecLen),SEEK_SET)){
      fprintf(stderr,"fseek error\n");
      exit(1);
    }
    if(fread(&(las->pointBuff[0]),sizeof(char),las->buffByte,las->ipoo)!=las->buffByte){
      fprintf(stderr,"Error reading point data, size %d\n",(int)las->buffByte);
      exit(1);
    }
  }/*read buffer*/


  if(nFilesOpen>=MAX_OPEN){
    if(las->ipoo){
      fclose(las->ipoo);
      las->ipoo=NULL;
    }
    nFilesOpen--;
  }

  offTo=(uint64_t)(j-las->buffStart)*(uint64_t)las->pRecLen;

  /*point format 3 and 4*/
  offset=0+offTo;
  memcpy(&las->x,&las->pointBuff[offset],4);
  offset=4+offTo;
  memcpy(&las->y,&las->pointBuff[offset],4);
  offset=8+offTo;
  memcpy(&las->z,&las->pointBuff[offset],4);
  offset=12+offTo;
  memcpy(&las->refl,&las->pointBuff[offset],2);
  offset=14+offTo;
  memcpy(&las->field,&las->pointBuff[offset],1);
  offset=15+offTo;
  memcpy(&las->classif,&las->pointBuff[offset],1);
  offset=16+offTo;
  memcpy(&las->scanAng,&las->pointBuff[offset],1);
  offset=18+offTo;
  memcpy(&las->psID,&las->pointBuff[offset],2);
  if(las->pointFormat==4){   /*full waveform data*/
    offset=28+offTo;
    memcpy(&las->packetDes,&las->pointBuff[offset],1);
    offset=29+offTo;
    memcpy(&las->waveMap,&las->pointBuff[offset],8);
    offset=37+offTo;
    memcpy(&las->waveLen,&las->pointBuff[offset],4);
    offset=41+offTo;
    memcpy(&las->time,&las->pointBuff[offset],4);
    offset=45+offTo;
    memcpy(&las->grad[0],&las->pointBuff[offset],4);
    offset=49+offTo;
    memcpy(&las->grad[1],&las->pointBuff[offset],4);
    offset=53+offTo;
    memcpy(&las->grad[2],&las->pointBuff[offset],4);
  }

  if(las->pointFormat==3){  /*there is RGB*/
    offset=28+offTo;
    memcpy(&las->RGB[0],&las->pointBuff[offset],3*2);
  }/*there is RGB*/

  return;
}/*readLasPoint*/


/*#########################################################################*/
/*read a waveform*/

unsigned char *readLasWave(uint64_t waveMap,int32_t waveLen,FILE *ipoo,uint64_t waveStart)
{
  unsigned char *wave=NULL;

  wave=uchalloc(waveLen,"waveform",(int)waveMap);
  if(fseeko(ipoo,(off_t)((uint64_t)waveMap+(uint64_t)waveStart),SEEK_SET)){
    printf("Error seeking through las file\n");
    exit(1);
  }
  if(fread(&(wave[0]),sizeof(unsigned char),waveLen,ipoo)!=waveLen){
    fprintf(stderr,"Error reading waveform at %ld for %ld\n",(long int)((uint64_t)waveMap+(uint64_t)waveStart),(long int)waveLen);
    exit(1);
  }
  return(wave);
}/*readWaveform*/


/*##############################################*/
/*read input file list*/

char **readInList(int *nFiles,char *inList)
{
  int i=0;
  char line[200];
  char **namen=NULL;
  FILE *ipoo=NULL;

  if((ipoo=fopen(inList,"r"))==NULL){
    fprintf(stderr,"Error opening input file list %s\n",inList);
    exit(1);
  }
  i=0;   /*count up files*/
  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1))i++;
  }/*count up files*/

  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }
  *nFiles=i;
  namen=chChalloc(*nFiles,"file names",0);
  i=0;
  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      namen[i]=challoc(strlen(line)+1,"file names",i+1);
      sscanf(line,"%s",namen[i]);
      i++;
    }
  }

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(namen);
}/*readInList*/



/*######################################*/
/*tidy las file*/

lasFile *tidyLasFile(lasFile *las)
{
  if(las){
    if(las->ipoo){
      fclose(las->ipoo);
      las->ipoo=NULL;
    }
    TIDY(las->wave);
    TIDY(las->pointBuff);
    TIDY(las);
  }

  return(las);
}/*tidyLasFile*/


/*#########################################################################*/
/*set coordinates*/

void setCoords(double *lon,double *lat,double *height,lasFile *lasIn)
{
  *lon=(double)lasIn->x*lasIn->posScale[0]+lasIn->posOffset[0];
  *lat=(double)lasIn->y*lasIn->posScale[1]+lasIn->posOffset[1];
  *height=(double)lasIn->z*lasIn->posScale[2]+lasIn->posOffset[2];

  return;
}/*setCoords*/


/*#########################################################################*/
/*determine coordinate of a bin*/

void binPosition(double *x,double *y,double *z,int bin,double xCent,double yCent,double zCent,float time,float *grad)
{
  double r=0;     /*distance along beam froma anchor point*/

  r=(double)bin*1000.0-(double)time;

  *x=xCent+(double)grad[0]*r;
  *y=yCent+(double)grad[1]*r;
  *z=zCent+(double)grad[2]*r;

  return;
}/*binPosition*/


/*############################################*/
/*check there is a waveform and just one per beam*/

char checkOneWave(lasFile *lasIn)
{
  if((lasIn->packetDes)&&(lasIn->field.nRet==lasIn->field.retNumb)&&(lasIn->waveLen>0))return(1);
  else                                                                                 return(0);
}/*checkOneWave*/


/*############################################*/
/*check file bounds*/

char checkFileBounds(lasFile *lasIn,double minX,double maxX,double minY,double maxY)
{
  if((lasIn->minB[0]<=maxX)&&(lasIn->minB[1]<=maxY)&&\
     (lasIn->maxB[0]>=minX)&&(lasIn->maxB[1]>=minY))return(1);
  else                                              return(0);
}/*checkFileBounds*/

/*the end*/
/*######################################*/

