#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"


/*###########################################*/
/*control structure*/

typedef struct{
  char metricNamen[200];   /*metric file input*/
  char trackNamen[200];    /*GEDI tracks to use*/
  char outNamen[200];      /*output filename*/
  double minSep;           /*minimum separation to accept*/
}control;


/*###########################################*/
/*track structure*/

typedef struct{
  double *x;
  double *y;
  int nTracks;
}trackStruct;


/*###########################################*/
/*main*/

int main(int argc,char **argv)
{
  control *dimage=NULL;
  control *readCommands(int,char **);
  trackStruct *tracks=NULL;
  trackStruct *readTracks(char *);
  void selectMetrics(trackStruct *,char *,char *,double);


  /*read commands*/
  dimage=readCommands(argc,argv);

  /*read track locations*/
  tracks=readTracks(dimage->trackNamen);

  /*read metris and output needed*/
  selectMetrics(tracks,dimage->metricNamen,dimage->outNamen,dimage->minSep);

  /*tidy up*/
  if(tracks){
    TIDY(tracks->x);
    TIDY(tracks->y);
    TIDY(tracks);
  }
  TIDY(dimage);
  return(0);
}/*main*/


/*###########################################*/
/*select and write metrics*/

void selectMetrics(trackStruct *tracks,char *metricNamen,char *outNamen,double minSep)
{
  int j=0,xCol=0,yCol=0,nUse=0;
  int64_t i=0,*useList=NULL;
  double sepSq=0,*minSepSq=NULL;
  double x=0,y=0,dx=0,dy=0;
  double thresh=0;
  char line[20000],*token=NULL;;
  char useIt=0,lastTok[100];
  char writtenHead=0;
  FILE *ipoo=NULL,*opoo=NULL;

  /*open metrics*/
  if((ipoo=fopen(metricNamen,"r"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",metricNamen);
    exit(1);
  }
  thresh=minSep*minSep;

  /*choose metrics*/
  if(!(useList=(int64_t *)calloc(tracks->nTracks,sizeof(int64_t)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  } 

  /*allocate distance array and set to blank*/
  minSepSq=dalloc(tracks->nTracks,"useList",0);
  for(i=0;i<tracks->nTracks;i++){
    minSepSq[i]=10000000.0;
    useList[i]=-1;
  }

  /*search for closest footprints to tracks*/
  i=0;
  while(fgets(line,10000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){  /*read data*/
      x=y=1000000000000.0;  /*silly values*/
      /*read coords*/
      j=1;
      token=strtok(line," ");
      while(token){
        if(j==yCol)y=atof(token);
        else if(j==xCol)x=atof(token);
        else if(j>yCol)break;
        token=strtok(NULL," ");
        j++;
      }

      /*loop over all footprints*/
      for(j=0;j<tracks->nTracks;j++){
        /*determine separation*/
        dx=x-tracks->x[j];
        dy=y-tracks->y[j];
        sepSq=dx*dx+dy*dy;
        /*record closes suitable*/
        if((sepSq<minSepSq[j])&&(sepSq<=thresh)){
          if(useList[j]<=0)nUse++;
          useList[j]=i;
          minSepSq[j]=sepSq;
        }
      }
      i++;
    }else{  /*read header*/
      j=1;
      token=strtok(line," ");
      while(token){
        if(!strncasecmp(token,"lon,",4))xCol=atoi(lastTok);
        else if(!strncasecmp(token,"lat,",4))yCol=atoi(lastTok);
        strcpy(lastTok,token);
        token=strtok(NULL," ");
        j++;
      }
      fprintf(stdout,"Metric coord cols %d %d\n",xCol,yCol);
    }
  }
  TIDY(minSepSq);
  fprintf(stdout,"Read %lld metrics and selected %d\n",i,nUse);

  /*rewind to start of file*/
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  /*open output*/
  if((opoo=fopen(outNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",outNamen);
    exit(1);
  }

  /*write out selected*/
  i=0;
  while(fgets(line,10000,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){  /*read data*/
      useIt=0;
      for(j=0;j<tracks->nTracks;j++){
        if(i==useList[j]){
          useIt=1;
          break;
        }
      }
      if(useIt)fprintf(opoo,"%s",line);
      i++;
    }else{  /*write out header*/
      if(!writtenHead){
        fprintf(opoo,"%s",line);
        writtenHead=1;
      }
    }
  }

  /*tidy up*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Written to %s\n",outNamen);
  TIDY(useList);
  return;
}/*selectMetrics*/


/*###########################################*/
/*read track locations*/

trackStruct *readTracks(char *namen)
{
  int i=0;
  trackStruct *tracks=NULL;
  char line[400],temp1[200],temp2[200];
  FILE *ipoo=NULL;

  /*allocate*/
  if(!(tracks=(trackStruct *)calloc(1,sizeof(trackStruct)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  tracks->x=NULL;
  tracks->y=NULL;
  tracks->nTracks=0;

  /*open data*/
  if((ipoo=fopen(namen,"r"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }

  /*count number of lines*/
  while(fgets(line,400,ipoo)!=NULL)if(strncasecmp(line,"#",1))tracks->nTracks++;
  tracks->x=dalloc(tracks->nTracks,"x tracks",0);
  tracks->y=dalloc(tracks->nTracks,"y tracks",0);

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
        tracks->x[i]=atof(temp1);
        tracks->y[i]=atof(temp2);
        i++;
      }else{
        fprintf(stderr,"Error reading %s\n",namen);
        exit(1);
      }
    }
  }/*data reading loop*/

  /*close file*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  fprintf(stdout,"There will be %d footprints\n",tracks->nTracks);
  return(tracks);
}/*readTracks*/


/*###########################################*/
/*read commands*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;

  /*allocate*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-metric",7)){
        checkArguments(1,i,argc,"-metric");
        strcpy(dimage->metricNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-tracks",7)){
        checkArguments(1,i,argc,"-tracks");
        strcpy(dimage->trackNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-minSep",7)){
        checkArguments(1,i,argc,"-minSep");
        dimage->minSep=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry chooseMetricPrints -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }

  return(dimage);
}/*readCommands*/

/*the end*/
/*###########################################*/
