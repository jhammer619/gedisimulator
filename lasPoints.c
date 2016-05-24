#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"


/*########################*/
/*# Outputs point clouds #*/
/*# from las files       #*/
/*########################*/




/*####################################*/
/*control structure*/

typedef struct{
  char **inList;
  int nFiles;
  char outRoot[200];
  char canNamen[200];  /*canopy filename*/
  char grNamen[200];   /*ground filename*/
  float thresh;        /*maximum intensity of area of interest*/
  uint64_t pBuffSize;  /*point buffer rading size in bytes*/

  /*output files*/
  FILE *opoo;  /*non-ground*/
  FILE *gPoo;  /*ground points*/

  /*options*/
  char ground;   /*separate ground*/

  /*footprint parameters*/
  float fSigma;      /*footprint sigma*/
  double coord[2];   /*footprint centre*/
  double maxSepSq;
  double maxSep;
}control;


/*####################################*/
/*main*/
  
int main(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  lasFile *las=NULL;
  lasFile *readLasHead(char *,uint64_t);
  lasFile *tidyLasFile(lasFile *);
  void openOutput(control *);
  void readWritePoints(control *,lasFile *);
  void setArea(control *);


  /*read command line*/
  dimage=readCommands(argc,argv);

  /*determine maximum area of interest*/
  setArea(dimage);

  /*open output*/
  openOutput(dimage);

  /*loop over files*/
  for(i=0;i<dimage->nFiles;i++){
    /*read lasFile*/
    las=readLasHead(dimage->inList[i],dimage->pBuffSize);

    /*read and output data*/
    readWritePoints(dimage,las);

    /*tidy lasFIle*/
    las=tidyLasFile(las);
  }/*file loop*/


  fprintf(stdout,"Points to %s\n",dimage->canNamen);
  if(dimage->ground)fprintf(stdout,"Ground to %s\n",dimage->grNamen);


  /*tidy up*/
  if(dimage){
    if(dimage->opoo){
      fclose(dimage->opoo);
      dimage->opoo=NULL;
    }
    if(dimage->gPoo){
      fclose(dimage->gPoo);
      dimage->gPoo=NULL;
    }
    TTIDY((void **)dimage->inList,dimage->nFiles);
    TIDY(dimage);
  }
  return(0);
}/*main*/


/*####################################*/
/*read and write points to file*/

void readWritePoints(control *dimage,lasFile *las)
{
  uint32_t i=0;
  double x=0,y=0,z=0;
  double dX=0,dY=0,sepSq=0;
  void readLasPoint(lasFile *,uint32_t);
  void setCoords(double *,double *,double *,lasFile *);
  char checkFileBounds(lasFile *,double,double,double,double);


  /*check file bounds*/
  if(checkFileBounds(las,dimage->coord[0]-dimage->maxSep,dimage->coord[0]+dimage->maxSep,\
                         dimage->coord[1]-dimage->maxSep,dimage->coord[1]+dimage->maxSep)){
    for(i=0;i<las->nPoints;i++){
      /*read one point*/
      readLasPoint(las,i);
      setCoords(&x,&y,&z,las);

      dX=x-dimage->coord[0];
      dY=y-dimage->coord[1];
      sepSq=dX*dX+dY*dY;
      if(sepSq<=dimage->maxSepSq){
        if(dimage->ground&&(las->classif==2))fprintf(dimage->gPoo,"%f %f %f %d\n",x,y,z,(int)(las->refl));
        else                                 fprintf(dimage->opoo,"%f %f %f %d\n",x,y,z,(int)(las->refl));
      }/*separation check*/
    }/*point loop*/
  }/*file bounds check*/

  return;
}/*readWritePoints*/


/*####################################*/
/*set area of interest*/

void setArea(control *dimage)
{
  float x=0,y=0,res=0;

  res=0.05;
  x=0.0;
  do{
    y=(float)gaussian((double)x,(double)dimage->fSigma,0.0);
    x+=res;
  }while(y>=dimage->thresh);

  dimage->maxSep=(double)x;
  dimage->maxSepSq=(double)(x*x);

  return;
}/*setArea*/


/*####################################*/
/*open output files*/

void openOutput(control *dimage)
{
  sprintf(dimage->canNamen,"%s.can.pts",dimage->outRoot);
  if((dimage->opoo=fopen(dimage->canNamen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",dimage->canNamen);
    exit(1);
  }

  if(dimage->ground){
    sprintf(dimage->grNamen,"%s.ground.pts",dimage->outRoot);
    if((dimage->gPoo=fopen(dimage->grNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",dimage->grNamen);
      exit(1);
    }
  }

  return;
}/*openOutput*/


/*####################################*/
/*read command line*/

control *readCommands(int argc,char **argv)
{
  int i=0;
  control *dimage=NULL;
  char **readInList(int *,char *);

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }

  /*defaults*/
  dimage->nFiles=1;
  dimage->inList=chChalloc(dimage->nFiles,"inList",0);
  dimage->inList[0]=challoc(200,"inList",0);
  strcpy(&(dimage->inList[0][0]),"/Users/stevenhancock/data/gedi/laselva/ALS/tile197.las");
  strcpy(dimage->outRoot,"gTest");
  dimage->opoo=NULL;
  dimage->gPoo=NULL;

  dimage->fSigma=11.0;  /*GEDI*/
  dimage->coord[0]=826367.0;
  dimage->coord[1]=1152271.0;
  dimage->ground=0;
  dimage->thresh=0.01;  /*1%*/

  dimage->pBuffSize=(uint64_t)200000000;


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
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-outRoot");
        strcpy(dimage->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=readInList(&dimage->nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-coord",6)){
        checkArguments(2,i,argc,"-coord");
        dimage->coord[0]=atof(argv[++i]);
        dimage->coord[1]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-LVIS",5)){
        dimage->fSigma=6.25;  /*LVIS*/
      }else if(!strncasecmp(argv[i],"-ground",7)){
        dimage->ground=1;
      }else if(!strncasecmp(argv[i],"-pBuff",6)){
        checkArguments(1,i,argc,"-pBuff");
        dimage->pBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0);
      }else if(!strncasecmp(argv[i],"-thresh",6)){
        checkArguments(2,i,argc,"-thresh");
        dimage->thresh=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to output ALS points within GEDI footprints\n#####\n\n-input name;     lasfile input filename\n-outRoot name;   output filename\n-inList list;    input file list for multiple files\n-coord lon lat;  footprint coordinate in same system as lasfile\n-LVIS;           use LVIS pulse length, sigma=6.25m\n-ground;         output canopy and ground separately\n-thresh t;       energy threshold to accept points\n-pBuff s;        point reading buffer size in Gbytes\n\nQuestions to svenhancock@gmail.com\n\n");
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

