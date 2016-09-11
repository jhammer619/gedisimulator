#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "tools.c"
#include "libLasRead.h"
#include "libLasProcess.h"


/*#########################*/
/*# Tests GEDI waveform   #*/ 
/*# generation   2015     #*/
/*# svenhancock@gmail.com #*/
/*#########################*/



/*#######################################*/
/*# Copyright 2015-2016, Steven Hancock #*/
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



uint64_t globPoints;

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
/*control structure*/

typedef struct{
  char **inList;
  int nFiles;
  char outRoot[200];
  char parNamen[200];     /*parameter filename*/
  FILE *parFile;     /*parameter file*/

  /*options*/
  uint64_t pBuffSize; /*point buffer rading size in bytes*/
  double footRes;     /*footprint grid resolution*/
  char writeWave;     /*write waveforms switch*/
  char writePar;      /*write parameter switch*/
  char roundOrigin;   /*round origin to ease file labelling*/
  char checkCover;    /*check footprint coverage*/
  float rhoRatio;     /*canopy reflectance over htound reflectance*/

  /*GEDI fotprint parameters*/
  char sideLobe;     /*side lobe switch*/
  float lobeAng;     /*lobe major axis, degrees*/
  int nLobes;        /*number of side lobes*/
  lobeStruct *lobe;  /*lobe structure*/
  double minX;       /*minimum latitude of interest offset*/
  double maxX;       /*maximum latitude of interest offset*/
  double minY;       /*minimum longitude of interest offset*/
  double maxY;       /*maximum longitude of interest offset*/
  
  float res;       /*range resolution*/
  float pFWHM;     /*pulse width in ns*/
  float pSigma;    /*pulse width in metres*/
  float fWidth;    /*footprint width*/
  float fSigma;    /*footprint width in sigma*/
  float pRes;      /*pulse range resolution*/
  pulseStruct *pulse;  /*pulse structure*/

  /*tolerances*/
  float iThresh;   /*intensity threshold*/

  /*grid to normalise sampling density*/
  int gX;
  int gY;
  float gridRes;
}control;


/*####################################*/
/*data structure*/

typedef struct{
  uint32_t nPoints;
  double *x;
  double *y;
  double *z;
  unsigned char *class;     /*point classification*/
  int *refl;
  char *retNumb;            /*discrete return number*/
}dataStruct;


/*####################################*/
/*grid structure*/

typedef struct{
  /*global grid*/
  uint64_t nX;
  uint64_t nY;
  uint64_t nDone;
  uint64_t *done;/*global list of done*/
  double minX;   /*corner for this section of data*/
  double minY;   /*corner for this section of data*/
  double globX;  /*global corner*/
  double globY;  /*global corner*/
  double res;   /*grid resolution*/
  float rRes;   /*range resolution*/

  /*file map*/
  int nWithin;   /*number of files intersecting*/
  int *map;      /*file map*/

  /*within file*/
  int nBins;        /*number of bins per waveform*/
  double minZ;      /*min elevation*/
  double maxZ;      /*max elevation*/
  int wX;           /*number of waves in x*/
  int wY;           /*number of waves in y*/
  int **nGrid;      /*grid cound to even weighting*/
  double **g0;      /*grid origins*/
  double **coord;   /*centres*/
  char *doIt;       /*local list of what needs doing*/

  /*data*/
  dataStruct *data;
  int *nFilesIn;    /*number of files contributing*/
}gridStruct;


/*####################################*/
/*waveform structure*/

typedef struct{
  float **wave;   /*waveforms*/
  float **canopy; /*canopy waveform*/
  float **ground; /*ground waveform*/
  double *z;     /*elevation*/
  double minZ;   /*elevation bounds*/
  double maxZ;   /*elevation bounds*/
  int nBins;     /*number of wave bins*/
  int nWaves;    /*number of different ways*/
  float cover;
  float slope;
}waveStruct;


/*####################################*/
/*main*/

int main(int argc,char **argv)
{
  int i=0,iX=0,iY=0,place=0;
  control *dimage=NULL;
  control *readCommands(int,char **);
  lasFile **las=NULL;
  gridStruct *grid=NULL;
  gridStruct *setUpGrid(lasFile **,control *);
  waveStruct *wave=NULL;
  waveStruct *makeGediWaves(control *,dataStruct *,gridStruct *,int,int);
  void writeGEDIwave(control *,waveStruct *);
  void setGediFootprint(control *);
  void setGediPulse(control *);
  void tidySMoothPulse();
  void checkThisFile(lasFile *,control *,int);
  void readGediGrid(control *,lasFile **,int,gridStruct *);
  void findPars(waveStruct *,control *);
  void writePars(waveStruct *,control *,double *);
  void writeWave(waveStruct *,control *,double *);

 
  /*read command line*/
  dimage=readCommands(argc,argv);

  /*open output if needed*/
  if(dimage->writePar){
    sprintf(dimage->parNamen,"%s.par",dimage->outRoot);
    if((dimage->parFile=fopen(dimage->parNamen,"w"))==NULL){
      fprintf(stderr,"Error opening output file %s\n",dimage->parNamen);
      exit(1);
    }
  }


  if(!(las=(lasFile **)calloc(dimage->nFiles,sizeof(lasFile *)))){
    fprintf(stderr,"error waveStruct allocation.\n");
    exit(1);
  }
  /*read headers and bounds*/
  for(i=0;i<dimage->nFiles;i++){
    las[i]=readLasHead(dimage->inList[i],dimage->pBuffSize);
    if(las[i]->ipoo){
      fclose(las[i]->ipoo);
      las[i]->ipoo=NULL;
    }
  }


  /*set up grid parameters*/
  grid=setUpGrid(las,dimage);

  /*set up the pulse*/
  setGediPulse(dimage);

  /*allocate space for footprint*/
  setGediFootprint(dimage);


  /*loop over las files*/
  for(i=0;i<dimage->nFiles;i++){
    fprintf(stdout,"File %d of %d\n",i+1,dimage->nFiles);
    fflush(stdout);
    /*read data*/
    readGediGrid(dimage,las,i,grid);

    /*loop over waveforms*/
    for(iX=0;iX<grid->wX;iX++){
      for(iY=0;iY<grid->wY;iY++){
        place=iY*grid->wX+iX;
        if(grid->doIt[place]==0)continue;

        /*make waveform*/
        wave=makeGediWaves(dimage,&(grid->data[place]),grid,iX,iY);

        /*find cover and slope*/
        findPars(wave,dimage);

        /*output*/
        if(dimage->writePar)writePars(wave,dimage,grid->coord[place]);
        if(dimage->writeWave)writeWave(wave,dimage,grid->coord[place]);

        /*tidy data*/
        TIDY(grid->data[place].x);
        TIDY(grid->data[place].y);
        TIDY(grid->data[place].z);
        TIDY(grid->data[place].class);
        TIDY(grid->data[place].refl);
        TIDY(grid->data[place].retNumb);
        TIDY(grid->map);
        grid->nWithin=0;
        /*tidy wave*/
        if(wave){
          TTIDY((void **)wave->wave,wave->nWaves);
          TTIDY((void **)wave->canopy,wave->nWaves);
          TTIDY((void **)wave->ground,wave->nWaves);
          TIDY(wave->z);
          TIDY(wave);
        }
      }
    }/*waveform loop*/

    /*tidy up*/
    TIDY(grid->data);
    TIDY(grid->doIt);
    TTIDY((void **)grid->nGrid,grid->wX*grid->wY);
    TTIDY((void **)grid->coord,grid->wX*grid->wY);
    TTIDY((void **)grid->g0,grid->wX*grid->wY);
  }/*file loop*/


  if(dimage->writePar){
    if(dimage->parFile){
      fclose(dimage->parFile);
      dimage->parFile=NULL;
    }
    fprintf(stdout,"Parameters written to %s\n",dimage->parNamen);
  }


  /*tidy up*/
  if(las){
    for(i=0;i<dimage->nFiles;i++)las[i]=tidyLasFile(las[i]);
    TIDY(las);
  }
  if(grid){
    TIDY(grid->done);
    TIDY(grid);
  }
  if(dimage){
    if(dimage->pulse){
      TIDY(dimage->pulse->y);
      TIDY(dimage->pulse->x);
      TIDY(dimage->pulse);
    }
    TIDY(dimage->lobe);
    TTIDY((void **)dimage->inList,dimage->nFiles);
    TIDY(dimage);
  }
  tidySMoothPulse();
  return(0);
}/*main*/


/*####################################*/
/*write waveform*/

void writeWave(waveStruct *wave,control *dimage,double *coord)
{
  int i=0,j=0;
  char namen[200];
  FILE *opoo=NULL;

  sprintf(namen,"%s.x.%d.y.%d.wave",dimage->outRoot,(int)coord[0],(int)coord[1]);
  if((opoo=fopen(namen,"w"))==NULL){
    fprintf(stderr,"Error opening output file %s\n",namen);
    exit(1);
  }

  for(i=0;i<wave->nBins;i++){
    fprintf(opoo,"%f",wave->z[i]);
    for(j=0;j<wave->nWaves;j++){
      fprintf(opoo," %f %f %f",wave->wave[j][i],wave->canopy[j][i],wave->ground[j][i]);
    }
    fprintf(opoo,"\n");
  }

  if(opoo){
    fclose(opoo);
    opoo=NULL;
  }
  fprintf(stdout,"Wave to %s\n",namen);
  return;
}/*writeWave*/


/*####################################*/
/*write land parameters*/

void writePars(waveStruct *wave,control *dimage,double *coord)
{
  fprintf(dimage->parFile,"%f %f %f %f\n",coord[0],coord[1],wave->cover,wave->slope);
  return;
}/*writePars*/


/*####################################*/
/*find surface parameters*/

void findPars(waveStruct *wave,control *dimage)
{
  int i=0;
  float totGr=0;
  float totCan=0;
  float var=0,meanG=0;
  float dz=0;

  totGr=totCan=0.0;
  meanG=0.0;
  for(i=0;i<wave->nBins;i++){
    totGr+=wave->ground[1][i];
    totCan+=wave->canopy[1][i];
    meanG+=(float)((double)wave->ground[1][i]*wave->z[i]);
  }

  /*slope*/
  if(totGr>0.0){
    meanG/=totGr;
    var=0.0;
    for(i=0;i<wave->nBins;i++){
      dz=((float)wave->z[i]-meanG);
      var+=dz*dz*wave->ground[1][i];
    }
    var/=totGr;

    if(var>(dimage->pSigma*dimage->pSigma)){
      wave->slope=atan2(sqrt(var-dimage->pSigma*dimage->pSigma),dimage->fSigma)*180.0/M_PI;
    }else wave->slope=0.0;
  }else wave->slope=-1.0;

  /*canopy cover*/
  if((totGr+totCan)>0.0)wave->cover=totCan/(totCan+totGr*dimage->rhoRatio);
  else                  wave->cover=-1.0;

  return;
}/*findPars*/


/*####################################*/
/*make a GEDI waveform*/

waveStruct *makeGediWaves(control *dimage,dataStruct *data,gridStruct *grid,int iX,int iY)
{
  int i=0,j=0,n=0;
  int place=0;
  int gX=0,gY=0;
  int bin=0;
  waveStruct *wave=NULL;
  float rScale=0;
  float tot=0;
  double dX=0,dY=0,sep=0;
  void cleanOutliers(waveStruct *,control *);


  if(!(wave=(waveStruct *)calloc(1,sizeof(waveStruct)))){
    fprintf(stderr,"error dataStruct allocation.\n");
    exit(1);
  }
  wave->nBins=grid->nBins;
  wave->nWaves=2;
  wave->canopy=fFalloc(wave->nWaves,"",0);
  wave->ground=fFalloc(wave->nWaves,"",0);
  wave->wave=fFalloc(wave->nWaves,"",0);
  for(i=0;i<wave->nWaves;i++){
    wave->canopy[i]=falloc(wave->nBins,"",0);
    wave->ground[i]=falloc(wave->nBins,"",0);
    wave->wave[i]=falloc(wave->nBins,"",0);
  }
  wave->z=dalloc(wave->nBins,"",0);
  wave->maxZ=grid->maxZ;
  wave->minZ=grid->minZ;

  place=iY*grid->wX+iX;

  for(i=0;i<wave->nBins;i++)wave->z[i]=wave->maxZ-(double)i*(double)dimage->res;

  for(n=0;n<dimage->nLobes;n++){

    for(i=0;i<data->nPoints;i++){
      dX=data->x[i]-(dimage->lobe[n].coord[0]+grid->coord[place][0]);
      dY=data->y[i]-(dimage->lobe[n].coord[1]+grid->coord[place][1]);
      sep=sqrt(dX*dX+dY*dY);

      rScale=(float)gaussian(sep,(double)dimage->lobe[n].fSigma,0.0);
      if(rScale>dimage->iThresh){  /*if bright enough to matter*/
        /*normalise for uneven sampling*/
        gX=(int)((data->x[i]-grid->g0[place][0])/(double)dimage->gridRes);
        gY=(int)((data->y[i]-grid->g0[place][1])/(double)dimage->gridRes);
        if((gX>=0)&&(gX<dimage->gX)&&(gY>=0)&&(gY<dimage->gY)){
          if(grid->nGrid[place][gY*dimage->gX+gX]>0)rScale/=(float)grid->nGrid[place][gY*dimage->gX+gX];
        }

        /*convolve by pulse shape*/
        for(j=0;j<dimage->pulse->nBins;j++){
          bin=(int)((wave->maxZ-data->z[i]+(double)dimage->pulse->x[j])/(double)dimage->res);
//fprintf(stdout,"%f %f %f %d\n",wave->minZ,wave->maxZ,data->z[i],bin);
          if((bin>=0)&&(bin<wave->nBins)){
            if(data->class[i]==2){
              wave->ground[0][bin]+=rScale*dimage->pulse->y[j]*(float)data->refl[i];
              wave->ground[1][bin]+=rScale*dimage->pulse->y[j];
            }else{
              wave->canopy[0][bin]+=rScale*dimage->pulse->y[j]*(float)data->refl[i];
              wave->canopy[1][bin]+=rScale*dimage->pulse->y[j];
            }
            wave->wave[0][bin]+=rScale*dimage->pulse->y[j]*(float)data->refl[i];
            wave->wave[1][bin]+=rScale*dimage->pulse->y[j];
          }/*bound check*/
        }/*pulse loop*/
      }/*intensity check*/
    }/*point loop*/
  }/*lobe loop*/

  /*tidy waveform if needed*/
  cleanOutliers(wave,dimage);

  /*normalise integral*/
  for(j=0;j<wave->nWaves;j++){
    tot=0.0;
    for(i=0;i<wave->nBins;i++)tot+=wave->wave[j][i]*dimage->res;
    if(tot>0.0){
      for(i=0;i<wave->nBins;i++){
        wave->wave[j][i]/=tot;
        wave->canopy[j][i]/=tot;
        wave->ground[j][i]/=tot;
      }
    }
  }

  return(wave);
}/*makeGediWaves*/


/*###################################*/
/*clean outlier points from waveform*/

void cleanOutliers(waveStruct *waves,control *dimage)
{
  int i=0,j=0,gStart=0;
  char pastGround=0;
  float gGap=0;  /*gap in ground return*/
  float maxGap=0;
  float maxGround=0,gThresh=0;
  float max=0,thresh=0;

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
/*set footprint grid for this file*/

void readGediGrid(control *dimage,lasFile **las,int i,gridStruct *grid)
{
  int j=0,l=0,m=0;
  int iX=0,iY=0,place=0;
  int gX=0,gY=0;
  int **doneFiles=NULL,*usedFiles=NULL;
  uint32_t k=0;
  float rScale=0;
  double x=0,y=0,z=0;
  double dx=0,dy=0;
  double xOff=0,yOff=0;
  double sepSq=0;
  double xOrigin=0,yOrigin=0;
  void recordPoint(dataStruct *,double,double,double,lasFile *);
  char checkFootCovered(int,int,int *,float,float);
  char checkDone(gridStruct *,int,int,double,double);
  char newFile=0;

  /*determine how many files overlap. Need to include a buffer*/
  TIDY(grid->map);
  grid->nWithin=0;
  grid->map=markInt(grid->nWithin,grid->map,i);
  grid->nWithin=1;
  for(j=0;j<dimage->nFiles;j++){
    if(i==j)continue;
    if((las[j]->minB[0]<=(las[i]->maxB[0]+dimage->maxX))&&(las[j]->minB[1]<=(las[i]->maxB[1]+dimage->maxY))&&\
       (las[j]->maxB[0]>=(las[i]->minB[0]-dimage->minX))&&(las[j]->maxB[1]>=(las[i]->minB[1]-dimage->minY))){
      grid->map=markInt(grid->nWithin,grid->map,j);
      grid->nWithin++;
    }
  }/*determine how many files overlap*/

  /*determine bounds rounding to match global*/
  if(las[i]->minB[0]>grid->globX){
    xOrigin=ceil((las[i]->minB[0]-grid->globX)/grid->res)*grid->res+grid->globX;
  }else xOrigin=grid->globX;
  if(las[i]->minB[1]>grid->globY){
    yOrigin=ceil((las[i]->minB[1]-grid->globY)/grid->res)*grid->res+grid->globY;
  }else yOrigin=grid->globY;


  /*determine elevation bounds*/
  grid->maxZ=-1000000000.0;
  grid->minZ=1000000000.0;
  for(j=0;j<grid->nWithin;j++){
    if(las[j]->minB[2]<grid->minZ)grid->minZ=las[j]->minB[2];
    if(las[j]->maxB[2]>grid->maxZ)grid->maxZ=las[j]->maxB[2];
  }/*determine bounds*/
  /*add buffers*/
  grid->minZ-=50.0;
  grid->maxZ+=50.0;

  grid->nBins=(int)((grid->maxZ-grid->minZ)/grid->rRes);

  /*allocate waveforms and gridding structures*/
  grid->nBins=(int)((grid->maxZ-grid->minZ)/(double)grid->rRes);
  grid->wX=(int)(ceil((las[i]->maxB[0]-xOrigin)/grid->res));
  grid->wY=(int)(ceil((las[i]->maxB[1]-yOrigin)/grid->res));
  grid->nGrid=iIalloc(grid->wX*grid->wY,"grid",0);
  grid->g0=dDalloc(grid->wX*grid->wY,"grid origin",0);
  grid->coord=dDalloc(grid->wX*grid->wY,"footprint centre",0);
  grid->doIt=challoc(grid->wX*grid->wY,"do it map",0);
  for(j=0;j<grid->wX;j++){
    xOff=(double)j*grid->res+xOrigin;
    for(k=0;k<grid->wY;k++){
      yOff=(double)k*grid->res+yOrigin;

      place=k*grid->wX+j;
      grid->doIt[place]=1;

      /*grid for normalising sampling*/
      grid->g0[place]=dalloc(2,"grid origin",place+1);
      grid->g0[place][0]=xOff+dimage->minX;
      grid->g0[place][1]=yOff+dimage->minY;
      grid->coord[place]=dalloc(2,"footprint centre",place+1);
      grid->coord[place][0]=xOff;
      grid->coord[place][1]=yOff;
      grid->nGrid[place]=ialloc(dimage->gX*dimage->gY,"nGrid",0);
      for(l=dimage->gX*dimage->gY-1;l>=0;l--)grid->nGrid[place][l]=0;
    }
  }

  /*allocate maps*/
  usedFiles=ialloc(grid->wX*grid->wY,"usedFile map",0);
  doneFiles=iIalloc(grid->wX*grid->wY,"done file map",0);
  

  if(!(grid->data=(dataStruct *)calloc(grid->wX*grid->wY,sizeof(dataStruct)))){
    fprintf(stderr,"error dataStruct allocation.\n");
    exit(1);
  }
  grid->nFilesIn=ialloc(grid->wX*grid->wY,"nFilesIn",0);
  for(l=grid->wX*grid->wY-1;l>=0;l--){
    grid->data[l].nPoints=0;
    grid->nFilesIn[l]=0;
    usedFiles[l]=0;
    doneFiles[l]=NULL;
  }

  /*loop files and points and assign to waveforms*/
  for(j=0;j<grid->nWithin;j++){
    for(k=0;k<las[grid->map[j]]->nPoints;k++){
      readLasPoint(las[grid->map[j]],k);
      setCoords(&x,&y,&z,las[grid->map[j]]);

      /*see which footprints it is within*/
      for(iX=0;iX<grid->wX;iX++){
        xOff=(double)iX*grid->res+xOrigin;
        for(iY=0;iY<grid->wY;iY++){
          yOff=(double)iY*grid->res+yOrigin;
          place=iY*grid->wX+iX;
    
          for(l=0;l<dimage->nLobes;l++){
            dx=(xOff+dimage->lobe[l].coord[0])-x;
            dy=(yOff+dimage->lobe[l].coord[1])-y;

            sepSq=dx*dx+dy*dy;
            if(sepSq<=dimage->lobe[l].maxSepSq){
              rScale=(float)gaussian(sqrt(sepSq),(double)dimage->lobe[l].fSigma,0.0);

              if(rScale>dimage->iThresh){  /*if bright enough to matter*/
                /*record the point*/
                recordPoint(&(grid->data[place]),x,y,z,las[grid->map[j]]);

                /*mark up the grid*/
                if(las[grid->map[j]]->field.retNumb==las[grid->map[j]]->field.nRet){  /*only once per beam*/
                  gX=(int)((x-grid->g0[place][0])/(double)dimage->gridRes);
                  gY=(int)((y-grid->g0[place][1])/(double)dimage->gridRes);
                  if((gX>=0)&&(gX<dimage->gX)&&(gY>=0)&&(gY<dimage->gY)){
                    grid->nGrid[place][gY*dimage->gX+gX]++;
                  }
                }

                /*see if this is a new file*/
                newFile=1;
                for(m=0;m<usedFiles[place];m++){
                  if(doneFiles[place][m]==grid->map[j]){
                    newFile=0;
                    break;
                  }
                }
                if(newFile){
                  doneFiles[place]=markInt(usedFiles[place],doneFiles[place],grid->map[j]);
                  usedFiles[place]++;
                  grid->nFilesIn[place]++;
                }
              }/*intensity test*/
            }/*within footprint check*/
          }/*lobe loop*/
        }/*y grid loop*/
      }/*x grid loop*/
    }/*point loop*/
  }/*file loop*/

  /*tidy maps*/
  TIDY(usedFiles);
  TTIDY((void **)doneFiles,grid->wX*grid->wY);


  /*don't do waveforms with insufficient data*/
  for(iX=0;iX<grid->wX;iX++){
    xOff=(double)iX*grid->res+xOrigin;
    for(iY=0;iY<grid->wY;iY++){
      yOff=(double)iY*grid->res+yOrigin;
      place=iY*grid->wX+iX;
      if(grid->data[place].nPoints==0)grid->doIt[place]=0;
      if(grid->nFilesIn[place]>1)grid->doIt[place]=checkDone(grid,iX,iY,xOrigin,yOrigin);
      else if(dimage->checkCover&&(grid->doIt[place]==1)){
       if(checkFootCovered(dimage->gX,dimage->gY,grid->nGrid[place],dimage->fSigma,dimage->gridRes))grid->doIt[place]=0;
      }
    }
  }

  /*tidy deleted footprints*/
  for(place=grid->wX*grid->wY-1;place>=0;place--){
    if(grid->doIt[place]==0){
      TIDY(grid->data[place].x);
      TIDY(grid->data[place].y);
      TIDY(grid->data[place].z);
      TIDY(grid->data[place].class);
      TIDY(grid->data[place].refl);
      TIDY(grid->data[place].retNumb);
    }
  }

  return;
}/*readGediGrid*/


/*####################################*/
/*check footprint is covered by ALS*/

char checkFootCovered(int gX,int gY,int *nGrid,float sigma,float gridRes)
{
  int i=0,j=0,nWithin=0;
  int thresh=0,nMissed=0;
  float dX=0,dY=0,sepSq=0;
  float radSq=0;

  radSq=sigma*sigma;

  for(i=0;i<gX;i++){
    dX=(float)(i-gX/2)*gridRes;
    for(j=0;j<gY;j++){
      dY=(float)(j-gY/2)*gridRes;
      sepSq=dX*dX+dY*dY;

      if(sepSq<radSq){
        if(nGrid[j*gX+i]==0)nMissed++;
        nWithin++;
      }
    }/*y loop*/
  }/*x loop*/

  thresh=(int)((float)nWithin*2.0/3.0);
  if(nMissed>thresh)return(1);
  else              return(0);
}/*checkFootCovered*/


/*####################################*/
/*check whether it's been done*/

char checkDone(gridStruct *grid,int iX,int iY,double xOrigin,double yOrigin)
{
  int i=0;
  uint64_t xOff=0,yOff=0;
  uint64_t place;
  char doIt=0;

  xOff=(uint64_t)(xOrigin/grid->res+0.5);
  yOff=(uint64_t)(yOrigin/grid->res+0.5);

  place=((uint64_t)iY+yOff)*grid->nX+(uint64_t)iX+xOff;

  doIt=1;
  for(i=0;i<grid->nDone;i++){
    if(grid->done[i]==place){
      doIt=0;
      break;
    }
  }

  if(doIt){
    if(grid->nDone>0){
      if(!(grid->done=(uint64_t *)realloc(grid->done,(grid->nDone+1)*sizeof(uint64_t)))){
        fprintf(stderr,"Balls\n");
        exit(1);
      }
    }else{
      if(!(grid->done=(uint64_t *)calloc(1,sizeof(uint64_t)))){
        fprintf(stderr,"error in array allocation\n");
        exit(1);
      }
    }
    grid->done[grid->nDone]=place;  
    grid->nDone++;
  }

  return(doIt);
}/*checkDone*/


/*####################################*/
/*record a las point*/

void recordPoint(dataStruct *data,double x,double y,double z,lasFile *las)
{
  data->x=markDo((int)data->nPoints,data->x,x);
  data->y=markDo((int)data->nPoints,data->y,y);
  data->z=markDo((int)data->nPoints,data->z,z);
  data->refl=markInt((int)data->nPoints,data->refl,las->refl);
  data->retNumb=markChar((int)data->nPoints,data->retNumb,(char)las->field.retNumb);
  data->class=markUchar((int)data->nPoints,data->class,(unsigned char)las->classif);
//fprintf(stdout,"Recorded %d %llu\n",data->nPoints,globPoints++);
  data->nPoints++;
  return;
}/*recordPoint*/ 


/*####################################*/
/*set up grid*/

gridStruct *setUpGrid(lasFile **las,control *dimage)
{
  int i=0;
  double maxX=0,maxY=0;
  gridStruct *grid=NULL;

  if(!(grid=(gridStruct *)calloc(1,sizeof(gridStruct)))){
    fprintf(stderr,"error waveStruct allocation.\n");
    exit(1);
  }

  /*global scene bounds*/
  grid->globX=grid->globY=1000000000000.0;
  maxX=maxY=-10000000000.0;
  for(i=0;i<dimage->nFiles;i++){
    if(las[i]->minB[0]<grid->globX)grid->globX=las[i]->minB[0];
    if(las[i]->minB[1]<grid->globY)grid->globY=las[i]->minB[1];
    if(las[i]->maxB[0]>maxX)maxX=las[i]->maxB[0];
    if(las[i]->maxB[1]>maxY)maxY=las[i]->maxB[1];
  }

  if(dimage->roundOrigin){
    grid->globX=ceil(grid->globX);
    grid->globY=ceil(grid->globY);
  }

  grid->done=NULL;
  grid->nDone=0;
  grid->res=dimage->footRes;
  grid->rRes=0.15;

  grid->nX=(uint64_t)((maxX-grid->minX)/grid->res);
  grid->nY=(uint64_t)((maxX-grid->minX)/grid->res);

  return(grid);
}/*setUpGrid*/


/*####################################*/
/*set GEDI footprint*/

void setGediFootprint(control *dimage)
{
  int i=0;
  float maxSep=0;
  float determineGaussSep(float,float);
  float totE=0;
  float az=0;
  double tX=0,tY=0;


  /*footprint width*/
  if(dimage->fSigma<0.0)dimage->fSigma=dimage->fWidth;
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
  dimage->lobe[i].coord[0]=0.0;
  dimage->lobe[i].coord[1]=0.0;
  dimage->lobe[i].fSigma=dimage->fSigma;
  maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*1.0);
  dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);
  if(dimage->sideLobe==0)dimage->lobe[i].E=1.0;
  else{  /*include side lobes*/
    totE=1.0+0.0599+0.0731+0.0317+0.0319+0.0167+0.0163;
    i=0;
    dimage->lobe[i].E=1.0/totE;

    /*first southern lobe*/
    i=1;
    dimage->lobe[i].E=0.0731/totE;
    dimage->lobe[i].coord[0]=0.0-20.0*sin(az);
    dimage->lobe[i].coord[1]=0.0-20.0*cos(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);

    /*first nothern lobe*/
    i=2;
    dimage->lobe[i].E=0.0599/totE;
    dimage->lobe[i].coord[0]=0.0+20.0*sin(az);
    dimage->lobe[i].coord[1]=0.0+20.0*cos(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);

    /*western lobe*/
    i=3;
    dimage->lobe[i].E=0.0319/totE;
    dimage->lobe[i].coord[0]=0.0-20.0*cos(az);
    dimage->lobe[i].coord[1]=0.0-20.0*sin(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);

    /*eastern lobe*/
    i=4;
    dimage->lobe[i].E=0.0317/totE;
    dimage->lobe[i].coord[0]=0.0+20.0*cos(az);
    dimage->lobe[i].coord[1]=0.0+20.0*sin(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);

    /*second southern lobe*/
    i=5;
    dimage->lobe[i].E=0.0167/totE;
    dimage->lobe[i].coord[0]=0.0-30.0*sin(az);
    dimage->lobe[i].coord[1]=0.0-30.0*cos(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);

    /*second northern lobe*/
    i=6;
    dimage->lobe[i].E=0.0163/totE;
    dimage->lobe[i].coord[0]=0.0+30.0*cos(az);
    dimage->lobe[i].coord[1]=0.0+30.0*sin(az);
    dimage->lobe[i].fSigma=dimage->fSigma;
    maxSep=determineGaussSep(dimage->lobe[i].fSigma,dimage->iThresh*dimage->lobe[i].E);
    dimage->lobe[i].maxSepSq=(double)(maxSep*maxSep);
  }/*side lobe test*/


  /*determine min and max bounds for a single footprint*/
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
  dimage->gridRes=1.5;
  dimage->gX=(int)((float)(dimage->maxX-dimage->minX)/dimage->gridRes)+2;
  dimage->gY=(int)((float)(dimage->maxY-dimage->minY)/dimage->gridRes)+2;

  return;
}/*setGediFootprint*/


/*####################################*/
/*maximum significant distance*/

float determineGaussSep(float fSigma,float thresh)
{
  float x=0,y=0;

  x=0.0;
  do{
    y=(float)gaussian((double)x,(double)fSigma,0.0);
    x+=0.2;
  }while(y>=thresh);

  return(x);
}/*determineGaussSep*/


/*####################################*/
/*set GEDI pulse*/

void setGediPulse(control *dimage)
{
  int i=0;
  float fwhm=0;   /*FWHM in metres*/
  float x=0,y=0;
  float max=0,tot=0;

  if(!(dimage->pulse=(pulseStruct *)calloc(1,sizeof(pulseStruct)))){
    fprintf(stderr,"error pulseStruct allocation.\n");
    exit(1);
  }

  /*pulse length*/
  /*calculate sigma from FWHM*/
  if(dimage->pSigma<0.0){  /*GEDI unless specificed*/
    fwhm=dimage->pFWHM*0.15;  /*time if for two way*/
    dimage->pSigma=fwhm/2.355;
  }

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

  return;
}/*setGediPulse*/


/*####################################*/
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
  strcpy(dimage->outRoot,"teast");
  dimage->pFWHM=12.0;   /*12 ns FWHM*/
  dimage->fWidth=7.5;  /*22m diameter contains 86% of energy*/
  dimage->res=0.15;
  dimage->pRes=0.01;
  dimage->footRes=50.0;
  dimage->writeWave=0;
  dimage->writePar=1;
  dimage->roundOrigin=1;
  dimage->checkCover=1;
  dimage->rhoRatio=0.57/0.4;

  /*switches*/
  dimage->sideLobe=0;   /*no side lobes*/
  dimage->lobeAng=0.0;
  dimage->pBuffSize=(uint64_t)200000000;

  dimage->iThresh=0.0006;
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
        dimage->inList[0]=challoc(strlen(argv[++i])+1,"input name list",0);
        strcpy(dimage->inList[0],argv[i]);
      }else if(!strncasecmp(argv[i],"-outRoot",8)){
        checkArguments(1,i,argc,"-ouRoot");
        strcpy(dimage->outRoot,argv[++i]);
      }else if(!strncasecmp(argv[i],"-gRes",5)){
        checkArguments(1,i,argc,"-gRes");
        dimage->footRes=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-writeWaves",11)){
        dimage->writeWave=1;
      }else if(!strncasecmp(argv[i],"-inList",7)){
        checkArguments(1,i,argc,"-inList");
        TTIDY((void **)dimage->inList,dimage->nFiles);
        dimage->inList=readInList(&dimage->nFiles,argv[++i]);
      }else if(!strncasecmp(argv[i],"-LVIS",5)){
        dimage->pSigma=0.6893/2.0;  /*two way trip*/
        dimage->fSigma=6.25;
      }else if(!strncasecmp(argv[i],"-pSigma",7)){
        checkArguments(1,i,argc,"-pSigma");
        dimage->pSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-fSigma",7)){
        checkArguments(1,i,argc,"-fSigma");
        dimage->fSigma=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-sideLobe",9)){
        dimage->sideLobe=1;
      }else if(!strncasecmp(argv[i],"-lobeAng",8)){
        checkArguments(1,i,argc,"-lobeAng");
        dimage->lobeAng=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-pBuff",6)){
        checkArguments(1,i,argc,"-pBuff");
        dimage->pBuffSize=(uint64_t)(atof(argv[++i])*1000000000.0);
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n#####\nProgram to create GEDI waveforms from ALS las files\n#####\n\n-input name;     lasfile input filename\n-outRoot root;   output filename root\n-inList list;    input file list for multiple files\n-gRes res;       grid resolution\n-writeWave;      write waveforms\n-LVIS;           use LVIS pulse length, sigma=6.25m\n-pSigma sig;     set pulse width\n-fSigma sig;     set footprint width\n-sideLobe;       use side lobes\n-lobeAng ang;    lobe axis azimuth\n-pBuff s;        point reading buffer size in Gbytes\n\nQuestions to svenhancock@gmail.com\n\n");
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

