#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include "libLasProcess.h"


/*########################*/
/*# Functions to process #*/
/*# waveform lidar data  #*/
/*########################*/

#define TOL 0.00001


/*##################################################*/
/*process waveform to denoise and deconvolve*/

float *processWave(unsigned char *wave,int waveLen,denPar *decon,float gbic)
{
  int i=0;
  float *temp=NULL;
  float *mediated=NULL;
  float *preSmoothed=NULL;
  float *denoised=NULL;
  float *smoothed=NULL;
  float *processed=NULL;
  float *denoise(float,float,int,int,float *,float,char);
  float *medianFloat(float *,int,int);
  float *deconvolve(float *,int,float **,int,float,int,double,char);
  float thisTail=0;        /*tail threshold to use here*/
  float *gaussWave=NULL;
  float *fitGaussians(float *,int,denPar *);
  float *hardHitWave(denPar *,int);
  void medNoiseStats(float *,uint32_t,float *,float *,float *,float,float,char);
  void meanNoiseStats(float *,uint32_t,float *,float *,float *,float,float,int);
  char checkHardEnergy(int *,float *,float);
  char hardTarget=0;


  /*convert to a float array for ease*/
  temp=falloc(waveLen,"presmoothed",0);
  for(i=0;i<waveLen;i++)temp[i]=(float)wave[i];


  /*determine noise statistics*/
  if(decon->varNoise){
    if(decon->medStats)medNoiseStats(temp,waveLen,&(decon->meanN),&(decon->thresh),&thisTail,decon->tailThresh,decon->threshScale,decon->bitRate);    
    else               meanNoiseStats(temp,waveLen,&(decon->meanN),&(decon->thresh),&thisTail,decon->tailThresh,decon->threshScale,(int)(decon->statsLen/decon->res));
  }else               thisTail=decon->tailThresh;
  if(thisTail<0)thisTail=decon->thresh;


  /*median filter if needed*/
  if(decon->medLen>0){
    mediated=medianFloat(temp,decon->medLen,waveLen);
    TIDY(temp);
  }else{
    mediated=temp;
    temp=NULL;
  }

  /*smooth before denoising*/
  if(decon->psWidth>0.0){
    preSmoothed=smooth(decon->psWidth,waveLen,mediated,decon->res);
    TIDY(mediated);
  }else{
    preSmoothed=mediated;
    mediated=NULL;
  }

  /*remove background noise*/
  if((decon->meanN>0.0)||(decon->thresh>0.0)){
    denoised=denoise(decon->meanN,decon->thresh,decon->minWidth,waveLen,preSmoothed,thisTail,decon->noiseTrack);
    TIDY(preSmoothed);
  }else{
    denoised=preSmoothed;
    preSmoothed=NULL;
  }

  /*smooth if required. Note that pulse is smoothed in readPulse()*/
  if(decon->sWidth>0.0){
    smoothed=smooth(decon->sWidth,waveLen,denoised,decon->res);
    TIDY(denoised);
  }else{
    smoothed=denoised;
    denoised=NULL;
  }

  /*scale by GBIC*/
  if((gbic>0.0)&&(gbic!=1.0))for(i=0;i<waveLen;i++)smoothed[i]/=gbic;


  /*Gaussian fitting*/
  if(decon->fitGauss||decon->gaussFilt){
    gaussWave=fitGaussians(smoothed,waveLen,decon);
    /*test for hard target*/
    if(decon->gaussFilt){
      if((decon->nGauss==1)&&(decon->gPar[2]<=decon->hardWidth))hardTarget=1;
      else hardTarget=checkHardEnergy(&decon->nGauss,decon->gPar,decon->hardWidth);
    }else hardTarget=0;

    /*copy Gaussian if to be used*/
    if(hardTarget){ /*single hit*/
      TIDY(gaussWave);
      TIDY(smoothed);
      gaussWave=hardHitWave(decon,waveLen);
    }else if(decon->fitGauss){ /*pass on Gaussian waveform*/
      TIDY(smoothed);
    }else{                /*delete fitting and pass original*/
      TIDY(gaussWave);
      gaussWave=smoothed;
      smoothed=NULL;
    }
  }else{   /*don't fit Gaussians*/
    gaussWave=smoothed;
    smoothed=NULL;
    hardTarget=0;
  }

  /*deconvolve if required*/
  if((decon->deconMeth>=0)&&(hardTarget==0)){
    processed=deconvolve(gaussWave,waveLen,decon->pulse,decon->pBins,\
                  decon->res,decon->maxIter,decon->deChang,decon->deconMeth);
    TIDY(gaussWave);
  }else{
    processed=gaussWave;    /*otherwise just use the denoised array*/
    gaussWave=NULL;
  }

  return(processed);
}/*processWave*/


/*##################################################*/
/*process floating point waveform to denoise and deconvolve*/

float *processFloWave(float *wave,int waveLen,denPar *decon,float gbic)
{
  int i=0;
  float *temp=NULL;
  float *mediated=NULL;
  float *preSmoothed=NULL;
  float *denoised=NULL;
  float *smoothed=NULL;
  float *processed=NULL;
  float *denoise(float,float,int,int,float *,float,char);
  float *medianFloat(float *,int,int);
  float *deconvolve(float *,int,float **,int,float,int,double,char);
  float thisTail=0;        /*tail threshold to use here*/
  float *gaussWave=NULL;
  float *fitGaussians(float *,int,denPar *);
  float *hardHitWave(denPar *,int);
  float *sampled=NULL;
  float *digitise(float *,int,char);
  void medNoiseStats(float *,uint32_t,float *,float *,float *,float,float,char);
  void meanNoiseStats(float *,uint32_t,float *,float *,float *,float,float,int);
  char checkHardEnergy(int *,float *,float);
  char hardTarget=0;


  /*determine noise statistics*/
  if(decon->varNoise){
    if(decon->medStats){
      /*convert to a set bit rate*/
      sampled=digitise(wave,waveLen,decon->bitRate);
      medNoiseStats(sampled,waveLen,&(decon->meanN),&(decon->thresh),&thisTail,decon->tailThresh,decon->threshScale,decon->bitRate);
      TIDY(sampled);
    }else meanNoiseStats(wave,waveLen,&(decon->meanN),&(decon->thresh),&thisTail,decon->tailThresh,decon->threshScale,(int)(decon->statsLen/decon->res));
  }else thisTail=decon->tailThresh;
  if(thisTail<0)thisTail=decon->thresh;

  /*convert to a float array for ease*/
  temp=falloc(waveLen,"presmoothed",0);
  for(i=0;i<waveLen;i++)temp[i]=wave[i];

  /*median filter if needed*/
  if(decon->medLen>0){
    mediated=medianFloat(temp,decon->medLen,waveLen);
    TIDY(temp);
  }else{
    mediated=temp;
    temp=NULL;
  }

  /*smooth before denoising*/
  if(decon->psWidth>0.0){
    preSmoothed=smooth(decon->psWidth,waveLen,mediated,decon->res);
    TIDY(mediated);
  }else{
    preSmoothed=mediated;
    mediated=NULL;
  }

  /*remove background noise*/
  if((decon->meanN>0.0)||(decon->thresh>0.0)){
    denoised=denoise(decon->meanN,decon->thresh,decon->minWidth,waveLen,preSmoothed,thisTail,decon->noiseTrack);
    TIDY(preSmoothed);
  }else{
    denoised=preSmoothed;
    preSmoothed=NULL;
  }

  /*smooth if required. Note that pulse is smoothed in readPulse()*/
  if(decon->sWidth>0.0){
    smoothed=smooth(decon->sWidth,waveLen,denoised,decon->res);
    TIDY(denoised);
  }else{
    smoothed=denoised;
    denoised=NULL;
  }

  /*scale by GBIC*/
  if((gbic>0.0)&&(gbic!=1.0))for(i=0;i<waveLen;i++)smoothed[i]/=gbic;


  /*Gaussian fitting*/
  if(decon->fitGauss||decon->gaussFilt){
    gaussWave=fitGaussians(smoothed,waveLen,decon);
    /*test for hard target*/
    if(decon->gaussFilt){
      if((decon->nGauss==1)&&(decon->gPar[2]<=decon->hardWidth))hardTarget=1;
      else hardTarget=checkHardEnergy(&decon->nGauss,decon->gPar,decon->hardWidth);
    }else hardTarget=0;

    /*copy Gaussian if to be used*/
    if(hardTarget){ /*single hit*/
      TIDY(gaussWave);
      TIDY(smoothed);
      gaussWave=hardHitWave(decon,waveLen);
    }else if(decon->fitGauss){ /*pass on Gaussian waveform*/
      TIDY(smoothed);
    }else{                /*delete fitting and pass original*/
      TIDY(gaussWave);
      gaussWave=smoothed;
      smoothed=NULL;
    }
  }else{   /*don't fit Gaussians*/
    gaussWave=smoothed;
    smoothed=NULL;
    hardTarget=0;
  }

  /*deconvolve if required*/
  if((decon->deconMeth>=0)&&(hardTarget==0)){
    processed=deconvolve(gaussWave,waveLen,decon->pulse,decon->pBins,\
                  decon->res,decon->maxIter,decon->deChang,decon->deconMeth);
    TIDY(gaussWave);
  }else{
    processed=gaussWave;    /*otherwise just use the denoised array*/
    gaussWave=NULL;
  }

  return(processed);
}/*processFloWave*/


/*####################################################*/
/*digitse*/

float *digitise(float *wave,int nBins,char bitRate)
{
  int i=0;
  int nDN=0;
  float *sampled=NULL;
  float resDN=0,max=0;
  float tot=0,newTot=0;

  sampled=falloc(nBins,"sampled wave",0);

  /*number of bins*/
  nDN=1;
  for(i=0;i<bitRate;i++)nDN*=2;

  /*find total*/
  max=-10000.0;
  for(i=0;i<nBins;i++){
    if(wave[i]>max)max=wave[i];
    tot+=wave[i];
  }
  if(max<=0.0)max=1.0;

  resDN=max/(float)nDN;
  newTot=0.0;
  for(i=0;i<nBins;i++){
    sampled[i]=floor(wave[i]/resDN+0.5)*resDN;
    newTot+=sampled[i];
  }

  /*rescale energy*/
  if(newTot>0.0){
    for(i=0;i<nBins;i++){
      sampled[i]*=tot/newTot;
    }
  }

  return(sampled);
}/*digitise*/


/*##########################################*/
/*see if one return dominates*/

char checkHardEnergy(int *nGauss,float *gPar,float hardWidth)
{
  int i=0,maxInd=0;
  float total=0,max=0;
  float energy=0,thresh=0;
  char hardHit=0;

  max=-100.0;
  total=0.0;
  for(i=0;i<(*nGauss);i++){
    energy=sqrt(2.0*M_PI)*gPar[3*i+1]*gPar[3*i+2];
    if(energy>max){
      max=energy;
      maxInd=i;
    }
    total+=energy;
  }

  thresh=0.97;   /*97% of energy in a single hit*/
  if((max>=(thresh*total))&&(gPar[3*maxInd+2]<=hardWidth)){
    hardHit=1;
    gPar[0]=gPar[3*maxInd];   /*overwrite Gaussian parameters with brightest*/
    gPar[1]=gPar[3*maxInd+1];
    gPar[2]=gPar[3*maxInd+2];
    (*nGauss)=1;
  }else hardHit=0;

  return(hardHit);
}/*checkHardEnergy*/


/*##########################################*/
/*single hard return*/

float *hardHitWave(denPar *decon,int numb)
{
  int i=0;
  float *hardWave=NULL;

  /*set to 0*/
  hardWave=falloc(numb,"hard waveform",0);
  for(i=0;i<numb;i++)hardWave[i]=0.0;

  /*the single return*/
  if(decon->gPar[0]>=0.0){
    i=(int)(decon->gPar[0]/decon->res);
    if((i>=0)&&(i<numb))hardWave[i]=1.0; /*sqrt(2.0*M_PI)*decon->gPar[1]*decon->gPar[2];*/
    else{
      fprintf(stderr,"Beyond bounds %f %f\n",decon->gPar[0],(float)i*decon->res);
    }
  }

  return(hardWave);
}/*hardHitWave*/


/*##########################################*/
/*correct for attenuation*/

float *correctAttenuation(float *denoised,int numb)
{
  int i=0;
  float tot=0,gap=0;
  float *trueArea=NULL;

  /*count up energy to nromalise visible to unity*/
  tot=0.0;
  for(i=0;i<numb;i++)tot+=denoised[i];

  /*count up gap and apply correction*/
  trueArea=falloc(numb,"",0);
  gap=1.0;
  for(i=0;i<numb;i++){
    if(gap>0.0)trueArea[i]=denoised[i]/gap;
    else       trueArea[i]=0.0;
    gap-=denoised[i]/tot;
  }

  for(i=0;i<numb;i++){
    trueArea[i]/=tot;
    if(trueArea[i]>1.0)trueArea[i]=1.0;
  }

  return(trueArea);
}/*correctAttenuation*/


/*##################################################*/
/*fit Gaussians*/

float *fitGaussians(float *wave,int waveLen,denPar *decon)
{
  int i=0;
  float *fitMultiGauss(float *,float *,int,float,int *,float);
  float *gaussWave=NULL;
  float *x=NULL;

  /*fit Gaussians. Parameters are packed at array end */
  x=falloc(waveLen,"x",0);
  for(i=0;i<waveLen;i++)x[i]=(float)i*decon->res;  /*Gaussian fitting is base 1*/
  decon->nGauss=0;
  gaussWave=fitMultiGauss(x,wave,waveLen,decon->gWidth,&(decon->nGauss),decon->minGsig);
  TIDY(x);

  /*transfer Gaussians parameters*/
  decon->gPar=falloc(3*decon->nGauss,"Gaussian parameters",0);
  for(i=3*decon->nGauss-1;i>=0;i--)decon->gPar[i]=gaussWave[i+waveLen];

  /*trim the fitted wave*/
  if(!(gaussWave=(float *)realloc(gaussWave,waveLen*sizeof(float)))){
    fprintf(stderr,"Error reallocating %lu\n",waveLen*sizeof(float));
    exit(1);
  }

  return(gaussWave);
}/*fitGaussians*/


/*##################################################*/
/*waveform stats*/

void medNoiseStats(float *wave,uint32_t waveLen,float *meanN,float *thresh,float *tailThresh,float tailOff,float threshScale,char bitRate)
{
  uint32_t i=0;
  int *hist=NULL;
  int maxHist=0;
  int nDN=0,ind=0;
  float hRes=0;
  float min=0,max=0;
  float modQuart=0;
  void modalDeviation(float,float *,uint32_t,float *,int,float);

  nDN=1;
  for(i=0;i<bitRate;i++)nDN*=2;

  min=10000000.0;
  max=-10000000.0;
  for(i=0;i<waveLen;i++){
    if(wave[i]>max)max=wave[i];
    if(wave[i]<min)min=wave[i];
  }
  hRes=(max-min)/(float)(nDN-1);
  if(hRes<=0.0){
    min=0.0;
    max=1.0;
    hRes=1.0;   /*it is a blank wave*/
  }

  /*calculate a histogram and work out the modal noise value*/
  hist=ialloc(nDN+1,"noise histogram",0);
  maxHist=-1000;
  for(i=0;i<waveLen;i++){
    ind=(int)((wave[i]-min)/hRes);
    hist[ind]++;
    if(hist[ind]>maxHist){
      maxHist=hist[ind];
      (*meanN)=wave[i];  /*modal noise*/
    }
  }

  /*modal quartile*/
  modalDeviation(*meanN,wave,waveLen,&modQuart,nDN,hRes);
  if(modQuart<2)modQuart=2;

  (*thresh)=(*meanN)+threshScale*(float)modQuart;
  if((*thresh)<14.0)(*thresh)=14.0;
  (*tailThresh)=(*thresh)+tailOff;

  TIDY(hist);
  return;
}/*medNoiseStats*/


/*############################################*/
/*mean noise statistics*/

void meanNoiseStats(float *sampled,uint32_t waveLen,float *meanN,float *thresh,float *thisTail,float tailThresh,float threshScale,int statBins)
{
  int i=0,start=0;
  float stdev=0;

  start=2;
  if((uint32_t)statBins>waveLen){
    fprintf(stderr,"Not enough bins for this statistics length\n");
    exit(1);
  }else if((statBins-start)<=0){
    fprintf(stderr,"What are you doing!?!\n");
    exit(1);
  }

  /*NOTE, subtract 2 from noise bins to avoid low smoothing at signal start*/
  (*meanN)=0.0;
  for(i=start;i<statBins;i++)(*meanN)+=sampled[i];
  (*meanN)/=(float)(statBins-start);
  stdev=0.0;
  for(i=start;i<statBins;i++)stdev+=((*meanN)-sampled[i])*((*meanN)-sampled[i]);
  stdev=sqrt(stdev/(float)(float)(statBins-start));


  (*thresh)=(*meanN)+threshScale*stdev;
  if(tailThresh<=0.0)(*thisTail)=(*thresh);
  else               (*thisTail)=(*thresh)+tailThresh;

  return;
}/*meanNoiseStats*/



/*############################################*/
/*deviation from the mode histogram*/

void modalDeviation(float modal,float *wave,uint32_t waveLen,float *modQuart,int nDN,float hRes)
{
  int i=0;
  int *devHist=NULL;
  int bin=0,total=0;
  int quartThresh=0;
  int cumul=0,subTot=0;

  total=subTot=0;
  devHist=ialloc(nDN,"",0);
  for(i=0;i<nDN;i++)devHist[i]=0;
  for(i=0;i<waveLen;i++){
    bin=(int)((wave[i]-(int)modal)/hRes+0.5);
    if(bin>=0){
      devHist[bin]++;
      subTot++;
    }
    total++;
  }

  /*calculate mode and top quartile*/
  quartThresh=(int)((float)total*3.0/4.0);
  cumul=total-subTot;
  for(i=0;i<nDN;i++){
    cumul+=devHist[i];
    if(cumul>=quartThresh){
      (*modQuart)=(unsigned char)i;
      break;
    }
  }
  TIDY(devHist);
  return;
}/*modalDeviation*/


/*##################################################*/
/*smooth waveform*/

float *smooth(float sWidth,int nBins,float *data,float res)
{
  int i=0,j=0;
  int tP=0;
  int p1=0,p2=0;         /*array places*/
  int nPulse=0;          /*number of pulse bins*/
  float *smoothed=NULL;
  float energy=0;
  float *setPulse(float,int *,float);
  float *markFloat(int,float *,float);
  int *markInt(int,int *,int);
  char newPulse=0;

  /*smooth as required*/
  smoothed=falloc(nBins,"",0);
  if(sWidth<TOL)memcpy(smoothed,data,sizeof(float)*nBins);
  else{
    /*test to see if a new pulse array is needed*/
    newPulse=1;
    for(tP=0;tP<smooPulse.nPulses;tP++){
      if(fabs(sWidth-smooPulse.sWidth[tP])<TOL){
        newPulse=0;
        break;
      }
    }/*new pulse test*/

    /*create new pulse if needed*/
    if(newPulse){
      tP=smooPulse.nPulses;
      if(!(smooPulse.pulse=(float **)realloc(smooPulse.pulse,(smooPulse.nPulses+1)*sizeof(float *)))){
        fprintf(stderr,"Error reallocating %lu\n",(smooPulse.nPulses+1)*sizeof(float *));
        exit(1);
      }
      smooPulse.pulse[tP]=setPulse(sWidth,&nPulse,res);
      smooPulse.sWidth=markFloat(smooPulse.nPulses,smooPulse.sWidth,sWidth);
      smooPulse.nBins=markInt(smooPulse.nPulses,smooPulse.nBins,nPulse);
      smooPulse.nPulses++;
    }

    for(i=0;i<nBins;i++){
      smoothed[i]=smooPulse.pulse[tP][0]*data[i];
      energy=smooPulse.pulse[tP][0];
      for(j=1;j<smooPulse.nBins[tP];j++){
        p1=i-j;
        p2=i+j;
        if((p1>=0)&&(p1<nBins)){
          smoothed[i]+=smooPulse.pulse[tP][j]*data[p1];
          energy+=smooPulse.pulse[tP][j];
        }
        if((p2>=0)&&(p2<nBins)){
          smoothed[i]+=smooPulse.pulse[tP][j]*data[p2];
          energy+=smooPulse.pulse[tP][j];
        }
      }
      smoothed[i]/=energy;
    }/*bin loop*/
  }/*smooth check*/

  return(smoothed);
}/*smooth*/


/*##################################################*/
/*tidy up smoothing pulse arrays*/

void tidySMoothPulse()
{
  TIDY(smooPulse.nBins);
  TIDY(smooPulse.sWidth);
  TTIDY((void **)smooPulse.pulse,smooPulse.nPulses);

  return;
}/*tidySMoothPulse*/


/*##################################################*/
/*set smoothing pulse*/

float *setPulse(float sWidth,int *nPulse,float res)
{
  int i=0;
  float *pulse=NULL;
  float y=0,x=0;
  float minY=0;
  double gaussian(double,double,double);

  minY=0.0001;
  for(i=0;;i++){
    x=(float)i*res;
    y=(float)gaussian((double)x,(double)sWidth,0.0);
    if(y<=minY)break;
  }

  (*nPulse)=i;
  pulse=falloc(*nPulse,"smoothing pulse",0);
  for(i=0;i<(*nPulse);i++){
    x=(float)i*res;
    pulse[i]=(float)gaussian((double)x,(double)sWidth,0.0);
  }

  return(pulse);
}/*setPulse*/


/*##################################################*/
/*deconvolve using Gold's method*/

float *deconvolve(float *data,int nBins,float **pulse,int pBins,float res,int maxIter,double minChange,char meth)
{
  int i=0,numb=0;
  float *decon=NULL;
  double *deconDo=NULL;
  double *pulseDo=NULL;
  double *dataDo=NULL;
  double *resamplePulse(int,float **,float,int);
  double *goldMeth(double *,double *,int,int,double);
  double *richLucy(double *,double *,int,int,double);
  double energy=0,newEn=0;   /*to balance energies*/


  /*arrays have to be of base 2 length*/
  numb=pow(2.0,(float)((int)(log((double)nBins)/log(2.0)+0.5)+1));
  dataDo=dalloc(numb,"dataDo",0);

  energy=0.0;
  for(i=0;i<nBins;i++){
    dataDo[i]=(double)data[i];
    energy+=dataDo[i];
  }

  /*pulse needs resampling to the correct resolution*/
  pulseDo=resamplePulse(numb,pulse,res,pBins);

  /*call deconvolution method*/
  if(meth==0)     deconDo=goldMeth(dataDo,pulseDo,numb,maxIter,minChange);
  else if(meth==1)deconDo=richLucy(dataDo,pulseDo,numb,maxIter,minChange);
  else{
    fprintf(stderr,"Deconvolution method not defined\n");
    exit(1);
  }
  TIDY(pulseDo);  /*tidy as we go along*/
  TIDY(dataDo);

  /*transfer data*/
  decon=falloc(nBins,"",0);
  newEn=0;
  for(i=0;i<nBins;i++)newEn+=deconDo[i];
  for(i=0;i<nBins;i++)decon[i]=(float)(deconDo[i]*energy/newEn);

  TIDY(deconDo);
  return(decon);
}/*deconvolve*/


/*##################################################*/
/*Richardson-Lucy deconvolution*/

double *richLucy(double *data,double *pulse,int numb,int maxIter,double minChange)
{
  int i=0,j=0;
  int real=0,imag=0;
  double *o=NULL,*work=NULL;
  double *smooth=NULL,*denom=NULL;
  double tot=0,changeSq=0,minChangeSq=0;
  double new=0;
  int gsl_fft_complex_radix2_forward(gsl_complex_packed_array,size_t,size_t);
  int gsl_fft_complex_radix2_backward(gsl_complex_packed_array, size_t,size_t);

  minChangeSq=minChange*minChange;

  /*perform FFT*/
  gsl_fft_complex_radix2_forward((gsl_complex_packed_array)pulse,1,numb);

  /*real arrays*/
  o=dalloc(numb,"o",0);
  denom=dalloc(numb,"denominator",0);
  /*complex arrays*/
  work=dalloc(2*numb,"workspace",0);
  smooth=dalloc(2*numb,"workspace",0);

  /*initial guess*/
  for(i=0;i<numb;i++){
    o[i]=data[i];
    tot+=data[i];
  }

  do{  /*the iterative deconvolution*/
    /*convolve pulse and estimate*/
    for(i=0;i<numb;i++){
      work[2*i]=o[i];
      work[2*i+1]=0.0;
    }
    gsl_fft_complex_radix2_forward((gsl_complex_packed_array)work,1,numb);
    for(i=0;i<numb;i++){
      real=i*2;
      imag=i*2+1;
      smooth[real]=pulse[real]*work[real]-pulse[imag]*work[real];
      smooth[imag]=pulse[imag]*work[real]+pulse[real]*work[imag];
    }
    gsl_fft_complex_radix2_backward((gsl_complex_packed_array)smooth,1,numb);

    /*calculate new estimate*/
    for(i=0;i<numb;i++){
      real=i*2;
      imag=i*2+1;
      denom[i]=sqrt(smooth[real]*smooth[real]+smooth[imag]*smooth[imag]);
      if(denom[i]>0.0)work[real]=data[i]/denom[i];
      else            work[real]=0.0;
      work[imag]=0.0;
    }
    gsl_fft_complex_radix2_forward((gsl_complex_packed_array)work,1,numb);
    for(i=0;i<numb;i++){
      real=2*i;
      imag=2*i+1;
      smooth[real]=pulse[real]*work[real]-pulse[imag]*work[real];
      smooth[imag]=pulse[imag]*work[real]+pulse[real]*work[imag];
    }
    gsl_fft_complex_radix2_backward((gsl_complex_packed_array)smooth,1,numb);

    tot=0.0;
    changeSq=0.0;
    for(i=0;i<numb;i++){
      real=2*i;
      imag=2*i+1;
      new=o[i]*sqrt(smooth[real]*smooth[real]+smooth[imag]*smooth[imag]);
      changeSq+=(o[i]-new)*(o[i]-new);
      o[i]=new;
      tot+=o[i];
    }
    changeSq/=(double)numb;
    /*fprintf(stdout,"Iter %d change %.20f\n",j,changeSq);*/
    if((minChange>=0.0)&&(changeSq<=minChangeSq))break;
    j++;
  }while(j<maxIter);

  TIDY(work);
  TIDY(smooth);
  TIDY(denom);
  return(o);
}/*richLucy*/


/*##################################################*/
/*Gold's method*/

double *goldMeth(double *data,double *pulse,int numb,int maxIter,double minChange)
{
  int i=0,j=0;
  int real=0,imag=0;
  double *o=NULL,*work=NULL;
  double *smooth=NULL,*denom=NULL;
  double tot=0,changeSq=0,minChangeSq=0;
  double new=0;
  int gsl_fft_complex_radix2_forward(gsl_complex_packed_array,size_t,size_t);
  int gsl_fft_complex_radix2_backward(gsl_complex_packed_array, size_t,size_t);

  minChangeSq=minChange*minChange;

  /*perform FFT*/
  gsl_fft_complex_radix2_forward((gsl_complex_packed_array)pulse,1,numb);

  /*real arrays*/
  o=dalloc(numb,"o",0);
  denom=dalloc(numb,"denominator",0);
  /*complex arrays*/
  work=dalloc(2*numb,"workspace",0);
  smooth=dalloc(2*numb,"workspace",0);

  /*initial guess*/
  for(i=0;i<numb;i++)o[i]=data[i];

  /*iterate over Gold's method*/
  do{
    /*fourier current estimate*/
    for(i=0;i<numb;i++){
      work[2*i]=o[i];
      work[2*i+1]=0.0;
    }
    gsl_fft_complex_radix2_forward((gsl_complex_packed_array)work,1,numb);

    /*blur with pulse*/
    for(i=0;i<numb;i++){
      real=i*2;
      imag=i*2+1;
      smooth[real]=pulse[real]*work[real]-pulse[imag]*work[real];
      smooth[imag]=pulse[imag]*work[real]+pulse[real]*work[imag];
    }
    gsl_fft_complex_radix2_backward((gsl_complex_packed_array)smooth,1,numb);

    /*reblur deoniminator*/
    tot=0.0;
    changeSq=0.0;
    for(i=0;i<numb;i++){
      real=i*2;
      imag=i*2+1;
      denom[i]=sqrt(smooth[real]*smooth[real]+smooth[imag]*smooth[imag]);
      if(denom[i]>0.0){
        new=o[i]*data[i]/denom[i];
        tot+=new;
      }else new=0.0;
      changeSq+=(o[i]-new)*(o[i]-new);
      o[i]=new;
    }
    changeSq/=(double)numb;
    /*fprintf(stdout,"Iter %d change %.20f\n",j,changeSq);*/
    if((minChange>=0.0)&&(changeSq<=minChangeSq))break;
    j++;
  }while(j<maxIter);

  TIDY(work);
  TIDY(smooth);
  TIDY(denom);

  return(o);
}/*goldMeth*/


/*##################################################*/
/*denoise by thresholding*/

float *denoise(float meanN,float thresh,int minWidth,int nBins,float *data,float tailThresh,char noiseTrack)
{
  int i=0,j=0;
  int start=0;
  float *denoised=NULL;
  float thisMean=0;
  char waveStart=0;

  denoised=falloc(nBins,"denoised",0);
  for(i=0;i<nBins;i++)denoised[i]=0.0;
  waveStart=0;

  if(noiseTrack)thisMean=meanN;
  else          thisMean=thresh;

  /*threshold*/
  start=-1;
  for(i=0;i<nBins;i++){
    if((data[i]>=tailThresh)||((waveStart==0)&&(data[i]>=thresh))){        /*check whether we're within a feature*/
      if(start<0)start=i;       /*mark start*/
      waveStart=1;
    }else if(start>=0){         /*left a feature*/
      if((i-start)>=minWidth){  /*check feature width*/
        for(j=start;j>=0;j--){      /*noise tracking*/
          if(data[j]<=thisMean)break;
          denoised[j]=data[j]-meanN;
          if(denoised[j]<0.0)denoised[j]=0.0;  /*force positive to avoid nan in Gold*/
        }
        for(i=start+1;i<nBins;i++){
          if(data[i]<=thisMean)break;
          denoised[i]=data[i]-meanN;
          if(denoised[i]<0.0)denoised[i]=0.0;  /*force positive to avoid nan in Gold*/
        }
      }
      start=-1;  /*reset counter*/
    }/*within feature check*/
  }/*threshold*/

  return(denoised);
}/*denoise*/


/*##################################################*/
/*resample pulse to complex array at correct res*/

double *resamplePulse(int numb,float **pulse,float res,int pBins)
{
  int i=0,bin=0,step=0,maxBin=0;
  float max=0,maxRange=0;
  double *pulseDo=NULL,total=0;

  pulseDo=dalloc(2*numb,"pulseDo",0);
  for(i=2*numb-1;i>=0;i--)pulseDo[i]=0.0;

  /*find the peak to set at zero as pulse is aligned by CofG*/
  max=-1000.0;
  for(i=0;i<pBins;i++){
    if(pulse[1][i]>max){
      max=pulse[1][i];
      maxRange=pulse[0][i];
      maxBin=i;
    }
  }/*peak finding*/

  /*find nearest pulse point to bin, start from centre*/
  step=(int)(res/(pulse[0][1]-pulse[0][0]));
  total=0.0;
  for(i=maxBin;i<pBins;i+=step){
    bin=(int)((pulse[0][i]-maxRange)/res);
    if(bin<0)bin+=numb;
    pulseDo[2*bin]=pulse[1][i];
    total+=pulse[1][i];
  }/*from centre up*/
  for(i=maxBin-1;i>=0;i-=step){  /*then work from the centre backwards*/
    bin=(int)((pulse[0][i]-maxRange)/res);
    if(bin<0)bin+=numb;
    pulseDo[2*bin]=pulse[1][i];
    total+=pulse[1][i];
  }/*from centre back*/

  /*normalise pulse and prevent zeroes*/
  for(i=0;i<numb;i++){
    if(pulseDo[2*i]<0.0)pulseDo[2*i]=0.0;
    else                pulseDo[2*i]/=total;
  }/*normalisation*/

  return(pulseDo);
}/*reamplePulse*/


/*##################################################*/
/*read pulse shape*/

void readPulse(denPar *denoise)
{
  int i=0,nGauss=0;
  char line[200];
  char temp[2][100];
  FILE *ipoo=NULL;
  /*smoothing*/
  float *smoothed=NULL;
  float *tempPulse=NULL;
  float *fitSingleGauss(float *,float *,int,float,int *,float **);
  float *gaussPar=NULL;

  if((ipoo=fopen(denoise->pNamen,"r"))==NULL){
    fprintf(stderr,"Error opening pulse file %s\n",denoise->pNamen);
    exit(1);
  }

  denoise->pBins=0;
  /*count number of bins*/
  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1))denoise->pBins++;
  }/*line counting*/


  /*rewind*/
  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error\n");
    exit(1);
  }

  denoise->pulse=fFalloc(2,"",0);
  for(i=0;i<2;i++)denoise->pulse[i]=falloc(denoise->pBins,"",i+1);

  /*read data*/
  i=0;
  while(fgets(line,200,ipoo)!=NULL){
    if(strncasecmp(line,"#",1)){
      if(sscanf(line,"%s %s",temp[0],temp[1])==2){
        if(i>denoise->pBins){
          fprintf(stderr,"Error\n");
          exit(1);
        }
        denoise->pulse[0][i]=atof(&(temp[0][0]))*denoise->pScale;
        denoise->pulse[1][i]=atof(&(temp[1][0]));
        i++;
      }
    }/*comment check*/
  }/*data reading*/

  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  /*smooth if required*/
  if(denoise->sWidth>0.0){
    smoothed=smooth(denoise->sWidth,denoise->pBins,&(denoise->pulse[1][0]),denoise->pulse[0][1]-denoise->pulse[0][0]);
    TIDY(denoise->pulse[1])
    denoise->pulse[1]=smoothed;
    smoothed=NULL;
  }/*smoothing*/

  /*fit a Gaussian or determine width for hard limits if required*/
  if((denoise->gaussPulse)||(denoise->gaussFilt)){
    nGauss=0;
    tempPulse=fitSingleGauss(denoise->pulse[0],denoise->pulse[1],denoise->pBins,0.5,&(nGauss),&gaussPar);

    /*set width if needed*/
    if(denoise->gaussFilt){
      if(nGauss>1){
        fprintf(stderr,"Multiple Gaussians fitted to pulse\n");
        exit(1);
      }
      denoise->hardWidth=gaussPar[2]*denoise->hardTol;
    }
    TIDY(gaussPar);

    /*copy Gaussian over pulse if needed*/
    if(denoise->gaussPulse){
      TIDY(denoise->pulse[1]);
      denoise->pulse[1]=tempPulse;
      tempPulse=NULL;
    }else TIDY(tempPulse);
  }

  return;
}/*readPulse*/


/*########################################################################*/
/*set default denoising parameters*/

void setDenoiseDefault(denPar *denoise)
{
  /*denoising*/
  denoise->meanN=12.0;
  denoise->tailThresh=-1.0;
  denoise->thresh=15.0;
  denoise->minWidth=6;
  denoise->sWidth=0.0;
  denoise->psWidth=0.0;
  denoise->medLen=0;
  denoise->varNoise=0;
  denoise->medStats=0;
  denoise->statsLen=30.0;
  denoise->noiseTrack=1;
  denoise->threshScale=1.0;
  denoise->bitRate=10;

  /*deconvolution*/
  denoise->deconMeth=-1;     /*do not deconvolve*/
  denoise->pScale=1.0;      /*scale pulse length by*/
  denoise->maxIter=2000;     /*maximum number of iterations*/
  denoise->deChang=pow(10,0-7.0);  /*change between decon iterations to stop*/
  strcpy(denoise->pNamen,"/home/sh563/data/bess/leica_shape/leicaPulse.dat");  /*pulse filename*/
  denoise->pulse=NULL;       /*pulse to deconvolve by*/
  denoise->pBins=0;          /*number of pulse bins*/
  denoise->res=0.15;

  /*Gaussian fitting*/
  denoise->gWidth=1.5;
  denoise->fitGauss=0;       /*do not fit Gaussians*/
  denoise->gaussPulse=0;     /*do not turn pulse to Gaussian*/
  denoise->minGsig=0.00001;  /*minimum Gaussian fitting width*/

  /*Gaussian hard target identification*/
  denoise->gaussFilt=0;    /*switch*/
  denoise->hardWidth=0.0;  /*maxWidth of hard feature*/
  denoise->hardTol=1.0;    /*tolerance to scale width by*/

  return;
}/*setDenoiseDefault*/


/*########################################################################*/
/*rotate about x axis*/

void rotateX(double *vect,double theta)
{
  int i=0;
  double temp[3];

  temp[0]=vect[0];
  temp[1]=vect[1]*cos(theta)+vect[2]*sin(theta);
  temp[2]=vect[2]*cos(theta)-vect[1]*sin(theta);

  for(i=0;i<3;i++)vect[i]=temp[i];
  return;
}/*rotateX*/


/*########################################################################*/
/*rotate about y axis*/

void rotateY(double *vect,double theta)
{
  int i=0;
  double temp[3];

  temp[0]=vect[0]*cos(theta)-vect[1]*sin(theta);
  temp[1]=vect[1];
  temp[2]=vect[0]*sin(theta)+vect[2]*cos(theta);

  for(i=0;i<3;i++)vect[i]=temp[i];
  return;
}/*rotateY*/


/*########################################################################*/
/*rotate about z axis*/

void rotateZ(double *vect,double theta)
{
  int i=0;
  double temp[3];

  temp[0]=vect[0]*cos(theta)+vect[1]*sin(theta);
  temp[1]=vect[1]*cos(theta)-vect[0]*sin(theta);
  temp[2]=vect[2];

  for(i=0;i<3;i++)vect[i]=temp[i];
  return;
}/*rotateZ*/


/*########################################################################*/
/*check point bounds*/

char boundsCheck(double x,double y,double z,double *bounds)
{
  if((x>=(bounds[0]))&&(y>=(bounds[1]))&&(z>=(bounds[2]))&&\
     (x<=(bounds[3]))&&(y<=(bounds[4]))&&(z<=(bounds[5])))return(1);
  else                                                    return(0);
}/*boundsCheck*/

/*the end*/
/*################################################*/

