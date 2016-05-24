
/*#############################################*/
/*strucure of denoising variables*/

typedef struct{
  /*denoising*/
  float thresh;      /*noise threshold*/
  float meanN;       /*mean noise level*/
  float tailThresh;  /*trailing edge threshold, meanN by default*/
  float sWidth;      /*smoothing length*/
  float psWidth;     /*pre-smoothing width*/
  int minWidth;      /*min width above noise to accept*/
  int medLen;        /*median filter width*/
  char varNoise;     /*variable noise threshold switch*/
  char medStats;     /*use median rather than mean statistics (ARSF)*/
  float statsLen;    /*distance to calculate noise statistics over*/
  char noiseTrack;   /*noise tracking switch*/
  float threshScale; /*scale variable noise threshold*/
  char bitRate;      /*bit rate if using variable noise with floats*/

  /*deconvolution*/
  float pScale;      /*scale pulse length by*/
  int maxIter;       /*maximum number of iterations*/
  double deChang;    /*change between decon iterations to stop*/
  char deconMeth;    /*deconvolution method*/
  char pNamen[200];  /*pulse filename*/
  float **pulse;     /*pulse to deconvolve by*/
  int pBins;         /*number of pulse bins*/
  float res;         /*waveform resolution*/

  /*Gaussian fitting*/
  char fitGauss;     /*switch*/
  float gWidth;      /*smoothing width before picking features*/
  int nGauss;        /*number of Gaussians*/
  float *gPar;       /*Gaussian parameters*/
  char gaussPulse;   /*convert pulse to a Gaussian*/
  float minGsig;     /*minimum Gaussian width to fit*/

  /*Gaussians to identify hard targets*/
  char gaussFilt;    /*switch*/
  float hardWidth;   /*maxWidth of hard feature*/
  float hardTol;     /*tolerance to scale width by*/
}denPar;


/*##################################################*/
/*global structure to save reallocating smoothing pulse*/

typedef struct{
  int nPulses;
  int *nBins;
  float *sWidth;
  float **pulse;
}smoothPulse;


/*#############################################*/
/*common function definitions*/

float *smooth(float,int,float *,float);
smoothPulse smooPulse;   /*global structure to save reallocation*/

/*the end*/
/*#############################################*/

