

/*######################*/
/*# A library for      #*/
/*# handling las files #*/
/*# S Hancock, 2015    #*/
/*######################*/



/*###########################################*/
/*lasfile point bit field*/

typedef struct{
  unsigned char retNumb: 3;
  unsigned char nRet: 3;
  unsigned char sDir: 1;
  unsigned char edge: 1;
}b_field;


/*###########################################*/
/*a las ARSF file, version 1.0 to 1.3*/

typedef struct{
  /*file pointer*/
  FILE *ipoo;
  char namen[300];   /*filename*/

  /*public header*/
  uint8_t vMajor;    /*version major*/
  uint8_t vMinor;    /*version minor*/
  int headLen;
  uint16_t doy;
  uint16_t year;
  uint16_t headSize;
  uint32_t offsetToP; /*offset to point data*/
  uint32_t nVarRec;   /*number of variable length records*/
  unsigned char pointFormat;   /*point data format ID*/
  uint16_t pRecLen;   /*point record length*/
  uint32_t nPoints;   /*number of point records*/
  uint32_t nPbyRet[7];/*number of point records by return*/
  double posScale[3];
  double posOffset[3];
  double minB[3];
  double maxB[3];
  uint64_t waveStart;

  /*geolocation*/
  uint16_t epsg;       /*projection code*/

  /*point*/
  int32_t x;           /*index to calculate coordinate*/
  int32_t y;           /*index to calculate coordinate*/
  int32_t z;           /*index to calculate coordinate*/
  uint16_t refl;       /*point intensity*/
  b_field field;         /*bit field containing all sorts*/
  unsigned char classif;  /*point classification*/
  unsigned char scanAng;  /*scan angle*/
  unsigned char userData; /*user data*/
  uint16_t psID;          /*point source ID, used for Leica's AGC*/
  double gpsTime;         /*GPS time*/
  unsigned char packetDes;/*waveform packed description*/
  uint64_t waveMap;       /*pointer to waveform in file*/
  uint32_t waveLen;       /*length of waveform in bins*/
  float time;             /*time in picoseconds of this wave*/
  float grad[3];          /*waveform gradient*/
  uint16_t RGB[3];        /*RGB image*/
  /*uint64_t counter;*/    /*I do not know why this is here?*/

  /*buffer to read multiple points at a time*/
  char *pointBuff;       /*buffer with multiple points*/
  uint32_t buffStart;    /*buffer start in number of points*/
  uint32_t buffLen;      /*buffer length in number of points*/
  uint64_t buffByte;     /*buffer length in number of bytes*/
  uint32_t maxBuffLen;   /*max buffer length in number of points*/
  uint64_t maxBuffSize;  /*max buffer length in number of bytes*/

  /*waveform*/
  unsigned char *wave;
}lasFile;


/*the end*/
/*###########################################*/

