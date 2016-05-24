#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "stdint.h"
#include "geotiffio.h"
#include "xtiffio.h"


/*library for handling my geotiff bits*/


/*########################################*/
/*write a geoTIFF file*/

void drawTiff(char *namen,double *geoL,int *geoI,double res,unsigned char *image,int nX,int nY,double scale,uint16_t epsg)
{
  void SetUpTIFFDirectory(TIFF *,int,int,int,int,double,double,float,double);
  void SetUpGeoKeys(GTIF *,uint16_t);
  void WriteImage(TIFF *,int,int,unsigned char *);
  TIFF *tif=(TIFF*)0;  /* TIFF-level descriptor */
  GTIF *gtif=(GTIF*)0; /* GeoKey-level descriptor */

  tif=XTIFFOpen(namen,"w");
  if(!tif){
    fprintf(stderr,"Error opening %s\n",namen);
    exit(1);
  }
  gtif=GTIFNew(tif);
  if(!gtif){
    fprintf(stderr,"failed in GTIFNew\n");
    exit(1);
  }
  SetUpTIFFDirectory(tif,nX,nY,geoI[0],geoI[1],geoL[0],geoL[1],res,scale);
  SetUpGeoKeys(gtif,epsg);
  WriteImage(tif,nX,nY,image);
  GTIFWriteKeys(gtif);
  GTIFFree(gtif);
  XTIFFClose(tif);
  fprintf(stdout,"Written to %s\n",namen);
  tif=NULL;
  gtif=NULL;

  return;
}/*drawTiff*/


/*##############################################*/
/*set up tiff difrectory*/

void SetUpTIFFDirectory(TIFF *tif,int nX,int nY,int xI,int yI,double xL,double yL,float res,double scale)
{
  double tiepoints[6];
  double pixscale[3];


  pixscale[0]=(double)res;
  pixscale[1]=(double)res;
  pixscale[2]=scale;

  tiepoints[0]=(double)xI;
  tiepoints[1]=(double)yI;
  tiepoints[2]=0.0;
  tiepoints[3]=xL;
  tiepoints[4]=yL;
  tiepoints[5]=0.0;



  TIFFSetField(tif,TIFFTAG_IMAGEWIDTH,    nX);
  TIFFSetField(tif,TIFFTAG_IMAGELENGTH,   nY);
  TIFFSetField(tif,TIFFTAG_COMPRESSION,   COMPRESSION_NONE);
  TIFFSetField(tif,TIFFTAG_PHOTOMETRIC,   PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif,TIFFTAG_PLANARCONFIG,  PLANARCONFIG_CONTIG);
  TIFFSetField(tif,TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(tif,TIFFTAG_ROWSPERSTRIP,  nX);

  TIFFSetField(tif,TIFFTAG_GEOTIEPOINTS, 6,tiepoints);
  TIFFSetField(tif,TIFFTAG_GEOPIXELSCALE, 3,pixscale);
  return;
}/*SetUpTIFFDirectory*/


/*##############################################*/
/*set up geokeys*/

void SetUpGeoKeys(GTIF *gtif,uint16_t epsg)
{
  GTIFKeySet(gtif,GTModelTypeGeoKey,TYPE_SHORT, 1, ModelTypeProjected);
  //GTIFKeySet(gtif,ProjectedCSTypeGeoKey,TYPE_SHORT,1,PCS_British_National_Grid );
  GTIFKeySet(gtif,ProjectedCSTypeGeoKey,TYPE_SHORT,1,epsg); //26917);
  //GTIFKeySet(gtif,33922, TYPE_DOUBLE,6,geokey);   /*tie point coordindates*/
  //GTIFKeySet(gtif,33550, TYPE_FLOAT,3,reskey);         /*scale between image and world*/
  //GTIFKeySet(gtif,33550, TYPE_DOUBLE,1,1,1);         /*scale between image and world*/


  //GTIFKeySet(gtif,GeogSemiMajorAxisGeoKey, TYPE_DOUBLE,1, (double)6377298.556);
  //GTIFKeySet(gtif,GeogInvFlatteningGeoKey, TYPE_DOUBLE,1, (double)300.8017);
  return;
}/*SetUpGeoKeys*/


/*##############################################*/
/*write image out*/

void WriteImage(TIFF *tif,int nX,int nY,unsigned char *image)
{
  int i;
  char buffer[nX];

  memset(buffer,0,(size_t)nX);
  for(i=0;i<nY;i++){
    if(!TIFFWriteScanline(tif,&(image[i*nX]),i,0))TIFFError("WriteImage","failure in WriteScanline\n");
  }
  return;
}/*WriteImage*/


/*##############################################*/

void WriteImageFlo(TIFF *tif,int nX,int nY,float *image)
{
  int i=0; //,j=0;

  for(i=0;i<nY;i++){
    if(!TIFFWriteScanline(tif,&(image[i*nX]),i,0))TIFFError("WriteImage","failure in WriteScanline\n");
  }
  return;
}/*WriteImageFlo*/


/*########################################*/
/*write a float32 geoTIFF file*/

void drawTiffFlo(char *outRoot,double *geoL,int *geoI,double res,float *image,int nX,int nY,double scale,uint16_t epsg)
{
  void SetUpTIFFDirectoryFlo(TIFF *,int,int,int,int,double,double,float,double);
  void SetUpGeoKeys(GTIF *,uint16_t);
  void WriteImageFlo(TIFF *,int,int,float *);
  TIFF *tif=(TIFF*)0;  /* TIFF-level descriptor */
  GTIF *gtif=(GTIF*)0; /* GeoKey-level descriptor */
  char namen[200];

  sprintf(namen,"%s.tif",outRoot);
  tif=XTIFFOpen(namen,"w");
  if(!tif){
    fprintf(stderr,"Error opening %s\n",namen);
    exit(1);
  }
  gtif=GTIFNew(tif);
  if(!gtif){
    fprintf(stderr,"failed in GTIFNew\n");
    exit(1);
  }
  SetUpTIFFDirectoryFlo(tif,nX,nY,geoI[0],geoI[1],geoL[0],geoL[1],res,scale);
  SetUpGeoKeys(gtif,epsg);
  WriteImageFlo(tif,nX,nY,image);
  GTIFWriteKeys(gtif);
  GTIFFree(gtif);
  XTIFFClose(tif);
  fprintf(stdout,"Written to %s\n",namen);
  tif=NULL;
  gtif=NULL;

  return;
}/*drawTiffFlo*/


/*##############################################*/
/*set up tiff headers*/

void SetUpTIFFDirectoryFlo(TIFF *tif,int nX,int nY,int xI,int yI,double xL,double yL,float res,double scale)
{
  double tiepoints[6];
  double pixscale[3];


  pixscale[0]=(double)res;
  pixscale[1]=(double)res;
  pixscale[2]=scale;

  tiepoints[0]=(double)xI;
  tiepoints[1]=(double)yI;
  tiepoints[2]=0.0;
  tiepoints[3]=xL;
  tiepoints[4]=yL;
  tiepoints[5]=0.0;


  TIFFSetField(tif,TIFFTAG_IMAGEWIDTH,    nX);
  TIFFSetField(tif,TIFFTAG_IMAGELENGTH,   nY);
  TIFFSetField(tif,TIFFTAG_COMPRESSION,   COMPRESSION_NONE);
  TIFFSetField(tif,TIFFTAG_PHOTOMETRIC,   PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif,TIFFTAG_PLANARCONFIG,  PLANARCONFIG_CONTIG);
  TIFFSetField(tif,TIFFTAG_BITSPERSAMPLE, 32);
  TIFFSetField(tif,TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP); // 3
  TIFFSetField(tif,TIFFTAG_ROWSPERSTRIP,  nX);
  TIFFSetField(tif,TIFFTAG_GEOTIEPOINTS, 6,tiepoints);
  TIFFSetField(tif,TIFFTAG_GEOPIXELSCALE, 3,pixscale);
  return;
}/*SetUpTIFFDirectoryFlo*/

/*the end*/
/*##############################################*/

