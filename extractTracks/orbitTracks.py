import h5py
import optparse
import string
import numpy as np
from osgeo import gdal
import random
from pyproj import Proj, transform
from scipy import pi


# intersect MODIS pixels
def intersectLines(x1,x0,y1,y0,xOrigin,yOrigin,pixelWidth,pixelHeight):
  xi=[]
  yi=[]
  xi=np.append(xi,x0)
  yi=np.append(yi,y0)
  tx=x0
  ty=y0

  # loop along
  dx=x1-x0
  dy=y1-y0
  xStep=1 if(dx>0) else -1
  yStep=1 if(dy>0) else -1
  tol=0.00000001

  xDist=1
  yDist=1

  while((xDist>0.0)&(yDist>0.0)):
    xT=(int((tx-xOrigin)/pixelWidth)+xStep)*pixelWidth+xOrigin
    yT=(int((ty-yOrigin)/pixelHeight)+yStep)*pixelHeight+yOrigin
    if(dx!=0.0):
      xDist=(xT-tx)/dx
    else:
      xDist=1000.0
    if(dy!=0.0):
      yDist=(yT-ty)/dy
    else:
      yDist=10000.0;

    if(xDist<=yDist):
      tx=xT
      ty=ty+dy*xDist
    else:
      tx=tx+dx*yDist
      ty=yT

    xi=np.append(xi,tx)
    yi=np.append(yi,ty)

  xi=np.append(xi,x1)
  yi=np.append(yi,y1)

  return(xi,yi)


# read command line
def clParser():
    desc = """Script extracts the MODIS binary layers from the HDFs and dumps and envi (*hopefully*) hdr file"""
    clp = optparse.OptionParser("Usage: %prog [options] filename", description = desc)
    clp.add_option("--input",type="string",dest="inNamen",help="input file name",default='GEDI_2yr_BaselineSIM_beam10_v01.h5')
    clp.add_option("--modis",type="string",dest="modNamen",help="MODIS input file name",default='/gpfs/data1/vclgp/htang/GEDI/ISS/MODIS_ancillary_sets.tif')
    clp.add_option("--outRoot",type="string",dest="outRoot",help="output file name root",default='teast')
    clp.add_option("--power",type="int",dest="power",help="power beam or not",default=0)
    clp.add_option("--cloud",type="float",dest="cFrac",help="cloud fraction",default=0.0)
    clp.add_option("--seed",type="int",dest="seed",help="random number seed",default=0)
    clp.add_option("--bounds",type=float,nargs='*',dest="bounds",help="Bounds of interest in input EPSG",default=[-180,180,-1,1])
    clp.add_option("--iEPSG",type=int,dest="iEPSG",help="Input EPSG",default=4326)
    clp.add_option("--oEPSG",type=int,dest="oEPSG",help="Output EPSG",default=4326)
    return clp.parse_args()

## Main ##
if __name__ == '__main__':
  (options, args) = clParser()
  modNamen=options.modNamen
  inNamen=options.inNamen
  minX=options.bounds[0]
  maxX=options.bounds[1]
  minY=options.bounds[2]
  maxY=options.bounds[3]
  iEPSG=options.iEPSG
  oEPSG=options.oEPSG
  power=options.power
  cFrac=options.cFrac
  random.seed(options.seed)

  # resolution in degrees
  sRes=(60.0/(2*6374000*pi))*(2*pi)

  # read data
  f=h5py.File(inNamen,'r')
  y=np.array( f['/GT/latitude'])
  x=np.array(f['/GT/longitude'])
  sim_sun_el=np.array(f['/GT/sun_el'])
  delta_time=f['/delta_time'].value
  f.close()

  # tanslate bounds from ALS to orbit projection
  epsg="epsg:"+str(iEPSG)
  inProj=Proj(init=epsg)
  epsg="epsg:"+str(oEPSG)
  outProj=Proj(init=epsg)
  bX,bY=transform(outProj,inProj,bounds[0:2],bounds[2:4])

  # which coords to use
  buff_res=(0.066667 * np.sin(52 * np.pi /180.0))*4
  useInd=np.sort(np.where((y>=(bY[0]-buff_res))&(y<(bY[1]+buff_res))&(x>=(bX[0]-buff_res))&(x<(bX[1]+buff_res)))[0])

  # read VCF
  driver = gdal.GetDriverByName('GTiff')
  dataset_MOD = gdal.Open(modNamen)
  b_landcover = dataset_MOD.GetRasterBand(1)
  b_vcf = dataset_MOD.GetRasterBand(2)
  b_startgrowing = dataset_MOD.GetRasterBand(3)
  b_endgrowing = dataset_MOD.GetRasterBand(4)
  cols = dataset_MOD.RasterXSize
  rows = dataset_MOD.RasterYSize
  transMod = dataset_MOD.GetGeoTransform()
  
  yOrigin = transMod[3]
  pixelWidth = transMod[1]
  pixelHeight = transMod[5]

  #prepare raster data set
  data_lc = b_landcover.ReadAsArray(0, 0, cols, rows)
  data_vcf = b_vcf.ReadAsArray(0, 0, cols, rows)
  data_leafon = b_startgrowing.ReadAsArray(0, 0, cols, rows)
  data_leafoff = b_endgrowing.ReadAsArray(0, 0, cols, rows)

  # open output
  f=open(options.namen, 'w+')

  # mission 0 time
  t0=11*30.4375+1

  # loop over usable tracks
  for i in range(1,len(useInd)):
    ind=useInd[i]
    # only look at sections of lines. Reset counters as a line passes
    if(useInd[i]!=(useInd[i-1]+1)):
      tX=x[ind]
      tY=y[ind]
      lon=np.empty(shape=(0))
      lat=np.empty(shape=(0))
      waveID=[]
      numb=0
      continue

    # cloud cover
    if(random.random()>=cFrac):
      continue

    # loop over intersecting pixels
    crossing_lon,crossing_lat=intersectLines(x[ind],x[ind-1],y[ind],y[ind-1],xOrigin,yOrigin,pixelWidth,pixelHeight)

    # loop over intersected cells
    for j in range(0,len(crossing_lon)):
      # determine modis properties
      xInd=int((crossing_lon[j]-xOrigin)/pixelWidth)
      yInd=int((crossing_lat[j]-yOrigin)/pixelHeight)
      dX=crossing_lon[j+1]-crossing_lon[j]
      dY=crossing_lon[j+1]-crossing_lon[j]
      if((xInd<0)|(xInd>=cols)|(yInd<0)|(yInd>=rows)):
        continue

      vcf=data_vcf[yInd][xInd]
      lc=data_lc[yInd][xInd]
      onDat=data_leafon[yInd][xInd]
      offDat=data_leafoff[yInd][xInd]
      sub_dtime=delta_time[ind]

      # check is forest
      if((lc>0)&(lc<10)&(vcf>10)):
        # check is leaf-on
        c_doy=(sub_dtime/(24*60*60))+t0
        while( c_doy > 365.25 ):
          c_doy-=365.25

        if (lc==1)|(lc==2):  # evergreen
          condition_doy=True
        elif onDat < offDat: #mostly nothern hemi
          condition_doy=(c_doy>onDat)&(c_doy<offDat)
        else: #growing period extending from year1 to the next year
          condition_doy=(c_doy>onDat)|(c_doy<offDat)

        # how many footprints?
        nIn=int(sqrt((crossing_lon[j]-crossing_lon[j-1])**2+(crossing_lat[j]-crossing_lat[j-1])**2)/sRes)

        # is the beam usable?
        if((condition_doy==True)&((power==1)||(sim_sun_el[ind]>=83)||(vcf<70))):
          # footprint labels
          if(sim_sun_el[ind]>=83):
            tim="night"
          else:
            tim="day"

          ye=int(((delta_time[ind]/(24*60*60))+t0)/365.25)+18
          doy=int(((delta_time[ind]/(24*60*60))+t0)%365.25)

          # if use, loop along length
          for p in range(0,nIn+1):
            tX=tX+dX
            tY=tY+dY
            lon=np.append(lon,tX)
            lat=np.append(lat,tY)
            waveID.append("%d.%d.%d.%s.zen.%d.y.%d.doy.%d"%(power,j,numb,tim,int(sim_sun_el[ind]),ye,doy))
            numb=numb+1  # increment footprint counter along
        else:  # increment to the end
          numb=numb+nIn
          tX=tX+nIn*dX
          tY=tY+nIn*dY

    # translate
    epsg="epsg:"+str(iEPSG)
    inProj=Proj(init=epsg)
    epsg="epsg:"+str(oEPSG)
    outProj=Proj(init=epsg)
    uX,uY=transform(inProj,outProj,lon,lat)
    # write data
    for p in range(0,len(lon)):
      line=str(lon[p])+" "+str(lat[p])+" "+waveID[p]+"\n"
      x.write(line)

  f.close()

  print("Written to ",options.namen)

