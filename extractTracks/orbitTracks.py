import h5py
import argparse
import string
import numpy as np
from osgeo import gdal
import random
from pyproj import Proj, transform
from scipy import pi
from math import sqrt


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
      xDist=100000.0
    if(dy!=0.0):
      yDist=(yT-ty)/dy
    else:
      yDist=100000.0;

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
def getCmdArgs():
    '''
    Get commdnaline arguments
    GEDI grid resolution
    latitude band
    '''
    p = argparse.ArgumentParser(description=("Print out GEDI footprint locations from orbital traks\n"))
    p.add_argument("--modis", dest ="modNamen", type=str, default='/gpfs/data1/vclgp/htang/GEDI/ISS/MODIS_ancillary_sets.tif', help=("MODIS ancillary filename."))
    p.add_argument("--output",type=str,dest="outNamen",help="output filename",default='teast.coords')
    p.add_argument("--cloud",type=float,dest="cFrac",help="cloud fraction",default=0.0)
    p.add_argument("--seed",type=int,dest="seed",help="random number seed",default=0)
    p.add_argument("--bounds",type=float,nargs='*',dest="bounds",help="Bounds of interest in input EPSG",default=[-180,180,-1,1])
    p.add_argument("--iEPSG",type=int,dest="iEPSG",help="Input EPSG",default=4326)
    p.add_argument("--oEPSG",type=int,dest="oEPSG",help="Output EPSG",default=4326)
    p.add_argument("--mission_length",dest = "miss_len", type=int, default=395, help=("Total mission length in days since 01 Nov 2018 00:00:00 UTC.\nDefault for the full nominal 2-yr mission and a valid input should be less than 365*2."))
    p.add_argument("--orbit_sim_dir",dest = "orbit_sim_dir", type=str, default='/gpfs/data1/vclgp/htang/GEDI/ISS/BaselineSIM', help=("Directory containing orbital simulations from cott Luthcke.\nDefault /gpfs/data1/vclgp/htang/GEDI/ISS/BaselineSIM."))#added on 2nd Apr 2018
    p.add_argument("-s","--skip_track",dest = "skip_track_list", type=int, nargs='*', default=[-1], help=("Track numbering to be dropped in order to assess its impact on mission level 1 requirement"))#added on Feb 23-2018
    p.add_argument("-p","--power_beams",dest = "pow_beam_list", type=int, nargs='*', default=[5,6], help=("Track numbers of power beams"))
    p.add_argument("--usePhen",dest = "usePhen", type=int,default=0, help=("Use phenology sqitch"))
    cmdargs = p.parse_args()

    return cmdargs


## Main ##
if __name__ == '__main__':
  options=getCmdArgs()
  modNamen=options.modNamen
  orbit_sim_dir=options.orbit_sim_dir
  bounds=options.bounds
  iEPSG=options.iEPSG
  oEPSG=options.oEPSG
  pow_beam_list=options.pow_beam_list
  cFrac=options.cFrac
  random.seed(options.seed)
  miss_len=options.miss_len
  outNamen=options.outNamen
  skip_track_list=options.skip_track_list
  usePhen=options.usePhen

  # mission 0 time in days
  t0=11*30.4375+1

  # tanslate bounds from ALS to orbit projection
  epsg="epsg:"+str(iEPSG)
  inProj=Proj(init=epsg)
  epsg="epsg:"+str(oEPSG)
  outProj=Proj(init=epsg)
  bX,bY=transform(outProj,inProj,bounds[0:2],bounds[2:4])

  # resolution in degrees
  #sRes=transform(outProj,inProj,0,60)[1]
  sRes=60.0 #(60/6371007.181)*180/pi

  # read MODIS
  driver = gdal.GetDriverByName('GTiff')
  dataset_MOD = gdal.Open(modNamen)
  b_landcover = dataset_MOD.GetRasterBand(1)
  b_vcf = dataset_MOD.GetRasterBand(2)
  b_startgrowing = dataset_MOD.GetRasterBand(3)
  b_endgrowing = dataset_MOD.GetRasterBand(4)
  cols = dataset_MOD.RasterXSize
  rows = dataset_MOD.RasterYSize
  transMod = dataset_MOD.GetGeoTransform()
  xOrigin = transMod[0]
  yOrigin = transMod[3]
  pixelWidth = transMod[1]
  pixelHeight = transMod[5]
  #prepare raster data set
  data_lc = b_landcover.ReadAsArray(0, 0, cols, rows)
  data_vcf = b_vcf.ReadAsArray(0, 0, cols, rows)
  data_leafon = b_startgrowing.ReadAsArray(0, 0, cols, rows)
  data_leafoff = b_endgrowing.ReadAsArray(0, 0, cols, rows)

  # open output
  fOut=open(outNamen, 'w')

  # loop over tracks
  for beamid in range(1,11):
    # do we need this beam?
    if beamid in skip_track_list:
      print("Skipping track",beamid)
      continue

    # read orbital data
    inNamen="%s/GEDI_2yr_BaselineSIM_beam%d_v01.h5" % (orbit_sim_dir,beamid)
    f=h5py.File(inNamen,'r')
    y=np.array( f['/GT/latitude'])
    x=np.array(f['/GT/longitude'])
    sim_sun_el=np.array(f['/GT/sun_el'])
    delta_time=f['/delta_time'].value
    f.close()
 
    # which coords to use
    buff_res=(0.066667 * np.sin(52 * np.pi /180.0))*4
    useInd=np.sort(np.where((y>=(bY[0]-buff_res))&(y<(bY[1]+buff_res))&(x>=(bX[0]-buff_res))&(x<(bX[1]+buff_res)))[0])

    # blank arrays
    uX=np.empty(shape=(0))
    uY=np.empty(shape=(0))
    waveID=[]
    numb=0
 
    # loop over usable tracks
    for i in range(1,len(useInd)):
      ind=useInd[i]
      # if over mission length, skip
      if((delta_time[ind]/(24*60*60))>miss_len):
        break
 
      # only look at sections of lines. Reset counters as a line passes
      if(useInd[i]!=(useInd[i-1]+1)):
        numb=0
        continue

      # WHich MODIS pixels are intersected
      crossing_lon,crossing_lat=intersectLines(x[ind],x[ind-1],y[ind],y[ind-1],xOrigin,yOrigin,pixelWidth,pixelHeight)
 
      # loop over intersected MODIS pixels
      for j in range(0,len(crossing_lon)-1):
 
        # cloud cover
        if(random.random()<cFrac):
          print("Cloudy")
          continue
 
        # determine modis properties
        xInd=int((crossing_lon[j]-xOrigin)/pixelWidth)
        yInd=int((crossing_lat[j]-yOrigin)/pixelHeight)
        if((xInd<0)|(xInd>=cols)|(yInd<0)|(yInd>=rows)):
          print("No modis")
          continue

        vcf=data_vcf[yInd][xInd]
        lc=data_lc[yInd][xInd]
        onDat=data_leafon[yInd][xInd]
        offDat=data_leafoff[yInd][xInd]
        power =1 if(beamid in pow_beam_list) else 0
 

        # check is leaf-on
        c_doy=(delta_time[ind]/(24*60*60))+t0
        while( c_doy > 365.25 ):
          c_doy-=365.25
        if(lc==1)|(lc==2):  # evergreen
          condition_doy=True
        elif onDat < offDat: #mostly nothern hemi
          condition_doy=(c_doy>onDat)&(c_doy<offDat)
        else: #growing period extending from year1 to the next year
          condition_doy=(c_doy>onDat)|(c_doy<offDat)
 
        # is the beam usable?
        if(((usePhen==0)|(condition_doy==True))&((power==1)|(sim_sun_el[ind]>=83)|(vcf<70))):
          # footprint labels
          if(sim_sun_el[ind]>=83):
            tim="night"
          else:
            tim="day"

          # date
          ye=int(((delta_time[ind]/(24*60*60))+t0)/365.25)+18
          doy=int(((delta_time[ind]/(24*60*60))+t0)%365.25)
 
          # if use, loop along length, keeping 60 m spacing
          if(numb>0):
            tX,tY=transform(inProj,outProj,crossing_lon[j],crossing_lat[j])
            diff=sqrt((tX-uX[len(uX)-1])**2+(tY-uY[len(uY)-1])**2)
            nSteps=int(diff/sRes+0.5)
            x00=uX[len(uX)-1]+dX*(nSteps+1)*sRes
            y00=uY[len(uY)-1]+dY*(nSteps+1)*sRes
            print(diff,nSteps,tX,tY,uX[len(uX)-1],uY[len(uY)-1])
          else:
            x00,y00=transform(inProj,outProj,crossing_lon[j],crossing_lat[j])

          # how many footprints within this MODIS cell?
          xS,yS=transform(inProj,outProj,crossing_lon[j],crossing_lat[j])
          xE,yE=transform(inProj,outProj,crossing_lon[j+1],crossing_lat[j+1])
          nIn=int(sqrt((xE-xS)**2+(yE-yS)**2)/sRes)

          # direction of track for footprints
          dX=xE-xS
          dY=yE-yS
          if((dX>2000)|(dY>2000)):
            print("A step too far",dX,dY,crossing_lon[j],crossing_lat[j],crossing_lon[j+1],crossing_lat[j+1])
          totLen=sqrt(dX**2+dY**2)
          dX=dX/totLen
          dY=dY/totLen

          # loop over footprints within cell
          for p in range(0,nIn+1):
            uX=np.append(uX,x00+dX*p*sRes)
            uY=np.append(uY,y00+dY*p*sRes)
            waveID.append("%d.%d.%d.%s.zen.%d.y.%d.doy.%d"%(power,j,numb,tim,int(sim_sun_el[ind]),ye,doy))
            numb=numb+1  # increment footprint counter along
        else:  # if not usable, increment to the end
          print("Unusable",condition_doy,power,sim_sun_el[ind],vcf)
 
    # write data
    useCoord=np.where((uX>=bounds[0])&(uX<=bounds[1])&(uY>=bounds[2])&(uY<=bounds[3]))[0]
    for p in range(0,len(uX)): #useCoord:
      line=str(uX[p])+" "+str(uY[p])+" "+waveID[p]+"\n"
      fOut.write(line)

  fOut.close()

  print("Written to ",outNamen)

