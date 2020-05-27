
'''
Handles ICESat-2 simulations
'''


#################
# Packages

import numpy as np
import h5py
from pyproj import Proj, transform
if __name__ == '__main__':
  import argparse



########################################

class ice2(object):
  '''
  Handles real ICESat-2 data
  '''

  #################################

  def __init__(self,namen,epsg=4326,minX=-100000000,maxX=100000000,minY=-1000000000,maxY=100000000):
    '''Class initialiser'''
    self.readPhotons(namen,epsg=epsg,minX=minX,maxX=maxX,minY=minY,maxY=maxY)

  #################################

  def readPhotons(self,namen,epsg=4326,minX=-100000000,maxX=100000000,minY=-1000000000,maxY=100000000):
    '''Read ICESat-2 HDF5 file'''
    f=h5py.File(namen,'r')
    lon=np.array(f['gt1l']['heights']['lon_ph'])
    lat=np.array(f['gt1l']['heights']['lat_ph'])
    z=np.array(f['gt1l']['heights']['h_ph'])
    # reproject
    if(epsg!=4326):
      inProj=Proj(init="epsg:4326")
      outProj=Proj(init="epsg:"+str(epsg))
      x,y=transform(inProj, outProj, lon, lat)
    else:
      x=lon
      y=lat
    # filter if needed
    useInds=np.where((x>=minX)&(x<=maxX)&(y>=minY)&(y<=maxY))
    if(len(useInds)>0):
      useInds=useInds[0]
      self.x=x[useInds]
      self.y=y[useInds]
      self.z=z[useInds]

  #################################

  def writeCoords(self,outNamen="test.pts"):
    '''Write out coordinates and photons'''
    f=open(outNamen,'w')
    for i in range(0,len(self.x)):
      line=str(self.x[i])+" "+str(self.y[i])+" "+str(self.z[i])+"\n"
      f.write(line)
    f.close()
    print("Written to",outNamen)


########################################

class iceSim(object):
  '''
  Reads and acts upon
  an ICEsat-2 simulation
  '''

  #################################

  def __init__(self,namen,epsg):
    '''Class initialiser'''
    temp=np.loadtxt(namen,unpack=True, dtype=float,comments='#',delimiter=' ')
    self.x=temp[0]
    self.y=temp[1]
    self.z=temp[2]
    self.minht=temp[3]
    self.WFGroundZ=temp[4]
    self.RH50=temp[5]
    self.RH60=temp[6]
    self.RH75=temp[7]
    self.RH90=temp[8]
    self.RH95=temp[9]
    self.CanopyZ=temp[10]
    self.canopycover=temp[11]
    self.shotN=temp[12]
    self.photonN=temp[13]
    self.iterationN=temp[14]
    self.refdem=temp[15]
    self.noiseInt=temp[16]
    self.signal=np.array(temp[17],dtype=np.int16)
    self.epsg=epsg
    return

  #################################

  def writeHDF(self,outNamen):
    '''Write the output HDF5 file'''
    # convert some parameters
    numb=self.x.shape[0]
    self.setDists()
    delta_time=self.dists*0.7/7599.68
    self.setPhtCount()
    inProj=Proj(init="epsg:"+str(self.epsg))
    outProj=Proj(init="epsg:"+str(4326))
    lon,lat=transform(inProj,outProj,self.x,self.y)
    segment_dist_x=np.remainder(self.dists,20)
    segment_id=np.around(np.array(self.dists/20))
    # open output
    f=h5py.File(outNamen,'w')
    # make top level groups
    f.create_group('#ref#')
    f.create_group('gt1l')
    f['gt1l'].create_group('bckgrd_atlas')
    f['gt1l'].create_group('geolocation')
    f['gt1l'].create_group('heights')
    f['gt1l'].create_group('orbit_info')
    f['gt1l'].create_group('veg_truth')
    # populate data
    f['#ref#'].create_dataset('a',data=[0,0])
    f['gt1l']['bckgrd_atlas'].create_dataset('bckgrd_rate',data=np.full(numb,1))
    f['gt1l']['geolocation'].create_dataset('segment_dist_x',data=segment_dist_x)
    f['gt1l']['geolocation'].create_dataset('segment_id',data=segment_id)
    f['gt1l']['geolocation'].create_dataset('sigma_h',data=np.full(numb,0.4))
    f['gt1l']['geolocation'].create_dataset('surf_type',data=np.full(numb,1))
    f['gt1l']['geolocation'].create_dataset('segment_length',data=np.full(numb,20))
    f['gt1l']['geolocation'].create_dataset('ph_index_beg',data=np.full(numb,1))
    f['gt1l']['geolocation'].create_dataset('segment_ph_cnt',data=self.seg_phtcount)
    f['gt1l']['geolocation'].create_dataset('solar_elevation',data=np.full(numb,20))
    f['gt1l']['geolocation'].create_dataset('delta_time',data=delta_time)
    f['gt1l']['heights'].create_dataset('delta_time',data=delta_time)
    f['gt1l']['heights'].create_dataset('h_ph',data=self.z)
    f['gt1l']['heights'].create_dataset('lon_ph',data=lon)
    f['gt1l']['heights'].create_dataset('lat_ph',data=lat)
    f['gt1l']['heights'].create_dataset('signal_conf_photon',data=self.signal)
    f['gt1l']['heights'].create_dataset('dist_ph_across',data=np.full(numb,5))
    f['gt1l']['heights'].create_dataset('dist_ph_along',data=np.full(numb,10))
    f['gt1l']['heights'].create_dataset('pce_mframe_cnt',data=np.full(numb,1))
    f['gt1l']['orbit_info'].create_dataset('rgt',data=(1))
    f['gt1l']['orbit_info'].create_dataset('cycle_number',data=np.full(numb,1))
    f['gt1l']['veg_truth'].create_dataset('x_utm',data=self.x)
    f['gt1l']['veg_truth'].create_dataset('y_utm',data=self.y)
    f['gt1l']['veg_truth'].create_dataset('refDEM',data=self.refdem)
    f['gt1l']['veg_truth'].create_dataset('rh50',data=self.RH50)
    f['gt1l']['veg_truth'].create_dataset('rh60',data=self.RH60)
    f['gt1l']['veg_truth'].create_dataset('rh75',data=self.RH75)
    f['gt1l']['veg_truth'].create_dataset('rh90',data=self.RH90)
    f['gt1l']['veg_truth'].create_dataset('rh95',data=self.RH95)
    f['gt1l']['veg_truth'].create_dataset('canopyz',data=self.CanopyZ)
    f['gt1l']['veg_truth'].create_dataset('canopy_cover',data=self.canopycover)

    # close up
    f.close()
    print("Written to",outNamen)
    return


  #################################

  def setPhtCount(self):
    '''Set photon count rates per 20 m'''
    minD=np.min(self.dists)
    maxD=np.max(self.dists)
    self.seg_phtcount=np.histogram(self.dists,int((maxD-minD)/20.0+1))[0]

  #################################

  def setDists(self):
    '''Set distances'''
    self.minX=np.min(self.x)
    self.minY=np.min(self.y)
    self.dists=np.sqrt((self.x-self.minX)**2+(self.y-self.minY)**2)

  #################################

  def padData(self,minLen):
    '''Pad data to make a minimum transect length'''

    # footprint spacing
    step=0.7

    # determine track length
    dx=self.x[-1]-self.x[0]
    dy=self.y[-1]-self.y[0]
    trackLen=sqrt(dx**2+dy**2)

    # see if 
    if(trackLen<minLen):  # then we need to pad
      vectX=dx/trackLen
      vectY=dy/trackLen

      nPer=self.x.shape[0]
      nExtras=int(int(minLen/step)/nPer+1)

      # make copies of originals
      x=np.copy(self.x)
      y=np.copy(self.y)
      z=np.copy(self.z)
      minht=np.copy(self.minht)
      WFGroundZ=np.copy(self.WFGroundZ)
      RH50=np.copy(self.RH50)
      RH60=np.copy(self.RH60)
      RH75=np.copy(self.RH75)
      RH90=np.copy(self.RH90)
      RH95=np.copy(self.RH95)
      CanopyZ=np.copy(self.CanopyZ)
      canopycover=np.copy(self.canopycover)
      shotN=np.copy(self.shotN)
      photonN=np.copy(self.photonN)
      iterationN=np.copy(self.iterationN)
      refdem=np.copy(self.refdem)
      noiseInt=np.copy(self.noiseInt)
      signal=np.copy(self.signal)

      for i in range(1,nExtras):
        # are we reversing?
        isOdd=i%2

        if(isOdd):
          self.x=
          self.y=
          self.z=

          self.minht=np.append(self.minht,minht)
          self.WFGroundZ=np.append(self.
          self.RH50=np.append(self.
          self.RH60=np.append(self.
          self.RH75=np.append(self.
          self.RH90=np.append(self.
          self.RH95=np.append(self.
          self.CanopyZ=np.append(self.
          self.canopycover=np.append(self.
          self.shotN=np.append(self.
          self.photonN=np.append(self.
          self.iterationN=np.append(self.
          self.refdem=np.append(self.
          self.noiseInt=np.append(self.
          self.signal=np.append(self.
        else:
 

    return


# iceSim() class end
########################################


def readCommands():
  '''Read the command line'''
  p = argparse.ArgumentParser(description=("Convert ICESat-2 .pts sims to HDF5"))
  p.add_argument("--input",dest="inNamen",type=str,help=("Input filename"))
  p.add_argument("--output",dest="outNamen",type=str,default='ice2.h5',help=("Output filename"))
  p.add_argument("--epsg",dest="epsg",type=int,default=32632,help=("Input EPSG"))
  p.add_argument("--minLen",dest="minLen",type=float,default=0,help=("Minimum acceptable length"))
  cmdargs = p.parse_args()
  return cmdargs


########################################
# pad ICESat-2 data

########################################
# Main block

if __name__ == '__main__':
  # read commands
  cmdargs=readCommands()
  # read data
  data=iceSim(cmdargs.inNamen,cmdargs.epsg)
  # pad data if needed
  data.padData(cmd.minLen)
  # write data
  data.writeHDF(cmdargs.outNamen)

# The end
########################################

