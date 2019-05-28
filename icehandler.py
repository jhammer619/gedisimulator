
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
    inProj=Proj(init="epsg:"+str(self.epsg))
    outProj=Proj(init="epsg:"+str(4326))
    lon,lat=transform(inProj,outProj,self.x,self.y)
    segment_dist_x=np.remainder(np.array(self.dists/20))
    segment_id=np.around(np.array(self.dists/20))
    # open output
    f=h5py.File(outNamen,'w')
    # make top level groups
    f.create_group('#ref#')
    f.create_group('GT1L')
    f['GT1L'].create_group('bckgrd_ATLAS')
    f['GT1L'].create_group('geolocation')
    f['GT1L'].create_group('heights')
    f['GT1L'].create_group('orbit_info')
    f['GT1L'].create_group('veg_truth')
    # populate data
    f['#ref#'].create_dataset('a',data=[0,0])
    f['GT1L']['geolocation'].create_dataset('segment_dist_x',data=segment_dist_x)
    f['GT1L']['geolocation'].create_dataset('segment_id',data=segment_id)
    f['GT1L']['geolocation'].create_dataset('sigma_h',data=np.full(numb,0.4))
    f['GT1L']['geolocation'].create_dataset('surf_type',data=np.full(numb,1))
    f['GT1L']['geolocation'].create_dataset('segment_length',data=np.full(numb,20))
    f['GT1L']['heights'].create_dataset('delta_time',data=delta_time)
    f['GT1L']['heights'].create_dataset('h_ph',data=self.z)
    f['GT1L']['heights'].create_dataset('lon_ph',data=lon)
    f['GT1L']['heights'].create_dataset('lat_ph',data=lat)
    f['GT1L']['heights'].create_dataset('signal_conf_photon',data=self.signal)
    f['GT1L']['orbit_info'].create_dataset('rgt',data=(1))
    f['GT1L']['veg_truth'].create_dataset('x_utm',data=self.x)
    f['GT1L']['veg_truth'].create_dataset('y_utm',data=self.y)
    f['GT1L']['veg_truth'].create_dataset('refDEM',data=self.refdem)
    f['GT1L']['veg_truth'].create_dataset('rh50',data=self.RH50)
    f['GT1L']['veg_truth'].create_dataset('rh60',data=self.RH60)
    f['GT1L']['veg_truth'].create_dataset('rh75',data=self.RH75)
    f['GT1L']['veg_truth'].create_dataset('rh90',data=self.RH90)
    f['GT1L']['veg_truth'].create_dataset('rh95',data=self.RH95)
    f['GT1L']['veg_truth'].create_dataset('canopyz',data=self.CanopyZ)
    f['GT1L']['veg_truth'].create_dataset('canopy_cover',data=self.canopycover)
    # close up
    f.close()
    print("Written to",outNamen)
    return


  #################################

  def setDists(self):
    '''Set distances'''
    self.minX=np.min(self.x)
    self.minY=np.min(self.y)
    self.dists=np.sqrt((self.x-self.minX)**2+(self.y-self.minY)**2)

# iceSim() class end
########################################


def readCommands():
  '''Read the command line'''
  p = argparse.ArgumentParser(description=("Convert ICESat-2 .pts sims to HDF5"))
  p.add_argument("--input",dest="inNamen",type=str,help=("Input filename"))
  p.add_argument("--output",dest="outNamen",type=str,default='ice2.h5',help=("Output filename"))
  p.add_argument("--epsg",dest="epsg",type=int,default=32632,help=("Input EPSG"))
  cmdargs = p.parse_args()
  return cmdargs


########################################
# Main block

if __name__ == '__main__':
  # read commands
  cmdargs=readCommands()
  # read data
  data=iceSim(cmdargs.inNamen,cmdargs.epsg)
  # write data


# The end
########################################

