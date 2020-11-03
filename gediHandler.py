
'''
Script to access simuated GEDI data
'''

##################################
import numpy as np
import h5py
from sys import exit
import matplotlib.pyplot as plt
if __name__ == '__main__':
  import argparse

###################################

class gediData(object):
  '''
  Simulated GEDI data handler
  '''

  def __init__(self,filename=None,minX=-100000000,maxX=100000000,minY=-1000000000,maxY=100000000,tocopy=None):
    '''
    Class initialiser. Calls a function
    to read waveforms between bounds
    '''
    if(filename):  # then read a file
      f=h5py.File(filename,'r')
      
      # check file format
      if(list(f)[0]!='BEAMDENSE'):   # then is a real file
        self.real=1
      else:                          # then is a simulated file
        self.real=0
      f.close()
      # read data
      if(self.real==1):
        self.readGEDI(filename,minX,maxX,minY,maxY)
      else:
        self.nWaves,self.lon,self.lat,self.waveID,self.wave,self.gWave,self.ZN,self.Z0,self.nBins,self.pSigma,self.fSigma,self.nTypes,self.idLen,self.slope,self.ZG,self.bDense,self.pDense,self.nPbins,self.zen=self.readSimGEDI(filename,minX,maxX,minY,maxY)
    else:          # create a blank space
      self.nWaves=0
      self.lon=None
      self.lat=None
      self.waveID=None
      self.wave=None
      self.gWave=None
      self.ZN=None
      self.Z0=None
      self.nBins=None
      self.nPbins=None
      self.pSigma=None
      self.fSigma=None
      self.nTypes=None
      self.idLen=None
      self.slope=None
      self.ZG=None
      self.bDense=None
      self.pDense=None
      self.zen=None


  ###########################################

  def readGEDI(self,filename,minX,maxX,minY,maxY):
    '''
    Read real GEDI data from file
    '''
    # open file for reading
    f=h5py.File(filename,'r')
    self.beamList=['BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011']
    self.nWaves=0




    # loop over beams
    for b in self.beamList:
      if((b in list(f))==False): # does this exist?
        continue                 # if not, skip it

      # read the coords and determine output
      allLat=(np.array(f[b]['geolocation']['latitude_bin0'])+np.array(f[b]['geolocation']['latitude_lastbin']))/2.0
      allLon=(np.array(f[b]['geolocation']['longitude_bin0'])+np.array(f[b]['geolocation']['longitude_lastbin']))/2.0
      useInd=np.where((allLat>=minY)&(allLat<=maxY)&(allLon>=minX)&(allLon<=maxX))

      if(len(useInd[0])>0):
        useInd=useInd[0]
      else:      # none in here
        continue

      # read lat lon
      self.lat=allLat[useInd]
      self.lon=allLon[useInd]

      # read Z
      self.Z0=np.array(f[b]['geolocation']['elevation_bin0'])[useInd]
      self.ZN=np.array(f[b]['geolocation']['elevation_lastbin'])[useInd]

      # read ID
      self.waveID=np.array(f[b]['shot_number'])[useInd]

      # read waveforms
      startInds=np.array(f[b]['rx_sample_start_index'])[useInd]
      self.lenInds=np.array(f[b]['rx_sample_count'])[useInd]
      totBins=np.sum(self.lenInds)
      self.waveform=np.empty(totBins,dtype=np.float32)

      # read raw data and repack
      self.nBins=np.max(self.lenInds)
      self.nWaves=len(useInd)
      jimlad=np.array(f[b]['rxwaveform'])
      self.wave=np.full((self.nWaves,self.nBins),np.median(jimlad),dtype=np.float32)
      lastInd=0
      for i in range(0,startInds.shape[0]):
        self.wave[i][0:self.lenInds[i]]=jimlad[startInds[i]:startInds[i]+self.lenInds[i]]
        lastInd=lastInd+self.lenInds[i]


      ## now read each element
      #for i in useInd:
      #  idx=int(f[b]['rx_sample_start_index'][i])-1
      #  cnt=nBins[i]
      #  rxdn=f[b]['rxwaveform'][idx:int(idx+cnt)]
      #  lon=f[b]["geolocation"]["longitude_bin0"][i]
      #  lat=f[b]["geolocation"]["latitude_bin0"][i]
      #  tot=f[b]['rx_sample_count'][i]
      #  if(np.max(rxdn)>300):
      #    #print(b,cnt,idx,'diff',tot-np.sum(rxdn),'stats',np.median(rxdn),np.min(rxdn),np.max(rxdn))
      #    plt.plot(rxdn)
      #    plt.ylabel('DN')
      #    plt.xlabel('Elevation (m)')
      #    outNamen=outRoot+"."+b+"."+str(i)+".png"
      #    plt.savefig(outNamen)
      #    plt.close()
      #    plt.clf()
      #    print("Written to",outNamen,lon,lat)
      ##self.nBins=int(len(f[b]['rxwaveform'])/nWaves)
      ##self.nWaves=self.nWaves+nWaves
    f.close()
    return

  ###########################################

  def readSimGEDI(self,filename,minX,maxX,minY,maxY):
    '''
    Read simulated GEDI data from file
    '''
    # open file for reading
    f=h5py.File(filename,'r')
    # extract region of interest
    lon=np.array(f['LON0'])
    lat=np.array(f['LAT0'])
    useInd=np.where((lon>=minX)&(lon<=maxX)&(lat>=minY)&(lat<=maxY))
    # if there are usable, read
    if(len(useInd)>0):
      useInd=np.ndarray.tolist(useInd[0])
      nWaves=len(useInd)
      lon=lon[useInd]
      lat=lat[useInd]
      # read data
      temp=np.array(f['WAVEID'])[useInd]
      # join up waveID characters
      if(temp.dtype!='int64'):
        waveID=[]
        for i in range(0,nWaves):
          waveID.append(''.join(np.array(temp[i], dtype=np.str)))
      else:
        waveID=temp
      # read all other data
      wave=np.array(f['RXWAVECOUNT'])[useInd]
      gWave=np.array(f['GRWAVECOUNT'])[useInd]
      ZN=np.array(f['ZN'])[useInd]
      Z0=np.array(f['Z0'])[useInd]
      nBins=np.array(f['NBINS'])[0]
      nPbins=np.array(f['NPBINS'])[0]
      pSigma=np.array(f['PSIGMA'])[0]
      fSigma=np.array(f['FSIGMA'])[0]
      nTypes=np.array(f['NTYPEWAVES'])[0]
      idLen=np.array(f['IDLENGTH'])[0]
      slope=np.array(f['SLOPE'])
      ZG=np.array(f['ZG'])
      bDense=np.array(f['BEAMDENSE'])
      pDense=np.array(f['POINTDENSE'])
      zen=np.array(f['INCIDENTANGLE'])
    else:
      nWaves=0
      lon=None
      lat=None
      waveID=None
      wave=None
      gWave=None
      ZN=None
      Z0=None
      nBins=None
      nPbins=None
      pSigma=None
      fSigma=None
      nTypes=None
      idLen=None
      slope=None
      ZG=None
      bDense=None
      pDense=None
      zen=None
    f.close()
    return(nWaves,lon,lat,waveID,wave,gWave,ZN,Z0,nBins,pSigma,fSigma,nTypes,idLen,slope,ZG,bDense,pDense,nPbins,zen)


  ###########################################

  def appendGEDI(self,tocopy,useInd=[]):
    '''Append another file to this one'''

    if(useInd==[]):  # copy all
      useInd=range(0,len(tocopy.lon))
    elif(len(useInd)==0): # none to copy
      return

    if(self.nWaves>0):  # if appending to existing data
      self.nWaves=self.nWaves+len(useInd)
      self.lon=np.append(self.lon,tocopy.lon)
      self.lat=np.append(self.lat,tocopy.lat)
      self.waveID=np.append(self.waveID,tocopy.waveID)
      self.wave=np.append(self.wave,tocopy.wave,axis=0)
      self.gWave=np.append(self.gWave,tocopy.gWave,axis=0)
      self.ZN=np.append(self.ZN,tocopy.ZN)
      self.Z0=np.append(self.Z0,tocopy.Z0)
      self.slope=np.append(self.slope,tocopy.slope)
      self.ZG=np.append(self.ZG,tocopy.ZG)
      self.bDense=np.append(self.bDense,tocopy.bDense)
      self.pDense=np.append(self.pDense,tocopy.pDense)
      self.zen=np.append(self.zen,tocopy.zen)
      # check for bin mismatch
      if(self.nBins!=tocopy.nBins):
        print("Bin number mismatch")
        exit(1)
    else:                  # if new data
      self.nWaves=tocopy.nWaves
      self.lon=tocopy.lon
      self.lat=tocopy.lat
      self.waveID=tocopy.waveID
      self.wave=tocopy.wave
      self.gWave=tocopy.gWave
      self.ZN=tocopy.ZN
      self.Z0=tocopy.Z0
      self.nBins=tocopy.nBins
      self.pSigma=tocopy.pSigma
      self.fSigma=tocopy.fSigma
      self.nTypes=tocopy.nTypes
      self.idLen=tocopy.idLen
      self.slope=tocopy.slope
      self.ZG=tocopy.ZG
      self.bDense=tocopy.bDense
      self.pDense=tocopy.pDense
      self.nPbins=tocopy.nPbins
      self.zen=tocopy.zen


  ###########################################

  def setOneZ(self,i):
    '''Set a single z array'''
    res=(self.Z0[i]-self.ZN[i])/self.nBins
    z=np.arange(self.Z0[i],self.ZN[i],-1*res)
    return(z)


  ###########################################

  def plotWaves(self,outRoot='teast',useInd=[]):
    '''Plot waveforms'''

    if(self.real==1):
      self.plotRealWaves(outRoot,useInd)
    else:
      self.plotSimWaves(outRoot,useInd)
    return


  ###########################################

  def plotRealWaves(self,outRoot,useInd):
    '''Plot waveforms from a real GEDI L1B file'''
    if(useInd==[]):
      useInd=range(0,len(self.lon))
    # loop over waves
    for i in useInd:
      # make z profile
      self.nBins=self.lenInds[i]
      self.res=(self.Z0[i]-self.ZN[i])/self.nBins
      self.z=np.linspace(self.Z0[i],self.ZN[i],num=self.nBins)

      # determine noise for scaling ground return
      reflScale,meanN,stdev=self.meanNoise(i)
      # find bounds
      minX,maxX=self.findBounds(meanN,stdev,i)
      # plot it
      #plt.plot(self.wave[i],self.z,label='Waveform')
      #plt.plot(self.gWave[i]*reflScale+meanN,z,label='Ground')
      plt.fill_betweenx(self.z,self.wave[i][0:self.nBins],meanN)
      # determine the x limit
      temp=np.sort(np.copy(self.wave[i][0:self.nBins]))
      minY=temp[int(temp.shape[0]*0.01)]-5  # to avoid artefacts
      plt.xlim(left=minY)
      #plt.legend()
      plt.ylim((minX,maxX))
      #plt.xlabel('DN')
      plt.ylabel('Elevation (m)')
      outNamen=outRoot+"."+str(self.waveID[i])+".x."+str(self.lon[i])+".y."+str(self.lat[i])+".png"
      plt.savefig(outNamen)
      plt.close()
      plt.clf()
      print("Written to",outNamen)
    return


  ###########################################

  def plotSimWaves(self,outRoot,useInd):
    '''Plot waveforms from a simulated GEDI file'''
    if(useInd==[]):
      useInd=range(0,len(self.lon))
    # loop over waves
    for i in useInd:
      # make z profile
      self.res=(self.Z0[i]-self.ZN[i])/self.nBins
      self.z=np.linspace(self.Z0[i],self.ZN[i],num=self.nBins)

      # determine noise for scaling ground return
      reflScale,meanN,stdev=self.meanNoise(i)
      # find bounds
      minX,maxX=self.findBounds(meanN,stdev,i)
      # plot it
      #plt.plot(self.wave[i],self.z,label='Waveform')
      #plt.plot(self.gWave[i]*reflScale+meanN,z,label='Ground')
      plt.fill_betweenx(self.z,self.wave[i],meanN)
      #plt.legend()
      #plt.xlim(left=0)
      plt.ylim((minX,maxX))
      #plt.xlabel('DN')
      plt.ylabel('Elevation (m)')
      outNamen=outRoot+"."+str(self.waveID[i])+".x."+str(self.lon[i])+".y."+str(self.lat[i])+".png"
      plt.savefig(outNamen)
      plt.close()
      plt.clf()
      print("Written to",outNamen)


  ###########################################

  def findBounds(self,meanN,stdev,i):
    '''Find the signal start and end'''
    thresh=3.5*stdev+meanN
    minWidth=3
    binList=np.where(self.wave[i]>thresh)
    buff=15

    topBin=0
    for j in range(0,len(binList[0])):
      if (binList[0][j]==(binList[0][j-1]+1))&(binList[0][j]==(binList[0][j-2]+2)):
        topBin=binList[0][j]
        break

    botBin=binList[len(binList)-1]
    for j in range(len(binList[0])-1,0,-1):
      if (binList[0][j]==(binList[0][j-1]+1))&(binList[0][j]==(binList[0][j-2]+2)):
        botBin=binList[0][j]
        break

    return(self.z[botBin]-buff,self.z[topBin]+buff)

  ###########################################

  def meanNoise(self,i):
    '''Calculate noise statistics'''
    statsLen=15
    noiseBins=int(statsLen/self.res)
    meanN=np.mean(self.wave[i][0:noiseBins])
    stdev=np.std(self.wave[i][0:noiseBins])
    totE=np.sum(self.wave[i]-meanN)*self.res
    return(totE,meanN,stdev)
 

  ###########################################

  def writeCoords(self):
    for i in range(0,len(self.lon)):
      print(self.lon[i],self.lat[i])


# end of gediData class
###########################################


###########################################
# read the command line

if __name__ == '__main__':
  def gediCommands():
    '''
    Read commandline arguments
    '''
    p = argparse.ArgumentParser(description=("Writes out properties of GEDI waveform files"))
    p.add_argument("--input",dest="inName",type=str,help=("Input GEDI HDF5 filename"))
    p.add_argument("--bounds", dest ="bounds", type=float,nargs=4,default=[-100000000,-100000000,100000000000,10000000000], help=("Bounds to plot between. minX minY maxX maxY"))
    p.add_argument("--outRoot",dest="outRoot",type=str,default='test',help=("Output graph filename root"))
    p.add_argument("--writeCoords",dest="writeCoords", action='store_true', default=False, help=("Write out coordinates insteda of plotting waveforms"))
    cmdargs = p.parse_args()
    return cmdargs


###########################################
# the main block

if __name__ == '__main__':
  # read the command line
  cmdargs=gediCommands()
  inName=cmdargs.inName
  bounds=cmdargs.bounds
  outRoot=cmdargs.outRoot

  # read data
  gedi=gediData(filename=inName,minX=bounds[0],maxX=bounds[2],minY=bounds[1],maxY=bounds[3])

  # mode switch
  if(cmdargs.writeCoords):
    gedi.writeCoords()
  else:
    print("Read",gedi.nWaves,"waveforms")
    # plot data
    gedi.plotWaves(outRoot=outRoot)

