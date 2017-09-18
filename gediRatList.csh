#!/bin/csh -f

#########################
# Runs GEDI simulations #
# over a list of coords #
#########################


set temp="/tmp/boundFiles.$$.dat"

# default settings
set pBuff=0.5
set outRoot="teast"
set LVIS=" "
set pSigma=" "
set pFWHM=" "
set fSigma=" "
set ground=" "
set sideLobe=" "
set lobeAng=" "
set topHat=" "
set noNorm=" "
set checkCove=" "
set maxScanAng=" "
set polyGround=" "
set pFile=" "
set res=" "
set hdf=" "
@ maxPer=40000
set ending="wave"


# read options
while ($#argv>0)
  switch("$argv[1]")

  case -inList
    set inList="$argv[2]"
  shift argv;shift argv
  breaksw

  case -outRoot
    set outRoot="$argv[2]"
  shift argv;shift argv
  breaksw

  case -coordList
    set coordList="$argv[2]"
  shift argv;shift argv
  breaksw

  case -pBuff
    set pBuff="$argv[2]"
  shift argv;shift argv
  breaksw

  case -LVIS
    set LVIS="-LVIS"
  shift argv
  breaksw

  case -pSigma
    set pSigma="-pSigma $argv[2]"
  shift argv;shift argv
  breaksw

  case -pFWHM
    set pFWHM="-pFWHM $argv[2]"
  shift argv;shift argv
  breaksw

  case -fSigma
    set fSigma="-fSigma $argv[2]"
  shift argv;shift argv
  breaksw

  case -ground
    set ground="-ground"
  shift argv
  breaksw

  case -sideLobe
    set sideLobe="-sideLobe"
  shift argv
  breaksw

  case -lobeAng
    set lobeAng="-lobeAng $argv[2]"
  shift argv;shift argv
  breaksw

  case -topHat 
     set topHat="-topHat"
  shift argv
  breaksw

  case -noNorm
    set noNorm="-noNorm"
  shift argv
  breaksw

  case -checkCover
    set checkCover="-checkCover"
  shift argv
  breaksw

  case -maxScanAng
    set maxScanAng="-maxScanAng $argv[2]"
  shift argv;shift argv
  breaksw

  case -pFile
    set pFile="-readPulse $argv[2]"
  shift argv;shift argv
  breaksw

  case -res
    set res="-res $argv[2]"
  shift argv;shift argv
  breaksw

  case -polyGround
    set polyGround="-polyGround"
  shift argv
  breaksw

  case -hdf
    set hdf="-hdf"
  shift argv
  breaksw

  case -maxPer
    @ maxPer=$argv[2]
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-inList name;      name of list with las file names"
    echo "-outRoot root;     output filename root"
    echo "-coordList name;   file containing list of coordinates and waveIDs (x y waveID)"
    echo "-pBuff size;       RAM buffer size"
    echo "-LVIS;             use LVIS pulse size"
    echo "-pSigma sigma;     pulse width, sigma in m"
    echo "-pFWHM fwhm;       pulse FWHM in ns"
    echo "-pFile name;       read the pulse from an ASCII file"
    echo "-fSigma sigma;     footprint width, sigma in m"
    echo "-res res;          output range resolution, in metres"
    echo "-ground;           output ground and canopy waveforms"
    echo "-sideLobe;         use the old side-lobes"
    echo "-lobeAng ang;      side lobe major axis azimuth, degrees"
    echo "-topHat;           use top hat rather than Gaussian footprint"
    echo "-noNorm;           do not normalise for footprint density"
    echo "-checkCover;       check that at least 2/3 of footprint is covered by ALS"
    echo "-maxScanAng ang;   maximimum scan angle to use, degrees"
    echo "-polyGround;       find the ground by fitting polynomial"
    echo "-hdf;              output in HDF5"
    echo "-maxPer n;         maximum number of runs per processor"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type gediRatList.csh -help"
    exit;

  endsw
end

set grabDir="butabe"
if( ! -e $grabDir )mkdir $grabDir

@ nCoords=`wc -l` < $coordList

# split into sub files
set tempRoot="/tmp/gediRatList.$$"
@ nReps=`echo "$nCoords $maxPer"|gawk '{print int($1/$2+1)}'`
gawk -v maxPer=$maxPer -f $GEDIRAT_ROOT/awk/splitCoords.awk root="$tempRoot" < $coordList

@ j=0
while( $j <= $nReps )
  set input="$tempRoot.$j.coords"
  if( ! -e $input )continue

  set output="$outRoot.$j"
  set grab="$grabDir/$output.grab"

  if( ! -e $grab )then
    touch $grab
    overlapLasFiles.csh -input $inList -coordList $input -rad 100 -output $temp
    gediRat -inList $temp -output $output -listCoord $input -pBuff $pBuff -waveID $waveID $LVIS $pSigma $pFWHM $fSigma $ground $sideLobe $lobeAng $topHat $noNorm $checkCove $maxScanAng $pFile $res $polyGround $hdf
  endif

  # delete zero sized files
  if( -e $temp )rm $temp
  if( -e $input )rm $input
  @ j++
end

echo "Ping"

