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
set pFile=" "
set res=" "


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
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type gediRatList.csh -help"
    exit;

  endsw
end


@ nCoords=`wc -l` < $coordList

@ i=1
while( $i <= $nCoords )
  set coord=`gawk -v i=$i '{if(NR==i){for(j=1;j<=NF;j++)print $j}}'` < $coordList
  set x=$coord[1]
  set y=$coord[2]
  set waveID=$coord[3]

  set output="$outRoot.$waveID.wave"

  if( ! -e $output )then
    touch $output
    overlapLasFiles.csh -input $inList -coord $x $y -rad 100 -output $temp
    gediRat -inList $temp -output $output -coord $x $y -pBuff $pBuff -waveID $waveID $LVIS $pSigma $pFWHM $fSigma $ground $sideLobe $lobeAng $topHat $noNorm $checkCove $maxScanAng $pFile $res
  endif

  if( -e $temp )rm $temp

  @ i++
end

echo "Ping"

