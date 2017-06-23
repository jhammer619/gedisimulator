#!/bin/csh -f


set bin="$GEDIRAT_ROOT"


# defaults
@ findLat=1
set epsg=32732
@ seed=0
set orbitDir="/gpfs/data1/vclgp/data/gedi/ancillary/orbits/test"
set orbitAng=51

# read command line
while ($#argv>0)
  switch("$argv[1]")

  case -input
    set input="$argv[2]"
  shift argv;shift argv
  breaksw

  case -lat
    set meanLat=$argv[2]
    @ findLat=0
  shift argv;shift argv
  breaksw

  case -epsg
    set epsg=$argv[2]
  shift argv;shift argv
  breaksw

  case -orbitDir
    set orbitDir="$argv[2]"
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-input name;     input metric file"
    echo "-lat y;          latitude in degrees EPSG:4326"
    echo "-epsg n;         EPSG code"
    echo "-orbitDir dir;   directory with GEDI track density files"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Type listGediWaves.csh -help"
    exit;

  endsw
end


# do we need to find the mean latitude
if( $findLat )then
  set meanLat=`gawk -f $bin/meanCoord.awk < $input|gdaltransform -s_srs EPSG:$epsg -t_srs EPSG:4326|gawk '{print $2}'`
endif

# GEDI track angle at this point
set trackAngle=`echo $meanLat $orbitAng|gawk 'BEGIN{pi=4*atan2(1,1);scale=pi/180}{lat=$1*scale;a=$2*scale;print a*cos(lat*pi/(2*a))}'`

# read track density
set trackRoot=`echo $meanLat|gawk '{if($1>=-0.5)printf("_lat%dN_",int($1+0.5));else printf("_lat%dS_",int($1+0.5))}'`
set trackFile=`ls $orbitDir/*$trackRoot*.csv`
set trackNumbs=`gawk -F, -f $bin/readGEDItrack.awk` < $trackFile

# sample footprints at appropriate density


