#!/bin/csh -f

# directory containing awk scripts needed
set bin="$GEDIRAT_ROOT/extractTracks"

# defaults
@ findLat=1
set epsg=32732
@ seed=0
set orbitDir="/gpfs/data1/vclgp/data/gedi/ancillary/orbits/test"
set orbitAng=51
set cloudFrac=0.5
set latRes=5

# temporary workspace files
set tempTrack="/tmp/gediTrackSpace.$$.dat"
set tempMetric="/tmp/gediMetricCoord.$$.dat"


# read command line
while ($#argv>0)
  switch("$argv[1]")

  case -metricFile
    set metricFile="$argv[2]"
    @ readMetric=1
  shift argv;shift argv
  breaksw

  case -alsFile
    set alsFile="$argv[2]"
    @ readMetric=0
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

  case -seed
    @ seed=$argv[2]
  shift argv;shift argv
  breaksw

  case -cloud
    set cloudFrac="$argv[2]"
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-metricFile name;  input metric file, if sampling from metrics"
    echo "-alsFile name;     input ALS bound file, if setting coordinates"
    echo "-lat y;            latitude in degrees EPSG:4326"
    echo "-epsg n;           EPSG code of ALS data"
    echo "-orbitDir dir;     directory with GEDI track density files"
    echo "-cloud frac;       cloud fraction"
    echo "-seed n;           random number seed"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Type listGediWaves.csh -help"
    exit;

  endsw
end


# do we need to find the mean latitude
if( $findLat && $readMetric)then   # from metric fle
  set meanLat=`gawk -f $bin/meanCoord.awk < $metricFile|gdaltransform -s_srs EPSG:$epsg -t_srs EPSG:4326|gawk '{print $2}'`
else if( $findLat )then            # from ALS bounds file
  set meanLat=`gawk 'BEGIN{x=y=0;n=0}($0&&($1!-"#")){x+=$2+$5;y+=$3+$6;n+=2}END{print x/n,y/n}' < $alsFile|gdaltransform -s_srs EPSG:$epsg -t_srs EPSG:4326|gawk '{print $2}'`
endif

# GEDI track angle at this point
set trackAngle=`echo $meanLat $orbitAng|gawk 'BEGIN{pi=4*atan2(1,1);scale=pi/180}{lat=$1*scale;a=$2*scale;print a*cos(lat*pi/(2*a))}'`

# read track density
set trackRoot=`echo $meanLat|gawk -v res=$latRes '{lat=(int($1/res+0.5+100)-100)*res;if(lat>=0)printf("_lat%dN_",lat);else printf("_lat%dS_",-1*lat)}'`
set trackFile=`ls $orbitDir/*$trackRoot*.csv`
#gawk -F, -v cFrac=$cloudFrac -v seed=$seed -f $bin/readGEDItrack.awk < $trackFile > $workSpace

# set out GEDI tracks
if( $readMetric )then # sample footprints at appropriate density from metric grid
  # prefilter data
  gawk -F, '($1!="id"){x[0]=x[1]=0;for(j=4;j<=NF;j++){k=(j-4)%2;x[k]+=$j};print $3,x[0],y[0]}' < $trackFile > $tempTrack
  #gawk '($0){if($1!="#"){}else{for(i=1;i<=NF;i++){if($i=="lon,")lonInd=($i-1);else if($i=="lat,")latInd=$(i-1)}}}' < $metricFile > $tempMetric

  # choose metric lines
  $bin/sampleTrackGrid -trackFile $tempTrack -metricFile $metricFile -cloudFrac $cloudFrac -seed $seed
else                  # choose coordinates for GEDI tracks

endif

# tidy up workspace
if( -e $tempTrack )rm $tempTrack
if( -e $tempMetric )rm $tempMetric

echo "Written to $output"

