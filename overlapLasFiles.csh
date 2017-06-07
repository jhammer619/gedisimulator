#!/bin/csh -f

set bin="$GEDIRAT_ROOT"

@ useCircle=0
set rad=80
set x=0
set y=0
set minX=-100000000
set maxX=100000000
set minY=-100000000
set maxY=100000000
set output="test.list"


while ($#argv>0)
  switch("$argv[1]")

  case -input
    set input="$argv[2]"
  shift argv;shift argv
  breaksw

  case -output
    set output="$argv[2]"
  shift argv;shift argv
  breaksw

  case -bounds
    set minX="$argv[2]"
    set minY="$argv[3]"
    set maxX="$argv[4]"
    set maxY="$argv[5]"
  shift argv;shift argv;shift argv;shift argv;shift argv
  breaksw

  case -rad
    set rad="$argv[2]"
  shift argv;shift argv
  breaksw

  case -coord
    @ useCircle=1
    set x="$argv[2]"
    set y="$argv[3]"
  shift argv;shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-input name;                   input filename"
    echo "-output name;                  output filename"
    echo "-bounds minX minY maxX maxY;   bounds of interest for square"
    echo "-coord x y;                    point of interest if footprint"
    echo "-rad r;                        radius of footprint"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type echi.rat -help"
    exit;

  endsw
end

if( $useCircle )
  set minX=`echo $x $rad|gawk '{printf("%.4f",$1-$2)}'`
  set minY=`echo $y $rad|gawk '{printf("%.4f",$1-$2)}'`
  set maxX=`echo $x $rad|gawk '{printf("%.4f",$1+$2)}'`
  set maxY=`echo $y $rad|gawk '{printf("%.4f",$1+$2)}'`
else
  set minX=`echo $minX $rad|gawk '{printf("%.4f",$1-$2)}'`
  set minY=`echo $minY $rad|gawk '{printf("%.4f",$1-$2)}'`
  set maxX=`echo $maxX $rad|gawk '{printf("%.4f",$1+$2)}'`
  set maxY=`echo $maxY $rad|gawk '{printf("%.4f",$1+$2)}'`
endif

gawk -f $bin/overlapLasFiles.awk -v minX=$minX -v maxX=$maxX -v minY=$minY -v maxY=$maxY < $input > $output
echo "Written to $output"

