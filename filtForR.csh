#!/bin/csh -f

set bin="$HOME/src/gediRat"
@ noName=1

while ($#argv>0)
  switch("$argv[1]")

  case -input
    set input="$argv[2]"
  shift argv;shift argv
  breaksw

  case -output
    set output="$argv[2]"
    @ noName=0
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-input name;    input filename"
    echo "-output name;   output filename"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type echi.rat -help"
    exit;

  endsw
end

if( $noName )then
  set output="$input:r.csv"
endif

gawk -f $bin/filtHeadForR.awk < $input |sed -e s%,%_%g -e s%" "%""%g|sed -e s%_%","%g > $output
gawk -f $bin/filtDataForR.awk < $input >> $output
echo "Written to $output"

