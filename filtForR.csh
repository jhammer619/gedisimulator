#!/bin/csh -f

set bin="$HOME/src/gediRat"

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

gawk -f $bin/filtHeadForR.awk < $input > $output
gawk '{if($1!="#")print $0}' < $input|sed -e s%" "%","% >> $output
echo "Written to $output"

