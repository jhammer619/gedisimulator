#!/bin/csh -f

while ($#argv>0)
  switch("$argv[1]")

  case -inDir
    set inDir="$argv[2]"
  shift argv;shift argv
  breaksw

  case -inRoot
    set inRoot="$argv[2]"
  shift argv;shift argv
  breaksw

  case -output
    set output="$argv[2]"
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-inDir name;     input directory"
    echo "-inRoot name;    input filename root"
    echo "-output name;    output filename"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type echi.rat -help"
    exit;

  endsw
end


ls $inDir/|gawk '{for(i=1;i<=NF;i++)print $i}'|gawk -F/ '{printf("%s/%s\n",dir,$NF)}' dir="$inDir"|grep $inRoot|grep wave|sed -e s%"*"%""% > $list
echo "Written to $list"


