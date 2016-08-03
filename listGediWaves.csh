#!/bin/csh -f

set inDir="."
set list="teast.list"

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
    set list="$argv[2]"
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-inDir name;      input directory"
    echo "-inRoot name;     input filename root"
    echo "-output name;     output filename"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Type listGediWaves.csh -help"
    exit;

  endsw
end


pushd $inDir/
ls|gawk '{for(i=1;i<=NF;i++)print $i}'|grep $inRoot|grep wave|sed -e s%"*"%""%|gawk '{printf("%s/%s\n",dir,$1)}' dir="$inDir" > $list
popd

echo "Written to $list"

