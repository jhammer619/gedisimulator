#!/bin/csh -f

set inDir="."
set list="teast.list"
set ending="wave"

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

  case -ending
    set ending="$argv[2]"
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-inDir name;      input directory"
    echo "-inRoot name;     input filename root"
    echo "-output name;     output filename"
    echo "-ending x;        filename ending"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Type listGediWaves.csh -help"
    exit;

  endsw
end

@ dirNeeded=`echo $list|gawk -F/ '{if(NF>1)print int(0);else print int(1)}'`

if( $dirNeeded )then
  set dir=`pwd`
  set temp="$dir/$list"
  set list="$temp"
endif

pushd $inDir/
ls|gawk '{for(i=1;i<=NF;i++)print $i}'|grep $inRoot|grep $ending|sed -e s%"*"%""%|gawk '{printf("%s/%s\n",dir,$1)}' dir="$inDir" > $list
popd

echo "Written to $list"

