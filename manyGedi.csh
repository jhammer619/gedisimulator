#!/bin/csh -f

set temp="/tmp/boundFiles.$$.dat"
set grabDir="butabe"
set pBuff=2


while ($#argv>0)
  switch("$argv[1]")

  case -coordList
    set coordList="$argv[2]"
  shift argv;shift argv
  breaksw

  case -inList
    set inList="$argv[2]"
  shift argv;shift argv
  breaksw

  case -outRoot
    set outRoot="$argv[2]"
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-inList name;        input las file list"
    echo "-outRoot name;       output filename root"
    echo "-coordList name;     wave coordinate list"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type echi.rat -help"
    exit;

  endsw
end

if( ! -e $grabDir )mkdir $grabDir
@ nCoords=`wc -l` < $coordList

@ i=1
while( $i <= $nCoords )
  set coord=`gawk -v i=$i '{if(NR==i){for(j=1;j<=NF;j++)print $j}}'` < $coordList
  set x=$coord[1]
  set y=$coord[2]
  set waveID=$coord[3]
 
  set output="$outRoot.$waveID.wave"
  set grab="$grabDir/$output:r.grab"

  if( ! -e $grab )then
    touch $grab
    overlapLasFiles.csh -input $inList -coord $x $y -rad 100 -output $temp
    gediRat -inList $temp -output $output -ground -checkCover -coord $x $y -pBuff $pBuff -waveID $waveID
  endif

  if( -e $temp )rm $temp

  @ i++
end

echo "Ping"

