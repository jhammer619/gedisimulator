BEGIN{
  numb=0;
  tol=0.000001;
}

{
  numb=0;
  xS[numb]=$1;
  yS[numb]=$2;
  zen[numb]=$3;
  waveID[numb]=$4;
  numb++;
}

END{
  # write them out
  for(i=0;i<nnumb;i++){
    y=yS[i];
    d=0.0;
    j=0;
    while((y<=(maxY+tol))&&(y>=(minY-tol))){
      x=xS[i]+d*sin(zen[i]);
      y=yS[i]+d*cos(zen[i]);

      if((x>=minX)&&(x<=maxX)&&(y>=minY)&&(y<=maxY)){
        print x,y,waveID[i];
      }
      d+=alongTrack;
      j++;
    }
  }
}

