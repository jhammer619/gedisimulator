BEGIN{
  readTrack=1;
  numb=nTrack=0;
  srand(int(seed));
  cumul[-1]=0.0;



}

($0){
  if($1=="###"){  # change reading mode
    readTrack=0;  
    # lay out tracks
    dx=maxX-minX;
    dy=maxY-minY;
    area=dx*dy;
    # number of tracks over area
    nX=int(dx/1000+1);
    totA=totD=0;


    # pick a track at random for the first cell
    tInd=int(rand()*nTrack-0.0000000001);
print nCross[tInd,0],nCross[tInd,1],nX;


    if(nX<=1){   # there is only a single cell

    }else{

    }

  }else if(readTrack==0){  # read metrics

    if($1!="#"){  # read data
      x=$lonInd;
      y=$latInd;
      # decide whether to record
      

    }else{  # read header
      for(i=1;i<=NF;i++){
        if($i=="lon,")lonInd=$(i-1);
        else if($i=="lat,")latInd=$(i-1);
      }
    }
  }else{      # read orbital tracks
    if($1!="id"){
      nCross[nTrack,0]=nCross[nTrack,1]=0;
      for(j=4;j<=NF;j++){
        k=(j-4)%2;
        nCross[nTrack,k]+=$j;
      }
      nTrack++;
    }
  }
}

END{
  # decide which saved points to output



}

