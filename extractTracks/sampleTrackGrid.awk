BEGIN{
  readHist=1;
  numb=nHist=0;
}

($0){
  if($1=="###"){  # change reading mode
    readHist=0;  
    # lay out tracks
    dx=maxX-minX;
    dy=maxY-minY;
    area=dx*dy;
    

  }else if(readHist==0){  # read metrics

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
  }else{      # read istogram
    aInd=$1;
    dInd=$2;
    hist[aInd,dInd]=$3;
    if(aInd>nHist)nHist]aInd;
    if(dInd>nHist)nHist]dInd;
  }
}

END{
  # decide which saved points to output



}

