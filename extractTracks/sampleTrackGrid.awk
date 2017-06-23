BEGIN{
  readHist=1;
  numb=nHist=0;
  srand(int(seed));
  cumul[-1]=0.0;



}

($0){
  if($1=="###"){  # change reading mode
    readHist=0;  
    # lay out tracks
    dx=maxX-minX;
    dy=maxY-minY;
    area=dx*dy;
    # number of tracks over area
    nX=int(dx/1000+1);
    totA=totD=0;
    for(i=0;i<nX;i++){
      n=rand();
      for(j=0;j<nHist;j++){{
        if((n>=cumul[j-1])&&(n<=cumul[j]))break;
      }
      totA+=nAsc[j];
      totD+=nDes[j];
    }
    # determine starting points
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
  }else{      # read histogram
    if(nHist>0)cumul[nHist]=$3+cumul[nHist-1];
    else       cumul[nHist]=0.0;
    nAsc[nHist]=$1;
    nDes[nHist]=$2;
    nHist++;
  }
}

END{
  # decide which saved points to output



}

