# Reads Hao's orbital sim output
# c0 cov - day - descend
# c1 cov - day - ascend
# c2 cov - night - descend
# c3 cov - night - ascend
# c4 pow - day   - descend
# c5 pow - day   - ascend
# c6 pow - night - descend
# c7 pow - night - ascend

BEGIN{
  nCases=8;
  numb=0;
}

($1!="id"){
  for(i=0;i<nCases;i++){
    nIn[i][0]=nIn[i][1]=$(i+4);
  }
  numb++;
}

END{
  # turn data into a PDF of ascend and descend for now
  maxIn=0;
  for(i=0;i<numb;i++){
    n[i][0]=n[i][1]=0;
    for(j=0;j<nCases;j+=2){
      k=j%2;
      n[i][k]+=nIn[i][j];
    }
    if(n[i][0]>maxIn)maxIn=n[i][0];
    if(n[i][1]>maxIn)maxIn=n[i][1];
  }

  for(i=0;i<=maxIn;i++)hist[i][0]=hist[i][1]=0;
  total[0]=total[1]=0;
  for(i=0;i<numb;i++){
    for(j=0;j<2;j++){
      hist[n[i][j]][j]++;
      total[j]++;
    }
  }

  # write out PDF
  for(i=0;i<=maxIn;i++){
    hist[i][0]/=total[0];
    hist[i][1]/=total[1];
    print i,hist[i][0],hist[i][1];
  }
}

