BEGIN{
  numb=n=0;
}

($0&&($1!="#")){
  namen=sprintf("%s.%d.coords",root,n);
  print $0 >> namen

  if((numb%maxPer)==0)n++;
  numb++;
}

