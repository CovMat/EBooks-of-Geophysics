BEGIN{
  iter = 10;
  ny = 14;
  nz = 11;
  i = 0;
  start = 0;
}
{
  if(start==1 && i==iter){
     j=j+1;

     if(j>=1 && j<=ny*nz){
        print $0 >name;
     }
     else{
        start=0;
     }
  }

  if($1=="DWS" && $3=="P-wave."){
     start=1;
     i=i+1;
     name="DWS_P";
     j=0;
     # output the next ny*nz lines.
  }
  if($1=="DWS" && $3=="S-wave."){
     start=1;
     #i=i+1;
     name="DWS_S";
     j=0;
     # output the next ny*nz lines.
  }
}
END{

}
