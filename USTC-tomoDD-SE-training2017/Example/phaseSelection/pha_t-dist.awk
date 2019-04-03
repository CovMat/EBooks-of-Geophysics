#
# Output travel time and hypocenter distance.
#
#    Input Files: station.dat; phase.dat
#    Output_Files: t_dist.dat
#
# Usage: awk -f pha_t-dist.awk
#

BEGIN{

  while( (getline < "station.dat")>0 ){
       i = i+1;
       staID[i] = $1;
       sta_lat[i] = $2;
       sta_lon[i] = $3;
       sta_dep[i] = -$4/1000;
  }
  nsta = i;

  n_P = 0;
  n_S = 0;
  while( (getline < "phase.dat_sichuan")>0 ){
       # event Info
       if ($1 == "#"){
          eveID = $15;
          lat1 = $8;
          lon1 = $9;
          dep1 = $10;
       }
       # phase Info
       if ($1 != "#"){
          j = 0;
          for (i=1;i<=nsta;i++){
              if (staID[i]==$1){
                 lat2 = sta_lat[i];
                 lon2 = sta_lon[i];
                 dep2 = sta_dep[i];
                 j = j+1;
              }
          }
          if (j==0){
             print "No this station:",$1
          }
          if (j>1){
             print $1," : This station exist more than once, Please Check!!"
          }
          if (j==1){
              # Output travel time and hypocenter distance.
              dlat = lat1-lat2
              dlon = lon1-lon2
              dist = sqrt( (dlat*111)^2 + (dlon*(cos(lat1*3.1415926/180)*111))^2 + (dep1-dep2)^2 );
              if ($4=="P"){
                 n_P = n_P +1;
                 print $2, dist, "1"  > "t_dist.dat"
              }
              else if ($4=="S"){
                 n_S = n_S +1;
                 print $2, dist, "2"  > "t_dist.dat"
              }
          }
       }
  }
  print "Sum of P-wave times:",n_P
  print "Sum of S-wave times:",n_S
}
