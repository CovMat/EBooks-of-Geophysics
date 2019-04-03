#
# Output travel time and hypocenter distance after excluding the outliers.
#
#    Input Files: station.dat; phase.dat
#    Input Parameters: slope_P, b_P, b1_P, b2_P;
#                      slope_S, b_S, b1_S, b2_S;
#          Please determine these values after running t_dist.m
#
#    Output_Files: t_dist.dat
#
# Usage: awk -f pha_t-dist_selection.awk
#

BEGIN{

  # Parameters
  slope_P = 0.151221;
  b_P = 3.291462;
  b1_P = 7;
  b2_P = 7;

  slope_S = 0.276095;
  b_S = 2.841760;
  b1_S = 7;
  b2_S = 7;

  # read station info
  while( (getline < "station.dat")>0 ){
       i = i+1;
       staID[i] = $1;
       sta_lat[i] = $2;
       sta_lon[i] = $3;
       sta_dep[i] = -$4/1000;
  }
  nsta = i;

  # read phase info
  n_P = 0;
  n_S = 0;
  while( (getline < "phase.dat_sichuan")>0 ){
       # event Info
       if ($1 == "#"){
          eveID = $15;
          lat1 = $8;
          lon1 = $9;
          dep1 = $10;
          n_eve_sta = 0;
          eveInfo = $0;
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

              # travel time and hypocenter distance.
              dlat = lat1-lat2
              dlon = lon1-lon2
              dist = sqrt( (dlat*111)^2 + (dlon*(cos(lat1*3.1415926/180)*111))^2 + (dep1-dep2)^2 );

              # Output the final phase data after exlcuding the outliers.
              if ($4=="P" &&  ($2 <= slope_P*dist + b_P + b1_P) && ($2 >= slope_P*dist + b_P - b2_P) ){
                 n_eve_sta = n_eve_sta + 1;
                 n_P = n_P +1;
                 if (n_eve_sta == 1){
                    print eveInfo > "phase.dat_selection"
                 }
                 print $0 > "phase.dat_selection"
              }
              else if ($4=="S" &&  ($2 <= slope_S*dist + b_S + b1_S) && ($2 >= slope_S*dist + b_S - b2_S) ){
                 n_eve_sta = n_eve_sta + 1;
                 n_S = n_S +1;
                 if (n_eve_sta == 1){
                    print eveInfo > "phase.dat_selection"
                 }
                 print $0 > "phase.dat_selection"
              }
          }
       }
  }
  print "Sum of final P-wave times:",n_P
  print "Sum of final S-wave times:",n_S
}
