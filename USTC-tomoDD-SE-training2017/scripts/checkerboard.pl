#!/usr/bin/perl -w 
#####################################################################
# This is for the tomoDD to do the checkerboard test
# When using this, please change the default values for your own need
# for ploting, please change some gmt command values by yourself for your
# need.
#
# Xin Zhang
#####################################################################

use strict;
use File::Copy;

#========== default values

my $tomoDD = "tomoSPDR2.1"; # tomoDD
my $tomoDD_syn="tomoSPDR2.1_syn"; # tomoDD_syn
my $anomaly_x=1;
my $anomaly_y=1;
my $anomaly_z=1;
my $tomoDD_inp="tomoSPDR.inp"; # tomoDD.inp
my $tomoDD_syn_inp="tomoSPDR_syn.inp"; # tomoDD_syn.inp
my $anomaly=0.05;
my $chk_vp=[];
my $chk_vs=[];
my $rec_vp=[];
my $rec_vs=[];
my $ini_vp=[];
my $ini_vs=[];
my $recover=[];
my $vp_model="Vp_model.dat";
my $vs_model="Vs_model.dat";
my $iteration=1;
my @head;
my @xarray;
my @yarray;
my @zarray;
my $NX=0;
my $NY=0;
my $NZ=0;
my $syn_dirnm="Syn";
my $vel_dirnm="Vel";
################### end!

@ARGV>=1 or die "Usage: checkerboard.pl MOD";
for(my $i=1;$i<=$iteration;$i++){
##### the out loop
mkdir ($syn_dirnm,0777) unless -d $syn_dirnm;
mkdir ($vel_dirnm) unless -d $vel_dirnm;
system("cp MOD ./$vel_dirnm/MOD");

open MOD,"<MOD" or die "can not open the file:$!";
open MOD2,">./$syn_dirnm/MOD" or die "can not open the file:$!";

####read the header
@head=split /\s+/,<MOD>;
print MOD2 "@head\n";
@xarray=split /\s+/,<MOD>;
print MOD2 "@xarray\n";
@yarray=split /\s+/,<MOD>;
print MOD2 "@yarray\n";
@zarray=split /\s+/,<MOD>;
print MOD2 "@zarray\n";
$NX=$head[1];
$NY=$head[2];
$NZ=$head[3];
#### end read header

# read the vp and vp/vs 
open VELP,">./Vp.1d" or die "can not open the file:$!";
open VELS,">./Vs.1d" or die "can not open the file:$!";
for(my $i=0;$i<$NY*$NZ;$i++){
   my $tmp=<MOD>;
   print VELP $tmp;
   chomp $tmp;
   push @$chk_vp,[split /\s+/, $tmp];
   push @$ini_vp,[split /\s+/,$tmp];
}
for(my $i=0;$i<$NY*$NZ;$i++){
   my $tmp=<MOD>;
   chomp $tmp;
   push @$chk_vs,[split /\s+/,$tmp];
}
   #print "s velocity\n";
   #print $chk_vp->[1][1];
####### write the vs intial model
   for(my $j=0;$j<$NY*$NZ;$j++){
       for(my $i=0;$i<$NX;$i++){
           if($chk_vs->[$j][$i] == 0){
              print "error:s velocity is zero: $j $i \n";
           }
           $ini_vs->[$j][$i]=$chk_vp->[$j][$i]/$chk_vs->[$j][$i]; 
           printf VELS "%5.3f ",$chk_vp->[$j][$i]/$chk_vs->[$j][$i];
       }
       print VELS "\n";
   }
close VELP;
close VELS;
#print $ini_vp->[10][10];
#print $ini_vs->[10][10];
######## end!

&create_checkerboard($chk_vp,$anomaly);
#&create_checkerboard($chk_vs,-$anomaly);

###### output checkerboard vp and vs to MOD
for(my $i=0;$i<$NY*$NZ;$i++){
    for(my $j=0;$j<$NX;$j++){
        printf MOD2 "%5.3f ",$chk_vp->[$i][$j];
    }
    print MOD2 "\n";
}
for(my $i=0;$i<$NY*$NZ;$i++){
    for(my $j=0;$j<$NX;$j++){
        printf MOD2 "%5.3f ",$chk_vs->[$i][$j];
    }
    print MOD2 "\n";
}
close MOD;
close MOD2;
##### end!

#### synthetic data
copy("$tomoDD_syn_inp","./$syn_dirnm/") or warn "can not copy the file:$!";
copy("ak135.15.SKS","./$syn_dirnm/") or warn "can not copy the file:$!";
copy("layer-16.dat","./$syn_dirnm/") or warn "can not copy the file:$!";
copy("$tomoDD_syn","./$syn_dirnm/") or warn "can not copy the file:$!";
chdir $syn_dirnm;
system("chmod 777 $tomoDD_syn");
#system("pwd");
system ("./$tomoDD_syn  $tomoDD_syn_inp");
&abs2abs();
&dt2dt();
&dc2dc();
chdir "../";
#### end

#### recover checkerboard
system("pwd");
system("cp $tomoDD_inp ./$vel_dirnm/");
copy("ak135.15.SKS","./$vel_dirnm/") or warn "can not copy the file:$!";
copy("layer-16.dat","./$vel_dirnm/") or warn "can not copy the file:$!";
copy("$tomoDD","./$vel_dirnm/") or warn "can not copy the file:$!";
chdir $vel_dirnm;
#system("pwd");
system("chmod 777 $tomoDD");
system("./$tomoDD  $tomoDD_inp");
#### end recover
close MOD;
close MOD2;

open VPM,"<./$vp_model" or die "can not open the file:$!";
open VSM,"<./$vs_model" or die "can not open the file:$!";
for(my $i=0;$i<$NY*$NZ;$i++){
   my $tmp=<VPM>;
   chomp $tmp;
   $tmp=~s/^\s+//;
   push @$rec_vp,[split /\s+/, $tmp];
}
for(my $i=0;$i<$NY*$NZ;$i++){
   my $tmp=<VSM>;
   chomp $tmp;
   $tmp=~s/^\s+//;
   push @$rec_vs,[split /\s+/,$tmp];
}
close VPM;
close VSM;
chdir "../";
#### end read

##### end!
}

#### sub create checkerboard
sub create_checkerboard() {
   #==========some initial value
   my $chk_vp=$_[0];
   my $ano_z=$_[1];
   my $stepz=0;
   for(my $k=1;$k<$NZ-1;$k++){
       $stepz=$stepz+1;
       if($stepz==$anomaly_z){
          $ano_z=-$ano_z;
          $stepz=0;
       }
       my $stepy=0;
       my $ano_y=$ano_z;
       for(my $j=1;$j<$NY-1;$j++){
           $stepy=$stepy+1;
           if($stepy==$anomaly_y){
              $ano_y=-$ano_y;
              $stepy=0;
           }
           my $stepx=0;
           my  $ano_x=$ano_y;
           my $m=$k*$NY+$j;
           for(my $i=1;$i<$NX-1;$i++){
               $stepx=$stepx+1;
               $chk_vp->[$m][$i]=$chk_vp->[$m][$i]*(1+$ano_x);
               if($stepx==$anomaly_x){
                  $ano_x=-$ano_x;
                  $stepx=0;
               }
            }
        }
    }
}
######### end!

##### sub abs2abs
sub abs2abs(){
    open ABS,"<absolute.syn" or die "can not open the file:$!";
    open ABS2,">syn.absolute.dat" or die "can not open the file:$!";
    my $i=0;my $event_old;
    while(<ABS>){
          $i=$i+1;
          (my $event,my $sta,my $time,my $weight,my $phase)=split;
          if($i==1){
               $event_old=$event;
               print ABS2 "#  $event\n";
           }
          if($event!=$event_old){
               print ABS2 "#  $event\n";
               print ABS2 "$sta  $time  $weight  $phase\n";
               $event_old=$event;
           }
          else{
                print ABS2 "$sta  $time  $weight  $phase\n";
          }
      }
close ABS;
close ABS2;
}
##### end!
##### sub dt2dt
sub dt2dt(){
    open DT,"<dt.syn" or die "can not open the file:$!";
    open DT2,">syn.dt.ct" or die "can not open the file:$!";
    my $i=0;my $ev1_old;my $ev2_old;
    while(<DT>){
          $i=$i+1;
          (my $ev1,my $ev2,my $sta,my $t1,my $t2,my $qual,my $pha)=split;
          if($i==1){
               $ev1_old=$ev1;
               $ev2_old=$ev2;
               print DT2 "# $ev1 $ev2\n";   
           }
          if($ev1!=$ev1_old || $ev2!=$ev2_old){
               $ev1_old=$ev1;
               $ev2_old=$ev2;
               print DT2 "# $ev1 $ev2\n";
               print DT2 "$sta $t1 $t2 $qual $pha\n";
          }
          else{
               print DT2 "$sta $t1 $t2 $qual $pha\n";
          }
    }
   close DT;
   close DT2;
}
##### end!
##### sub dc2dc
sub dc2dc(){
    open DC,"<dtcc.syn" or die "can not open the file:$!";
    open DC2,">syn.dt.cc" or die "can not open the file:$!";
    my $i=0;my $ev1_old;my $ev2_old;
    while(<DC>){
          $i=$i+1;
          (my $ev1,my $ev2,my $sta,my $t1,my $t2,my $pha)=split;
          if($i==1){
               $ev1_old=$ev1;
               $ev2_old=$ev2;
               print DC2 "# $ev1 $ev2 0.0\n";
           }
          if($ev1!=$ev1_old || $ev2!=$ev2_old){
               $ev1_old=$ev1;
               $ev2_old=$ev2;
               print DC2 "# $ev1 $ev2 0.0\n";
               print DC2 "$sta $t1 $t2 $pha\n";
          }
          else{
               print DC2 "$sta $t1 $t2 $pha\n";
          }
    }
    close DC;
    close DC2;
}
###### end!

#### end
