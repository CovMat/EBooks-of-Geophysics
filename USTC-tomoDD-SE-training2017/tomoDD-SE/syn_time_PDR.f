	subroutine syn_time(log, ndt, stdim, nsrc, nsta, 
     &	dt_dt, dt_idx, dt_qual, eve_sta,
     &	dt_ista, dt_ic1, dt_ic2,dt_sta,
     &	src_cusp, src_t, tmp_ttp, tmp_tts)

	implicit none

	include'tomoFDD.inc'

c	Parameters:
	integer		log		! Log-file identifier
	integer		ndt		! No. of data
	integer		stdim		! Column dimenson of arrays tmp_tt[ps]
	integer		nsrc		! No. of sources
	integer         nsta
	real		dt_dt(MAXDATA)	! [1..ndt] Observed time differences
	integer		dt_idx(MAXDATA)	! [1..ndt]
	integer		dt_ista(MAXDATA)! [1..ndt] Station indices
        character       dt_sta(MAXDATA)*7
	integer		dt_ic1(MAXDATA)	! [1..ndt] Event indices
	integer		dt_ic2(MAXDATA)	! [1..ndt] Event indices
	integer		src_cusp(MAXEVE)! [1..nsrc] Event keys
	real		src_t(MAXEVE)	! [1..nsrc] Event times
	real		tmp_ttp(MAXOBS,MAXEVE)! [1.., 1..MAXOBS]
	real		tmp_tts(MAXOBS,MAXEVE)! [1.., 1..MAXOBS]
	real		dt_cal(MAXDATA)	! [1..ndt] Theoretical time differences
	real		dt_res(MAXDATA)	! [1..ndt] Time-difference residuals
	real            dt_qual(MAXDATA)

	integer         eve_sta(MAXEVE,MAXOBS+1)
	integer         evID1, evID2, j1, j2

c	Local variables:
	integer		i,j

        integer         log1, log2
        real            tt1, tt2
        real            uniform_noise(MAXDATA)
        integer*4       k
        real            rand
        real            a,b,c
        real            tt1_ns, tt2_ns
        real            qual1, qual2, noise

	write(log,'("~ getting synthetic data files ...")')

        a=0.12 ! Control the constant noise level
        b=0.03 ! Control the random noise level for P wave
        c=0.05 ! Control the random noise level for S wave

        a=0.0 ! Control the constant noise level
        b=0.0 ! Control the random noise level for P wave
        c=0.0 ! Control the random noise level for S wave

        call freeunit(log1)
        open(log1, file='dt.syn', status='unknown')
        call freeunit(log2)
        open(log2, file='absolute.syn', status='unknown')
                
c       construct uniform noise for each station
        do i=1, nsta
           uniform_noise(i)=(rand(k)-0.5)*a
        enddo
c       Mulitple sources
        tt1 = 0.0
        tt2 = 0.0
        k=0

c       Mulitple sources
	do i=1,ndt

cz--- find the sequence number for each event
	   evID1 = dt_ic1(i)
	   evID2 = dt_ic2(i)

cz--- call the subroutine find_id to look for the sequence number
cz--- corresponding to each station-event pair
	   call find_id(eve_sta, dt_ista(i), evID1, j1)
	   call find_id(eve_sta, dt_ista(i), evID2, j2)
	   if (dt_idx(i).eq.1 .or. dt_idx(i).eq.3) then
c       P phase
	      tt1 = tmp_ttp(j1,dt_ic1(i))
	      tt2 = tmp_ttp(j2,dt_ic2(i)) 
c       add constant noise to each station and random noise
c       to all the arrival times

              noise=(rand(k)-0.5)*b
              tt1=tt1+uniform_noise(dt_ista(i))+noise

              tt1_ns=uniform_noise(dt_ista(i))+noise
              if(tt1_ns.le.0.008) then
                 qual1=1
              elseif(tt1_ns.le.0.016) then
                 qual1=0.5
              elseif(tt1_ns.le.0.032) then
                 qual1=0.25
              elseif(tt1_ns.le.0.050) then
                 qual1=0.125
              else
                 qual1=0.065
              endif

              noise=(rand(k)-0.5)*b
              tt2=tt2+uniform_noise(dt_ista(i))+noise
	      
              tt2_ns=uniform_noise(dt_ista(i))+noise
              if(tt2_ns.le.0.008) then
                 qual2=1
              elseif(tt2_ns.le.0.016) then
                 qual2=0.5
              elseif(tt2_ns.le.0.032) then
                 qual2=0.25
              elseif(tt2_ns.le.0.050) then
                 qual2=0.125
              else
                 qual2=0.065
              endif

              if(dt_ic1(i) .ne. dt_ic2(i) ) then ! difference time               
                 write(log1,*)src_cusp(dt_ic1(i)),src_cusp(dt_ic2(i)),' ',dt_sta(i),' ',
     &                       tt1,tt2,(qual1+qual2)/2.0,' ','P'
              else
                 write(log2,*)src_cusp(dt_ic1(i)),' ',dt_sta(i),' ',tt1,
     &                       qual1,' ','P'
              endif

	   elseif (dt_idx(i).eq.2 .or. dt_idx(i).eq.4) then
c       S phase
	      tt1 = tmp_tts(j1,dt_ic1(i))
	      tt2 = tmp_tts(j2,dt_ic2(i))

c       add uniform noise to each station and
c       Gaussion noise to all the arrival time
              noise=(rand(k)-0.5)*c
              tt1=tt1+uniform_noise(dt_ista(i))+noise
	      
              tt1_ns=uniform_noise(dt_ista(i))+noise
              if(tt1_ns.le.0.008) then
                 qual1=1
              elseif(tt1_ns.le.0.016) then
                 qual1=0.5
              elseif(tt1_ns.le.0.032) then
                 qual1=0.25
              elseif(tt1_ns.le.0.050) then
                 qual1=0.125
              else
                 qual1=0.065
              endif

              noise=(rand(k)-0.5)*c
              tt2=tt2+uniform_noise(dt_ista(i))+noise
              tt2_ns=uniform_noise(dt_ista(i))+noise
              if(tt2_ns.le.0.008) then
                 qual2=1
              elseif(tt2_ns.le.0.016) then
                 qual2=0.5
              elseif(tt2_ns.le.0.032) then
                 qual2=0.25
              elseif(tt2_ns.le.0.050) then
                 qual2=0.125
              else
                 qual2=0.065
              endif

              if(dt_ic1(i) .ne. dt_ic2(i) ) then ! difference time
                 write(log1,*)src_cusp(dt_ic1(i)),src_cusp(dt_ic2(i)),' ',dt_sta(i),' ',
     &                        tt1,tt2,(qual1+qual2)/2.0,' ','S'
              else
                 write(log2,*)src_cusp(dt_ic1(i)),' ',dt_sta(i),' ',tt1,
     &                        qual1,' ','S'
              endif

	   endif
	   
	enddo
	
	end			!of subroutine dtres


