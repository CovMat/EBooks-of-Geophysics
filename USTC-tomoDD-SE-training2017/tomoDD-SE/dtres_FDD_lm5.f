	subroutine dtres_FDD(log, ndt, stdim, nsrc, nsta, 
     &	dt_dt, dt_idx,
     &	dt_ista, dt_ic1, dt_ic2,
     &	src_cusp, src_t, tmp_ttp, tmp_tts,
     &	dt_cal, dt_res, eve_sta)

	use tomoFDD
	implicit none

c	Parameters:
	integer		log		! Log-file identifier
	integer		ndt		! No. of data
	integer		stdim		! Column dimenson of arrays tmp_tt[ps]
	integer		nsrc		! No. of sources
	integer         nsta
	real		dt_dt(MAXDATA)	! [1..ndt] Observed time differences
	integer		dt_idx(MAXDATA)	! [1..ndt]
	integer		dt_ista(MAXDATA)! [1..ndt] Station indices
	integer		dt_ic1(MAXDATA)	! [1..ndt] Event indices
	integer		dt_ic2(MAXDATA)	! [1..ndt] Event indices
	integer		src_cusp(MAXEVE)! [1..nsrc] Event keys
	real		src_t(MAXEVE)	! [1..nsrc] Event times
	real		tmp_ttp(MAXOBS,MAXEVE)! [1.., 1..MAXOBS]
	real		tmp_tts(MAXOBS,MAXEVE)! [1.., 1..MAXOBS]
	real		dt_cal(MAXDATA)	! [1..ndt] Theoretical time differences
	real		dt_res(MAXDATA)	! [1..ndt] Time-difference residuals

	integer         eve_sta(MAXEVE,MAXOBS+1)
	integer         evID1, evID2, j1, j2

c	Local variables:
	integer		i,j
	real		tt1, tt2

      write(log,'("~ getting residual vector...")')

      if (nsrc.eq.1) then
c        Single source
         do i=1,ndt
            dt_res(i) = dt_dt(i)
         enddo
      else
c        Mulitple sources
	 tt1 = 0.0
	 tt2 = 0.0
         do i=1,ndt

cz--- find the sequence number for each event
            evID1 = dt_ic1(i)
            evID2 = dt_ic2(i)

cz--- call the subroutine find_id to look for the sequence number
cz--- corresponding to each station-event pair
            call find_id(eve_sta, dt_ista(i), evID1, j1)
            call find_id(eve_sta, dt_ista(i), evID2, j2)
            if (dt_idx(i).eq.1 .or. dt_idx(i).eq.3) then
c              P phase
               tt1 = tmp_ttp(j1,dt_ic1(i)) -
     &              src_t(dt_ic1(i))/1000
               tt2 = tmp_ttp(j2,dt_ic2(i)) -
     &              src_t(dt_ic2(i))/1000
              if (tt1.eq.0 .or. tt2.eq.0) then
	       write(*,*)' P Phase'
	       write(*,*)  tmp_ttp(j1,dt_ic1(i)), tmp_ttp(j2,dt_ic2(i))
	       write(*,*)  dt_ic1(i), dt_ic2(i)
	       write(*,*)  dt_ista(i)
               write(*,'("FATAL ERROR (theor tt=0). ")')
c               stop
              endif
            elseif (dt_idx(i).eq.2 .or. dt_idx(i).eq.4) then
c              S phase
               tt1 = tmp_tts(j1,dt_ic1(i)) -
     &              src_t(dt_ic1(i))/1000
               tt2 = tmp_tts(j2,dt_ic2(i)) -
     &              src_t(dt_ic2(i))/1000
             if (tt1.eq.0 .or. tt2.eq.0) then
	       write(*,*)' S Phase'
	       write(*,*)  tmp_tts(j1,dt_ic1(i)), tmp_tts(j2,dt_ic2(i))
	       write(*,*)  dt_ic1(i), dt_ic2(i)
	       write(*,*)  dt_ista(i)
               write(*,'("FATAL ERROR (theor tt=0).")')
c               stop 	       
             endif

            endif
	    if (dt_ic1(i) .ne. dt_ic2(i) ) then ! difference time
                dt_cal(i) = tt1 - tt2
                dt_res(i) = dt_dt(i) - dt_cal(i)
            else
	        dt_res(i) = dt_dt(i) - tt1  ! absolute catalog time
	    endif 
         enddo
      endif

      end !of subroutine dtres


