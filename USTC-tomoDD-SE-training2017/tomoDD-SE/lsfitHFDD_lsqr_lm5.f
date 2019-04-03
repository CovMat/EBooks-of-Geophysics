	subroutine lsfitH_FDD_shot(log, iter, ndt, nev, nsrc,
     &	damp, eve_sta,
     &	idata, ev_cusp, src_cusp,
     &	dt_res, dt_wt,
     &	dt_ista, dt_ic1, dt_ic2, src_type,
     &	src_dx, src_dy, src_dz, src_dt, src_ex, src_ey, src_ez, src_et,
     &	exav, eyav, ezav, etav, dxav, dyav, dzav, dtav,
     &	rms_cc, rms_ct, rms_cc0, rms_ct0,
     &	rms_ccold, rms_ctold, rms_cc0old, rms_ct0old,
     &	tmp_xp, tmp_yp, tmp_zp, tmp_xs, tmp_ys, tmp_zs,dt_idx, acond)

c	implicit none
	
	use tomoFDD

c	Parameters:
	integer		log		! Log-file identifier
	integer		iter
	integer		ndt
	integer		nev
	integer		nsrc
	real		damp
	integer		idata
	integer		ev_cusp(MAXEVE)	! [1..nev] Event keys
	integer		src_cusp(MAXEVE)! [1..nev] Event keys
	real		dt_res(MAXDATA)	! [1..ndt]
	real		dt_wt(MAXDATA)	! [1..ndt]
	integer		dt_ista(MAXDATA)! [1..ndt]
	integer		dt_ic1(MAXDATA)	! [1..ndt]
	integer		dt_ic2(MAXDATA)	! [1..ndt]
	real		src_dx(MAXEVE)	! [1..nev]
	real		src_dy(MAXEVE)	! [1..nev]
	real		src_dz(MAXEVE)	! [1..nev]
	real		src_dt(MAXEVE)	! [1..nev]
	real		src_ex(MAXEVE)	! [1..nev]
	real		src_ey(MAXEVE)	! [1..nev]
	real		src_ez(MAXEVE)	! [1..nev]
	real		src_et(MAXEVE)	! [1..nev]
	integer         src_type(MAXEVE) ! [1..nev]
	real		exav
	real		eyav
	real		ezav
	real		etav
	real		dxav
	real		dyav
	real		dzav
	real		dtav
	real		rms_cc
	real		rms_ct
	real		rms_cc0
	real		rms_ct0
	real		tmp_xp(MAXOBS,MAXEVE)! [1.., 1..nev]
	real		tmp_yp(MAXOBS,MAXEVE)! [1.., 1..nev]
	real		tmp_zp(MAXOBS,MAXEVE)! [1.., 1..nev]
	real		tmp_xs(MAXOBS,MAXEVE)! [1.., 1..nev]
	real		tmp_ys(MAXOBS,MAXEVE)! [1.., 1..nev]
	real		tmp_zs(MAXOBS,MAXEVE)! [1.., 1..nev]
	integer		dt_idx(MAXDATA)	! [1..ndt]
	real		acond		! Condition number
        integer         eve_sta(MAXEVE,MAXOBS+1)

c	Local variables:
	real		anorm
	real		arnorm
	real		atol
	real		btol
	real		conlim
	character	dattim*25
	real		d(MAXDATA+4)	! Data vector
	real		dtavold
	real		dxavold
	real		dyavold
	real		dzavold
	real		etavold
	real		exavold
	real		eyavold
	real		ezavold
	real		factor
	integer		i, i1,i2, j
	integer		istop
	integer		itnlim
	integer		iw(2*(8*MAXDATA+4*MAXEVE)+1)	! lsqr index array
	integer		k1
	integer		k2
	integer		leniw
	integer		lenrw
	integer		m
	integer		n
	integer		nar
	integer		nndt
	real		norm(MAXEVE*4)
	real            norm1(MAXEVE*4)
	real		norm_test(MAXEVE*4)
	real		resvar1
	real		rms_cc0old
	real		rms_ccold
	real		rms_ct0old
	real		rms_ctold
	real		rnorm
	real		rw(8*MAXDATA+4*MAXEVE)
	real		se(MAXEVE*4)	! Solution error
	real		w1(MAXEVE*4)	! Work space
	real		w2(MAXEVE*4)	! Work space
	real		wtinv(MAXDATA+4)! +4 = mean shift constr
	real		wt(MAXDATA+4)
	real		x(MAXEVE*4)	! Solution vector
	real		xnorm
	integer 	unknowns
	integer         old_index(MAXEVE*4)
	integer         new_index(MAXEVE*4)
	integer         old_data_index(MAXDATA)
	integer         new_data_index(MAXDATA)
	integer         fix(MAXEVE*4)
	integer         new_loc
        


	write(log,'(/,"~ setting up G matrix.. ")')



	j=0 ! find # of earthquake data	
	do i=1,ndt
c---    remove blasts and shots, only relocate earthquakes

	   new_data_index(i)=0
	   old_data_index(i)=0
	   
	   if(src_type(dt_ic1(i)).eq.0.and.src_type(dt_ic2(i)).eq.0) then 
                                                      !Earthquake data
	      j=j+1
	   endif
	enddo
	
	nndt=j
c	write(*,*)'ndt=',ndt,'nndt=',nndt
	
c       If mean shift not contstrained
c       nar = 8*ndt
c       nndt = ndt
	
	nar = 8*nndt
	
	iw(1) = nar
	
c       Prepare sparse data and design vector
	j=0
	do i=1,ndt
c       Weight data first
	   wt(i) = dt_wt(i)
	   if (abs(wt(i)).gt.3) write(16,*)'Warning! Weight=',wt(i)
	   
c       Set up non-zero G matrix elements and apply weights
	   if (nsrc.eq.1) then
	      k1 = 1
	      k2 = 1
	   else
              call find_id(eve_sta,dt_ista(i),dt_ic1(i),k1)
              call find_id(eve_sta,dt_ista(i),dt_ic2(i),k2)
	   endif

	  
	   if(src_type(dt_ic1(i)).eq.0.and.src_type(dt_ic2(i)).eq.0) then 
                                                        ! earthquake data
	      j=j+1
	      
	      old_data_index(j)=i
	      new_data_index(i)=j
	      d(j) = dt_res(i)*1000.0 * wt(i)
	      if (wt(i).ne.0) then
		 wtinv(j) = 1.0/wt(i)
	      else
		 wtinv(j) = 1.0
	      endif	      

	      iw(1+j) = j
	      iw(1+nndt+j) = j
	      iw(1+2*nndt+j) = j
	      iw(1+3*nndt+j) = j
	      iw(1+4*nndt+j) = j
	      iw(1+5*nndt+j) = j
	      iw(1+6*nndt+j) = j
	      iw(1+7*nndt+j) = j
	      
	      
	      if (dt_idx(i).eq.2 .or. dt_idx(i).eq.4) then ! S-wave data
		 rw(j)        = tmp_xs(k1,dt_ic1(i)) * wt(i) 
		 rw(nndt+j)   = tmp_ys(k1,dt_ic1(i)) * wt(i) 
		 rw(2*nndt+j) = tmp_zs(k1,dt_ic1(i)) * wt(i) 
		 rw(3*nndt+j) = wt(i)
                 !write(*,*)i, wt(i),  rw(j), rw(nndt+j), rw(3*nndt+j)
		 if (dt_ic1(i).ne.dt_ic2(i)) then 
		    rw(4*nndt+j) = -tmp_xs(k2,dt_ic2(i)) * wt(i) 
		    rw(5*nndt+j) = -tmp_ys(k2,dt_ic2(i)) * wt(i) 
		    rw(6*nndt+j) = -tmp_zs(k2,dt_ic2(i)) * wt(i) 
		    rw(7*nndt+j) = -wt(i)
		 else
		    rw(4*nndt+j) = 0.0
		    rw(5*nndt+j) = 0.0
		    rw(6*nndt+j) = 0.0
		    rw(7*nndt+j) = 0.0	       
		 endif
	      else
		 rw(j)        = tmp_xp(k1,dt_ic1(i)) * wt(i)
		 rw(nndt+j)   = tmp_yp(k1,dt_ic1(i)) * wt(i)
		 rw(2*nndt+j) = tmp_zp(k1,dt_ic1(i)) * wt(i)
		 rw(3*nndt+j) = wt(i)
                 !write(*,*)i, wt(i),  rw(j), rw(nndt+j), rw(2*nndt+j)
		 if(dt_ic1(i).ne.dt_ic2(i)) then
		    rw(4*nndt+j) = -tmp_xp(k2,dt_ic2(i)) * wt(i)
		    rw(5*nndt+j) = -tmp_yp(k2,dt_ic2(i)) * wt(i)
		    rw(6*nndt+j) = -tmp_zp(k2,dt_ic2(i)) * wt(i)
		    rw(7*nndt+j) = -wt(i)
		 else
		    rw(4*nndt+j) = 0.0
		    rw(5*nndt+j) = 0.0
		    rw(6*nndt+j) = 0.0
		    rw(7*nndt+j) = 0.0
		 endif
	      endif
	      
c       Set up column indexes with non-zero elements
	      iw(1+nar+      j)  = 4*dt_ic1(i) - 3
	      iw(1+nar+  nndt+j) = 4*dt_ic1(i) - 2
	      iw(1+nar+2*nndt+j) = 4*dt_ic1(i) - 1
	      iw(1+nar+3*nndt+j) = 4*dt_ic1(i)
	      iw(1+nar+4*nndt+j) = 4*dt_ic2(i) - 3
	      iw(1+nar+5*nndt+j) = 4*dt_ic2(i) - 2
	      iw(1+nar+6*nndt+j) = 4*dt_ic2(i) - 1
	      iw(1+nar+7*nndt+j) = 4*dt_ic2(i)
	   endif
	enddo
c---    Check to make sure # of earthquake data is correct
	
	write(*,*)'j=',j
	
c       Scale G matrix so the L2 norm of each column is 1.
	write(log,'("~ scaling G columns ... ")')
	
c       G array scaling
	do i=1,4*nev
	   norm(i) = 0.0
	enddo

	do i=1,nar
	   if(abs(rw(i)).gt.2) then
	      write(log,*)nar, nev, iw(1+nar+i), i, rw(i), ndt
	      write(log,*)'Warning! rw(i)>2'	    
	   endif
	   norm(iw(1+nar+i)) = norm(iw(1+nar+i)) + rw(i)**2
	enddo
c	write(16,*)(norm(i),i=1,4*nev)
	
	do i=1,nev*4
	   norm(i) = sqrt(norm(i)/nndt)
c--- initialize old_index and new_index
	   old_index(i) = 0
	   new_index(i) = 0
	enddo
	
c---    Check which L2 norm is equal to zero
	j=0
	do i=1,4*nev
	   fix(i)=0 ! fix location or origin time
	   if(norm(i).gt.0) then
	      j=j+1
	      old_index(j)=i
	      new_index(i)=j
	      norm1(j)=norm(i)
	      fix(i)=1
	   endif
	enddo

	new_loc=j		! new effective # of unknowns

	write(*,*)'new_loc=:',new_loc,'  4*nev=',4*nev

	do i=1,new_loc
	   norm(i)=norm1(i)
	enddo
	
	do i=1, nar
	   i1=iw(1+nar+i)
	   if (fix(i1).eq.1) then ! not fixed
	      iw(1+nar+i)=new_index(i1)
	   else
	      iw(1+nar+i)=new_loc
	      rw(i)=0.0
c	      write(*,*)'Warning....Should not happen here!'
	   endif
	enddo

	do i=1,nar
	   if( iw(1+nar+i).ge.1.and.iw(1+nar+i).le.new_loc) then
	      rw(i) = rw(i) / norm(iw(1+nar+i))
	   else
	      write(16,*) 'Loc:Column index is not in the range'
	      write(16,*) iw(1+nar+i)
	   endif
	enddo

c     Testing...
	do i=1,new_loc
	   norm_test(i) = 0.0
	enddo
	do i=1,nar
	   if(iw(1+nar+i).ge.1.and.iw(1+nar+i).le.new_loc)
     &	   norm_test(iw(1+nar+i)) = norm_test(iw(1+nar+i)) + rw(i)**2
	enddo
	do i=1,new_loc
	   norm_test(i) = sqrt(norm_test(i)/nndt)
	   if (abs(norm_test(i)-1).gt.0.001) then
	      write(*,'("FATAL ERROR (lsqr: G scaling).")') 
	      stop
	   endif
	enddo
	
	
c       Least square fitting using the algorithm of
c       Paige and Saunders, acm-trans. math. software, vol.8, no. 2,
c       jun., 1982, p. 195.

c       Set up input parameter first
	m = nndt
	n = new_loc
	leniw = 2*nar+1
	lenrw = nar
	do i= 1,n
	   w1(i) = 0.0
	   w2(i) = 0.0
	   x(i) = 0.0
	   se(i) = 0.0
	enddo
	atol = 0.000001
	btol = 0.000001
	conlim = 100000.0
c	itnlim = 100*n
	itnlim = 4*n
	istop = 0
	anorm = 0.0
	acond = 0.0
	rnorm = 0.0
	arnorm = 0.0
	xnorm = 0.0
	
	call datetime(dattim)
	write(log,'("~ lsqr ...    ", a)') dattim
c	write(*,*)'nar,m,n'
c	write(*,*)nar,m,n
c       d= data vector; w1,w2 = workspace; x= solution vector; se=solution error
	call lsqr(m, n, damp, 
     & leniw, lenrw, iw, rw, 
     & d, w1, w2, x, se, 
     & atol, btol, conlim, itnlim, 
     & istop, anorm, acond, rnorm, arnorm, xnorm)

	write(log,'("  istop = ",i1,"; acond (CND)=",f8.1,"; anorm =",f8.1,
     & "; arnorm =",f8.1,"; xnorm =",f8.1)')
     & istop, acond, anorm, arnorm, xnorm

	if (nsrc.eq.1) nsrc = nev

c     Rescale model vector
	do i=1,new_loc
	   x(i) = x(i) / norm(i)
	   se(i) = se(i) / norm(i)
	enddo

c       Unweight and rescale G matrix
	do i=1,nndt
	   rw(i)        = rw(i)        * wtinv(i) * norm(iw(1+nar+i))
	   rw(nndt+i)   = rw(nndt+i)   * wtinv(i) * norm(iw(1+nar+nndt+i))
	   rw(2*nndt+i) = rw(2*nndt+i) * wtinv(i) * norm(iw(1+nar+2*nndt+i))
	   rw(3*nndt+i) = rw(3*nndt+i) * wtinv(i) * norm(iw(1+nar+3*nndt+i))
	   rw(4*nndt+i) = rw(4*nndt+i) * wtinv(i) * norm(iw(1+nar+4*nndt+i))
	   rw(5*nndt+i) = rw(5*nndt+i) * wtinv(i) * norm(iw(1+nar+5*nndt+i))
	   rw(6*nndt+i) = rw(6*nndt+i) * wtinv(i) * norm(iw(1+nar+6*nndt+i))
	   rw(7*nndt+i) = rw(7*nndt+i) * wtinv(i) * norm(iw(1+nar+7*nndt+i))
	enddo
	
c       Compute residuals from d = G*x
	do i=1,nndt
	   i1=old_data_index(i)
	   d(i) = -dt_res(i1)*1000.0
	enddo

	call aprod(1, m, n, x, d, leniw, lenrw, iw, rw)
c--- Map back residuals 
	do i=1,ndt
	   i1=new_data_index(i)
	   if(i1.gt.0) then !replace old reiduals, or keep it fixed 
	      dt_res(i) = -d(i1)/1000.0
	   endif   
	enddo

c     Get residual statistics (avrg, rms, var..)
	unknowns=new_loc
	call resstat_FDD(log, idata, ndt, unknowns, dt_res, wt, dt_idx, 
     & rms_cc, rms_ct, rms_cc0, rms_ct0, 
     & rms_ccold, rms_ctold, rms_cc0old, rms_ct0old, 
     &             resvar1)


c       Scale errors
c       The standard error estimates returned by LSQR increase monotonically
c       with the iterations.  If LSQR shuts down early because of loose tolerances,
c       or because the rhs-vector is special, the estimates will be too small.
c       (I think they are most likely to be accurate if the rhs is random.)
c
c       Remember that se(j) is covariance(j) / (m - n)
c       where m - n = 1000000.  I've never quite understood why we
c       divide by that number.

c       Errors for the 95% confidence level,
c       thus multiply the standard errors by 2.7955
	factor = 2.7955

c     initialize solution and error vectors
      do i=1, nev
	 src_dt(i)=0.0
	 src_dx(i)=0.0
	 src_dy(i)=0.0
	 src_dz(i)=0.0
	 src_et(i)=0.0
	 src_ex(i)=0.0
	 src_ey(i)=0.0
	 src_ez(i)=0.0
      enddo

c     Store solution and errors
      do j=1,new_loc
	 i=old_index(j)
	 i1=int(i/4)
	 i2=mod(i,4)
	 if(i2.eq.0) then !origin time
	    src_dt(i1) = -x(j)
	    src_et(i1) = sqrt(se(j)) * sqrt(resvar1) *factor
	 elseif(i2.eq.1) then ! x-coordinate
	    src_dx(i1+1) = -x(j)
	    src_ex(i1+1) = sqrt(se(j)) * sqrt(resvar1) *factor
	 elseif(i2.eq.2) then ! y-coordinate
	    src_dy(i1+1) = -x(j)
	    src_ey(i1+1) = sqrt(se(j)) * sqrt(resvar1) *factor
	 else
	    src_dz(i1+1) = -x(j) ! z-coordinate
	    src_ez(i1+1) = sqrt(se(j)) * sqrt(resvar1) *factor
	 endif
      enddo

      do i=1,nev
         src_cusp(i) = ev_cusp(i)
      enddo

c     Get average errors and vector changes
      exavold = exav
      eyavold = eyav
      ezavold = ezav
      etavold = etav
      dxavold = dxav
      dyavold = dyav
      dzavold = dzav
      dtavold = dtav
      exav = 0.0
      eyav = 0.0
      ezav = 0.0
      etav = 0.0
      dxav = 0.0
      dyav = 0.0
      dzav = 0.0
      do i=1,nev
         exav = exav + src_ex(i)
         eyav = eyav + src_ey(i)
         ezav = ezav + src_ez(i)
         etav = etav + src_et(i)
         dxav = dxav + abs(src_dx(i))
         dyav = dyav + abs(src_dy(i))
         dzav = dzav + abs(src_dz(i))
         dtav = dtav + abs(src_dt(i))
      enddo
      exav = exav/nev
      eyav = eyav/nev
      ezav = ezav/nev
      etav = etav/nev
      dxav = dxav/nev
      dyav = dyav/nev
      dzav = dzav/nev
      dtav = dtav/nev

      if (iter.eq.1) then
         exavold = exav
         eyavold = eyav
         ezavold = ezav
         etavold = etav
         dxavold = dxav
         dyavold = dyav
         dzavold = dzav
         dtavold = dtav
      endif

c     Output location statistics
      write(log,'(/,"Location summary:")')
      write(log,'(
     & " mean 2sig-error (x,y,z,t) [m,ms]: ",/,f7.1,f7.1,f7.1,f7.1,
     & " (",f7.1,f7.1,f7.1,f7.1")",/,
     & " mean shift (x,y,z,t) [m,ms] (DX,DY,DZ,DT): ",/,
     & f7.1,f7.1,f7.1,f7.1," (",f7.1,f7.1,f7.1,f7.1,")")')
     & exav, eyav, ezav, etav, exav-exavold, eyav-eyavold, 
     & ezav-ezavold, etav-etavold, 
     & dxav, dyav, dzav, dtav, dxav-dxavold, dyav-dyavold, 
     & dzav-dzavold, dtav-dtavold

      end !of subroutine lsfit_lsqr
