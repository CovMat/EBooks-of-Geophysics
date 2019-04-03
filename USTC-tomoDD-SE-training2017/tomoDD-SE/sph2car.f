	program sph2car
	real	wlat
	real	wlon
	real	xc
	real	yc
	real	zc
	real    tlat
	real	tlon
	real	tdep, dep
	real	theta

	open(10, file="borehole.sph")
	open(15, file="borehole.car")
        wlat=37.5
        wlon=-122
        theta=-36
	
	do i=1,106
	   read(10,*)tlat, tlon, dep
           tdep=0.0
           do j=1, 1000000
              call sph2car_ft(tlat,tlon,tdep,wlat,wlon,xc,yc,zc,theta)
              write(15,*)i,xc,yc,zc,tlat,tlon
              tdep=tdep+0.1
              if(tdep.gt.dep/1000) goto 999
           enddo
999        continue
        enddo

	end 
