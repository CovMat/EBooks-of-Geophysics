
	subroutine trialsrc_FDD_shot(istart, sdc0_lat, sdc0_lon, sdc0_dep,
     &	nev, ev_cusp, ev_lat, ev_lon, ev_dep, ev_type,
     &	nsrc, src_cusp, src_lat0, src_lon0,
     &	src_x0, src_y0, src_z0, src_t0,
     &	src_lat, src_lon, src_dep, src_type,
     &	src_x, src_y, src_z, src_t, wlat, wlon, theta)


	use tomoFDD
	implicit none


c	Parameters:
	integer		istart
	real		sdc0_lat	! Cluster center
	real		sdc0_lon	! Cluster center
	real		sdc0_dep	! Cluster center
	integer		nev		! No. of events
	integer		ev_cusp(MAXEVE)	! [1..nev]
	real		ev_lat(MAXEVE)	! [1..nev]
	real		ev_lon(MAXEVE)	! [1..nev]
	real		ev_dep(MAXEVE)	! [1..nev]
	integer         ev_type(MAXEVE)
	integer		nsrc		! No of trial sources
	integer		src_cusp(MAXEVE)! [1..nev]
	real		src_lat0(MAXEVE)! [1..nev]
	real		src_lon0(MAXEVE)! [1..nev]
	real		src_x0(MAXEVE)	! [1..nev]
	real		src_y0(MAXEVE)	! [1..nev]
	real		src_z0(MAXEVE)	! [1..nev]
	real		src_t0(MAXEVE)	! [1..nev]
	doubleprecision	src_lat(MAXEVE)	! [1..nev]
	doubleprecision	src_lon(MAXEVE)	! [1..nev]
	real		src_dep(MAXEVE)	! [1..nev]
	real		src_x(MAXEVE)	! [1..nev]
	real		src_y(MAXEVE)	! [1..nev]
	real		src_z(MAXEVE)	! [1..nev]
	real		src_t(MAXEVE)	! [1..nev]
	integer         src_type(MAXEVE)
	real            wlat
	real            wlon
	real            theta
c	Local variables:
	integer	i
	real	x, y, z
	real    lat, lon, dep

c     Set up parameters for initial inversion
      if (istart.eq.1) then
c        Cluster center as initial trial source
         nsrc = 1
         do i=1,nev
            src_cusp(i) = 0
            src_lon(i) = sdc0_lon
            src_lat(i) = sdc0_lat
            src_dep(i) = sdc0_dep
            src_x(i) = 0.0
            src_y(i) = 0.0
            src_z(i) = 0.0
            src_t(i) = 0.0
c           Store initial trial source
            src_lon0(i) = sdc0_lon
            src_lat0(i) = sdc0_lat
            src_x0(i) = 0.0
            src_y0(i) = 0.0
            src_z0(i) = 0.0
            src_t0(i) = 0.0
	    src_type(i) = ev_type(i)
         enddo
      else
c        Catalog sources as initial trial source
c        Add noise for synthetic data mode
         nsrc = nev
         do i=1,nev
            src_cusp(i) = ev_cusp(i)
            src_lon(i) = ev_lon(i)
            src_lat(i) = ev_lat(i)
            src_dep(i) = ev_dep(i)
            src_lon0(i) = ev_lon(i)
            src_lat0(i) = ev_lat(i)
	    lat=src_lat(i)
	    lon=src_lon(i)
	    dep=src_dep(i)
	    call sph2car_ft(lat,lon,dep,wlat,wlon,x,y,z,theta)
            src_x(i) = x *1000.0
            src_y(i) = y *1000.0
	    src_z(i) = z *1000.0
            src_t(i) = 0.0
            src_x0(i) = src_x(i)
            src_y0(i) = src_y(i)
            src_z0(i) = src_z(i)
            src_t0(i) = src_t(i)
	    src_type(i)= ev_type(i)
         enddo
      endif
      end !of subroutine trialsrc_FDD
