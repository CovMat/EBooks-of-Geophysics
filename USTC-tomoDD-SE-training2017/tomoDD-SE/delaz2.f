c Compute distance and azimuth on a sphere

	subroutine delaz2(alat, alon, blat, blon, del, dist, az)

	implicit none

	doubleprecision	alat, alon	! Coordinates of first point
	real		blat, blon	! Coordinates of second point
	real		del		! Sentral angle (degrees)
	real		dist		! Distance (km)
	real		az		! Azimuth from a to b (degrees)

C	Local variables:
	doubleprecision	acol, bcol
	doubleprecision	azr
	doubleprecision	blatr, blonr
	doubleprecision	colat
	doubleprecision	cosdel
	doubleprecision	delr
	doubleprecision	flat
	doubleprecision	geoa
	doubleprecision	geob
	doubleprecision	rad
	doubleprecision	radius
	doubleprecision alatr, alonr
	doubleprecision diflon
	doubleprecision pi2
	doubleprecision tana, tanb
	real		xtop, xden

c	Built-in functions: Declarations not needed
	doubleprecision dtan
	doubleprecision	datan
	doubleprecision	dsin
	doubleprecision	dcos
	doubleprecision	dacos

c	doubleprecision top,den

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/delaz2.f,v 1.6 2001/02/09 01:48:17 julian Exp julian $"/
	save rcsid

	data pi2/1.570796d0/
	data rad/1.745329d-02/
	data flat/.993231d0/

c-----convert to radians
	alatr=alat*rad
	alonr=alon*rad
	blatr=blat*rad
	blonr=blon*rad
c-----convert latitudes to geocentric colatitudes
	tana=flat*dtan(alatr)
	geoa=datan(tana)
	acol=pi2-geoa
	tanb=flat*dtan(blatr)
	geob=datan(tanb)
	bcol=pi2-geob
c-----calculate delta
	diflon=blonr-alonr
	cosdel=dsin(acol)*dsin(bcol)*dcos(diflon)+dcos(acol)*
     &	dcos(bcol)
	delr=dacos(cosdel)
c-----calculate azimuth from a to b

c*****	Note the use of single precision xtop and xden instead
c	of the double precision top and den in the original
c	program.
c*****	Note also the call to atan2 instead of datan2.
c	Both of these changes were made so that dyn.load
c	would work in Splus.  For some reason, the ld command
c	ld -r -d didn't find _d_atan2
c						WLE 10/16/91

	xtop = dsin(diflon)
	xden=(dsin(acol)/dtan(bcol))-dcos(diflon)*dcos(acol)
	azr=atan2(xtop,xden)
c----- convert to degrees
	del=delr/rad
	az=azr/rad
	if(az.lt.0.0) az=360.+az
c-----compute distance in kilometers
	colat=pi2-(alatr+blatr)/2.d0
	radius=6378.163d0*
     &	(1.d0+3.35278d-3*((1.d0/3.d0)-(dcos(colat)**2)))
	dist=delr*radius
	return
c  ***** end of subroutine delaz *****
	end
