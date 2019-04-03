c Convert latitude and longitude to kilometers relative
c to center of coordinates by short distance conversion.

	subroutine dist(xlat, xlon, xkm, ykm)

	implicit none

c	Parameters:
	doubleprecision	xlat, xlon	! (input)
	real		xkm, ykm	! (output)

c	Local variables:
	doubleprecision lat1, lat2, lat3
	real	q
	real	xx
	real	yp

	include "geocoord.inc"

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/dist.f,v 1.6 2001/02/06 03:04:35 julian Exp julian $"/
	save rcsid

c Set up short distance conversion by subr. SETORG
      q=60*xlat-olat
      yp=q+olat
      lat1=datan(rlatc*dtan(RAD*yp/60.0))
      lat2=datan(rlatc*dtan(RAD*OLAT/60.0))
      LAT3=(LAT2+LAT1)/2.
      xx=60*xlon-olon  !  - wegen LON E
      q=q*aa
      xx = xx*bb*dcos(LAT3)
      IF(rotate.ne.0.) then
c** rotate coordinate system anticlockwise
        yp=cost*q+sint*xx
        xx=cost*xx-sint*q
        q=yp
      ENDIF

      xkm=xx
      ykm=q

      return
      end
