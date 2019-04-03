c Find a free fortran i/o unit-number

	subroutine FREEUNIT(iunit)

	implicit none

c	Parameters:
	integer	iunit

c	Local variables:
	logical	lopen

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/freeunit.f,v 1.4 2001/02/03 01:19:53 julian Exp julian $"/
	save rcsid

      do iunit=10,999
         if(iunit.eq.999) stop'FREEUNIT>>> no free unit found!'
         inquire(unit=iunit,opened=lopen)
         if(.not.lopen) RETURN
      enddo
      RETURN
      end ! of subr. freeunit
