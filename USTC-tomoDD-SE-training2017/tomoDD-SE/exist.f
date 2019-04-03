	subroutine exist(fn)

	implicit none

c	Parameters:
	character	fn*80

c	Local variables
	logical	ex

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/exist.f,v 1.3 2001/02/03 01:13:20 julian Exp julian $"/
	save rcsid

      inquire(FILE=fn,exist=ex)
      if (.not.ex) then
          write(*,'("FILE DOES NOT EXIST / CHECK IDATA,IPHASE: ",a)')fn
          stop
      endif
      end
