c Length of character string, excluding trailing blanks

	integer function trimlen(t)

	implicit none

c	Parameter:
	character t*(*)

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/trimlen.f,v 1.1 2001/02/03 02:55:40 julian Exp julian $"/
	save rcsid

      do 1 trimlen=LEN(t),1,-1
    1    if(t(trimlen:trimlen).ne.' ')RETURN
      trimlen=1
      end ! of integer function trimlen
