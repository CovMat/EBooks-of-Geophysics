c Find specified value in ordered integer vector

	integer function ifindi(n, ia, iv)

	implicit none

c	Parameters:
	integer		n
	integer		ia(n)	! [1..n] Vector to search
	integer		iv	! Value to find

c	Local variables:
	integer		i
	integer		k	! 2^(no. of chops)

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/ifindi.f,v 1.6 2001/03/08 19:16:07 dinger Exp julian $"/
	save rcsid

      if (n.le.0) then
         ifindi = 0
         return
      endif

      if (iv.lt.ia(1) .or. iv.gt.ia(n)) then
c        Outside range of vector
         ifindi=0
         return
      endif

      k = 2
      i = nint(real(n)/k)
10    if (k.gt.2*n) then
c        Value not in vector
         ifindi = 0
         return
      endif
      k = k*2
      if (iv.lt.ia(i)) then
c        Value smaller:  Search below
         i = i-nint(real(n)/k)
         goto 10
      endif
      if (iv.gt.ia(i)) then
c        Value larger:  Search above
         i = i+nint(real(n)/k)
         goto 10
      endif

c     Value found: iv == ia[i]
      ifindi = i
      return
      end
