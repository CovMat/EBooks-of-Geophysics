c Matrix mutliply: c(m,n) = a(m) * b(m,n)
c overwrites matrix b

	subroutine matmult3(maxm, maxn, m, n, a, b)

	implicit none

c	Parameters:
	integer	maxm
	integer	maxn
	integer	m
	integer	n
	real	a(maxm)
	real	b(maxm,maxn)

c	Local variables:
	integer i, j		! Dummy loop indices

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/matmult3.f,v 1.4 2001/02/09 01:57:14 julian Exp julian $"/
	save rcsid

      do i=1,m
            do j=1,n
                  b(i,j) = a(i)*b(i,j)
            enddo
      enddo
      return
      end
