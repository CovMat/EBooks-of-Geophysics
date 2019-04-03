c Matrix mutliply: c(m,m) = a(m,m)*b(m,m)

	subroutine matmult1(maxm, m, a, b, c)

	implicit none

c	Parameters:
	integer	maxm
	integer	m
	real	a(maxm,maxm)	! (1..maxm, 1..maxm)
	real	b(maxm,maxm)	! (1..maxm, 1..maxm)
	real	c(maxm,maxm)	! (1..maxm, 1..maxm)

c	Local variables:
	integer	j, l, k		! Dummy loop indices
	real	sum

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/matmult1.f,v 1.4 2001/02/13 23:55:34 julian Exp julian $"/
	save rcsid

      do j=1,m
         do k=1,m
            sum=0.0
            do l=1,m
               sum=sum+a(j,l)*b(l,k)
            end do
            c(j,k)=sum
         end do
      end do
      return
      end
