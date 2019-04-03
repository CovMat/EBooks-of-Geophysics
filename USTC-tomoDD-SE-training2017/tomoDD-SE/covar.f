	subroutine covar(maxev0, v, n, q, cvm)

	implicit none

c	Parameters:
	integer	maxev0
	real	v(maxev0*4, maxev0*4)
	integer	n
	real	q(maxev0*4)
	real	cvm(maxev0*4, maxev0*4)

c	Local variables:
	integer	i, j, k		! Dummy loop indices
	real	sum

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/covar.f,v 1.2 2001/02/06 23:16:55 julian Exp julian $"/
	save rcsid

      do i=1,n
         do j=1,i
           sum= 0
           do k=1,n
              if(q(k).ne.0) sum= sum + v(i,k)*v(j,k) * (1/(q(k)*q(k)))
           enddo
           cvm(i,j)= sum
           cvm(j,i)= sum
         enddo
      enddo
      return
      end
