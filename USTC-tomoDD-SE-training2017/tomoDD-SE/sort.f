c Sort a float vector

	subroutine sort(n, ra)

	implicit none

c	Parameters:
	integer	n
	real	ra(n)	! (1..n)

c	Local variables:
	integer	i, j, l
 	integer	ir
	real	rra

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/sort.f,v 1.7 2001/03/08 19:24:15 dinger Exp julian $"/
	save rcsid

      if (n.le.1) then
         return
      endif


      l = n/2+1
      ir = n
10    continue
         if (l.gt.1) then
            l = l-1
            rra = ra(l)
         else
            rra = ra(ir)
            ra(ir) = ra(1)
            ir = ir-1
            if (ir.eq.1) then
               ra(1) = rra
               return
            endif
         endif
         i = l
         j = l+l
20       if (j.le.ir) then
            if (j.lt.ir) then
               if (ra(j).lt.ra(j+1)) j = j+1
            endif
            if (rra.lt.ra(j)) then
               ra(i) = ra(j)
               i = j
               j = j+j
            else
               j = ir+1
            endif
         go to 20
         endif
         ra(i) = rra
      go to 10
      end
