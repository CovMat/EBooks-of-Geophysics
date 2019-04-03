cccccccccccccccccccc Use the storage format used in LSQRccccccccccccccc
ccccccccccccccccccccccc subroutine ccccccccccccccccccccccccccccccccccccc
      subroutine lsqrmv5(transa,m,n,x,y,dparm,iparm)

c--- This also considers the damping factor for the normal equation
      use matvec

      character*1 transa
      integer m,n,iparm(*)
      real x(*),y(*),dparm(*)
c     Local variables:
      real t(mmax2)
      integer i1, i
      integer j1
      integer k
      integer kk
      logical lsame
      external lsame

      kk=iw(1)
      i1=1
      j1=kk+1

c     main iteration loop

c--- first calculate G'G
c--- y<-A*x

         do i=1,mmax2
            t(i) = 0.0
         enddo

         do k = 1,kk
            t(iw(i1+k)) = t(iw(i1+k)) + rw(k)*x(iw(j1+k))
         enddo

         do i=1,n
            y(i) = dparm(1)*dparm(1)*x(i)
         enddo
c     compute  y = a(transpose)*x
         do k = 1,kk
            y(iw(j1+k)) = y(iw(j1+k)) + rw(k)*t(iw(i1+k))
         enddo

      end
