cccccccccccccccccccc Use the dense storage formatccccccccccccccc
ccccccccccccccccccccccc subroutine ccccccccccccccccccccccccccccccccccccc

      subroutine densemv(transa,m,n,x,y,dparm,iparm)

      use matvec
      implicit none
      character*1 transa
      integer m,n,iparm(*)
      real x(*),y(*),dparm(*)

      call sgemv(transa,m,n,1.0,der2,densemmax,x,1,0.0,y,1)

      end

