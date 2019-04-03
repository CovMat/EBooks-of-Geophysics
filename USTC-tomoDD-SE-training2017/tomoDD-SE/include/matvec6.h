	module matvec

c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c
      integer*8 s_rw,s_iw
      integer mmax,nmax,kmax,nnzmax,mmax2
      integer densemmax, densenmax
      integer lwrk,liwrk
      parameter(densemmax=10000,densenmax=10000)
      parameter(mmax=9000,nmax=9000,kmax=9000,mmax2=700000)
      parameter(nnzmax=densemmax*densenmax*10)
      parameter(lwrk=mmax+nmax+13*kmax+8*kmax**2+32*mmax+8)
      parameter(liwrk = 8*kmax)

      real,allocatable :: rw(:)
      integer,allocatable :: iw(:)
      real,allocatable :: U(:,:)
      real,allocatable :: V(:,:)
      real,allocatable :: Sigma(:)
      real,allocatable :: bnd(:)
      real,allocatable :: work(:)
      integer,allocatable :: iwork(:)
      real,allocatable :: der2(:,:)
      real,allocatable :: U2(:,:)

      contains

c       allocate the arrays
        subroutine alloc_matvec()

        allocate(rw(nnzmax))
        if (.not.allocated(rw)) print *, "Memory allocation ERROR for rw"
        s_rw=4
        s_rw=(s_rw*nnzmax)
        print *, "Size of rw=",s_rw

        allocate(iw(2*nnzmax))
        if (.not.allocated(iw)) print *, "Memory allocation ERROR for iw"
        s_iw=4
        s_iw=(s_iw*2*nnzmax)
        print *, "Size of iw=",s_iw

        allocate(U(mmax,kmax+1))
        if (.not.allocated(U)) print *, "Memory allocation ERROR for U"

        allocate(V(nmax,kmax+1))
        if (.not.allocated(V)) print *, "Memory allocation ERROR for V"

        allocate(Sigma(kmax))
        if (.not.allocated(Sigma)) print *, "Memory allocation ERROR for Sigma"

        allocate(bnd(kmax))
        if (.not.allocated(bnd)) print *, "Memory allocation ERROR for bnd"

        allocate(work(lwrk))
        if (.not.allocated(work)) print *, "Memory allocation ERROR for work"

        allocate(iwork(liwrk))
        if (.not.allocated(iwork)) print *, "Memory allocation ERROR for iwork"
 
        allocate(der2(densemmax,densenmax))
        if (.not.allocated(der2)) print *, "Memory allocation ERROR for iwork"

	allocate(U2(densemmax,densenmax))
	if (.not.allocated(U2)) print *, "Memory allocation ERROR for iwork"

        end subroutine

        end module

