c---------------------------------------------------------------------
c Update by Haijiang Zhang on Jan. 19, 2005
c Change the program to make it more stable
c Part of codes is from Steve Roecker at RPI 
c-----------------------------------------------------------------         

	subroutine cover( nx, ny, nz, h, t, xs, ys, zs, x0, y0, 
     &                  z0,xr,yr,zr,ray,nray,travel_time,isgood)

      use tomoFDD
      integer      nx,ny,nz,i,j,k,is,js,ks,ish,jsh,ksh,iiiii,
     *             nseg,
     *             iscell,jscell,kscell,
     *             iseg,nk,nk2,nj,nj2,nray
      real         t(MAXGRID),du(MAXGRID/100),
     *             ray(3,MAXSEG)
      real x0,y0,z0,x,y,z,h,xi,yj,zk,xs,ys,zs,
     *             gradt(3),dd(3),d,length,
     *             dt,fx,fy,fz, gradtm
      real xr,yr,zr
      integer md, isgood
      

      real travel_time
c     character*80 sfile,tfile,rfile,ufile,pfile

      x=xr
      y=yr
      z=zr
c      write(16,*)'nx=,ny=,nz='
c      write(16,*)nx, ny, nz

      if  (nx*ny*nz.gt.MAXGRID)  then
         write (*,*) '***** ERROR:  dimensions are too big'
         stop
      end if

c      write(16,*)'xs,ys,zs',xs,ys,zs
c      write(16,*)'x0,y0,z0',x0,y0,z0
c      write(16,*)'h=',h
c      write(16,*)'x,y,z=',x,y,z
c
c     TEST FOR 2D MODELS
      if  (nx.eq.1)  then
         x0 = x0 - h/2.
         nx = 2
         write (*,*) '2d model encountered (nx=1)'
         if  (nx*ny*nz.gt.MAXGRID)  then
            write (*,*) '***** ERROR:  dimensions are too big (for 2d)'
            goto 800
         end if
         do 310, i=1,ny*nz
            du(i) = t(i)
 310     continue
         do 320, k=1,nz
            nk = nx*ny*(k-1)
            nk2 = 1*ny*(k-1)
            do 325, j=1,ny
               iiiii = nk+nx*(j-1)+1
               t(iiiii) = du(nk2+1*(j-1)+1)
               iiiii = iiiii+1
               t(iiiii) = du(nk2+1*(j-1)+1)
 325        continue
 320     continue
      endif
      if  (ny.eq.1)  then
         y0 = y0 - h/2.
         ny = 2
         write (*,*) '2d model encountered (ny=1)'
         if  (nx*ny*nz.gt.MAXGRID)  then
            write (*,*) '***** ERROR:  dimensions are too big (for 2d)'
            goto 800
         end if

         do 330, i=1,nx*nz
            du(i) = t(i)
 330     continue
         do 340, k=1,nz
            nk = nx*ny*(k-1)
            nk2 = nx*1*(k-1)
            do 345, i=1,nx
               iiiii = nk+i
               t(iiiii) = du(nk2+i)
               iiiii = iiiii+nx
               t(iiiii) = du(nk2+i)
 345        continue
 340     continue
      endif
      if  (nz.eq.1)  then
         z0 = z0 - h/2.
         nz = 2
         write (*,*) '2d model encountered (nz=1)'
         if  (nx*ny*nz.gt.MAXGRID)  then
            write (*,*) '***** ERROR:  dimensions are too big (for 2d)'
            goto 800
         end if
         nk = nx*ny
         do 350, i=1,nx
            do 360, j=1,ny
               iiiii = nk+nx*(j-1)+i
               t(iiiii) = t(nx*(j-1)+i)
 360        continue
 350     continue
      endif

c
c     FIND THE INTEGER SOURCE POINT
      is = nint((xs-x0)/h) + 1
      js = nint((ys-y0)/h) + 1
      ks = nint((zs-z0)/h) + 1
c      write(*,*)
c     FIND THE CELL CONTAINING THE SOURCE POINT
      iscell = int((xs-x0)/h) + 1
      jscell = int((ys-y0)/h) + 1
      kscell = int((zs-z0)/h) + 1
c
      nseg = 0
      iseg = 0
      isegm1 = 0
      isegm2 = 0

c
C     FIND CELL
      i = int((x-x0)/h) + 1
      if  (((x-x0).lt.0.) .or. (i.ge.nx))  then
        write (*,*) 'x out of bounds: x0,x,nx,i',x0,x,nx,i
	goto 800
      endif
      j = int((y-y0)/h) + 1
      if  (((y-y0).lt.0.) .or. (j.ge.ny))  then
        write (*,*) 'y out of bounds: y0,y,ny,j',y0,y,ny,j
	goto 800
      endif
      k = int((z-z0)/h) + 1
      if  (((z-z0).lt.0.) .or. (k.ge.nz))  then
        write (*,*) 'z out of bounds: z0,z,nz,k',z0,z,nz,k
	goto 800
      endif
      ish = i
      jsh = j
      ksh = k
      xi = h*(i-1) + x0
      yj = h*(j-1) + y0
      zk = h*(k-1) + z0
c
c     CALCULATE THE TRAVELTIME AT THE SHOTPOINT (TRILINEAR INTERPOLATION)
      nk = nx*ny*(k-1)
      nj = nx*(j-1)
      nk2 = nx*ny*k
      nj2 = nx*j
c     DON'T USE SHOT IF IT'S IN REGION OF MODEL WHERE TIMES WERE NOT CALCULATED
C        (SEE PARAMETER maxoff IN punch.c)
      if (t(nk+nj+i).gt.1.e9   .or. t(nk+nj+i+1).gt.1.e9  .or. 
     *    t(nk+nj2+i).gt.1.e9  .or. t(nk+nj2+i+1).gt.1.e9 .or. 
     *    t(nk2+nj+i).gt.1.e9  .or. t(nk2+nj+i+1).gt.1.e9 .or. 
     *    t(nk2+nj2+i).gt.1.e9 .or. t(nk2+nj2+i+1).gt.1.e9 ) goto 800
      fx = (x-xi)/h
      fy = (y-yj)/h
      fz = (z-zk)/h
      dt = (1.-fx)*(1.-fy)*(1.-fz)*t(nk+nj+i)
      dt = dt + fx*(1.-fy)*(1.-fz)*t(nk+nj+i+1)
      dt = dt + (1.-fx)*fy*(1.-fz)*t(nk+nj2+i)
      dt = dt + (1.-fx)*(1.-fy)*fz*t(nk2+nj+i)
      dt = dt + fx*fy*(1.-fz)*t(nk+nj2+i+1)
      dt = dt + fx*(1.-fy)*fz*t(nk2+nj+i+1)
      dt = dt + (1.-fx)*fy*fz*t(nk2+nj2+i)
      dt = dt + fx*fy*fz*t(nk2+nj2+i+1)
c      dist = sqrt((xs-x)**2+(ys-y)**2+(zs-z)**2)
c
      travel_time=dt

      length = 0.
c
C     *** BACK-PROPAGATE RAY ***
C     *** Propogate ray from earthquake to station
        nray=1
	ray(1,1)=x
	ray(2,1)=y
	ray(3,1)=z

 100  continue
c
c     FIND LOCATION OF CELL
      xi = h*(i-1) + x0
      yj = h*(j-1) + y0
      zk = h*(k-1) + z0
c
c     IF RAY IN SOURCE CELL, FINISH RAY BY ADDING TO COVERAGE
c     This is the original judgement from J. Hole.
c      if  ((i.eq.iscell .and. j.eq.jscell .and. k.eq.kscell) .or.
c     *     (((x-xs)**2+(y-ys)**2+(z-zs)**2).lt.1.e-6))  then
       if  ((i.eq.iscell .and. j.eq.jscell .and. k.eq.kscell) .or.
     *     (((x-xs)**2+(y-ys)**2+(z-zs)**2).lt.(h*h*1.e-6)))  then

c        *** LENGTH OF RAY IN SOURCE CUBE APPROXIMATED BY STRAIGHT LINE
         length = length + sqrt((x-xs)**2+(y-ys)**2+(z-zs)**2)
         
         if(ray(1,nray).eq.xs.and.ray(2,nray).eq.ys.and.
     &      ray(3,nray).eq.zs) then
            isgood=1
            goto 800
         else
	    ray(1,nray+1)=xs
	    ray(2,nray+1)=ys
	    ray(3,nray+1)=zs
            nray=nray+1
            isgood=1 
            goto 800
         endif
      endif
c
c     FIND THE RAY GRAD(T)
      nk = nx*ny*(k-1)
      nj = nx*(j-1)
      nk2 = nx*ny*k
      nj2 = nx*j
      gradt(1) = (t(nk+nj+i+1)+t(nk+nj2+i+1)
     *     +t(nk2+nj+i+1)+t(nk2+nj2+i+1)
     *     -t(nk+nj+i)-t(nk+nj2+i)
     *     -t(nk2+nj+i)-t(nk2+nj2+i)) / (4.*h)
      gradt(2) = (t(nk+nj2+i)+t(nk+nj2+i+1)
     *     +t(nk2+nj2+i)+t(nk2+nj2+i+1)
     *     -t(nk+nj+i)-t(nk+nj+i+1)
     *     -t(nk2+nj+i)-t(nk2+nj+i+1)) / (4.*h)
      gradt(3) = (t(nk2+nj+i)+t(nk2+nj+i+1)
     *     +t(nk2+nj2+i)+t(nk2+nj2+i+1)
     *     -t(nk+nj+i)-t(nk+nj+i+1)
     *     -t(nk+nj2+i)-t(nk+nj2+i+1)) / (4.*h)
c     IF RAY IN SOURCE (seismometer) CUBE, USE STRAIGHT RAY FROM SOURCE
      if  (i.ge.(is-2) .and. i.lt.(is+2) .and. 
     *     j.ge.(js-2) .and. j.lt.(js+2) .and. 
     *     k.ge.(ks-2) .and. k.lt.(ks+2))  then
         gradt(1) = x - xs
         gradt(2) = y - ys
         gradt(3) = z - zs
      endif
c      write (29,*) 'gradt',gradt
c
       gradtm = sqrt(gradt(1)*gradt(1)+ gradt(2)*gradt(2)+
     &           gradt(3)*gradt(3))
c     FIND THE POSSIBLE STEP SCALARS
      if  (gradt(1).gt.0.)  then
         dd(1) = (xi-x)/gradt(1)
         dlen = dd(1)*gradtm
c         write(*,*)'dlen=',dlen
         if (dlen.gt.-h/1000.)  then

c---The following tests prevent the ray from doubling back on itself.  The
c         first is simply to keep the ray from going back to where it's just been.
c         The second prevents a loop over 4 cells.   In this particular blockif,
c         iseg will be -1, so a previous iseg of 1 means that the ray is simply
c         trying to return to the cell it just left.  An example of a 4 cell loop
c         would be isegs of 20, 1, -20, where an additional -1 would bring it
c         back to the original cell, and hence set up the possibility of a closed
c         loop.

            if  ((iseg.eq.1) .or.
     *           ((iseg+isegm1+isegm2.eq.1).and.
     *            (nseg.ge.3)))  then
                  gradt(1) = 0.
                  dd(1) = -1e20    
            endif
          endif
      else if  (gradt(1).lt.0.)  then
         dd(1) = (xi+h-x)/gradt(1)
         dlen = dd(1)*gradtm
c         write(*,*)'dlen=',dlen
         if (dlen.gt.-h/1000.)  then
            if  ((iseg.eq.-1) .or.
     *           ((iseg+isegm1+isegm2.eq.-1).and.
     *            (nseg.ge.3)))  then
                  gradt(1) = 0.
                  dd(1) = -1e20    
            endif
         end if
      else
         dd(1) = -1e20
      end if
      if  (gradt(2).gt.0.)  then
         dd(2) = (yj-y)/gradt(2)
         dlen = dd(2)*gradtm
c         write(*,*)'dlen=',dlen
         if (dlen.gt.-h/1000.)  then
            if((iseg.eq.20) .or. 
     *          ((iseg+isegm1+isegm2.eq.20).and.
     *            (nseg.ge.3))) then
               gradt(2) = 0.
               dd(2) = -1e20
            endif
         end if
      else if  (gradt(2).lt.0.)  then
         dd(2) = (yj+h-y)/gradt(2)
         dlen = dd(2)*gradtm
         if (dlen.gt.-h/1000.)  then
             if( (iseg.eq.-20) .or.
     *          ((iseg+isegm1+isegm2.eq.-20).and.
     *            (nseg.ge.3))) then
               gradt(2) = 0.
               dd(2) = -1e20
            endif
         end if
      else
         dd(2) = -1e20
      end if
      if  (gradt(3).gt.0.)  then
         dd(3) = (zk-z)/gradt(3)
         dlen = dd(3)*gradtm
c         write(*,*)'dlen=',dlen
         if (dlen.gt.-h/1000.)  then
             if( (iseg.eq.300) .or.
     *          ((iseg+isegm1+isegm2.eq.300).and.
     *            (nseg.ge.3))) then
                  gradt(3) = 0.
                  dd(3) = -1e20    
            endif            
         end if
      else if  (gradt(3).lt.0.)  then
         dd(3) = (zk+h-z)/gradt(3)
         dlen = dd(3)*gradtm
c         write(*,*)'dlen=',dlen
         if (dlen.gt.-h/1000.)  then
             if( (iseg.eq.-300) .or.
     *          ((iseg+isegm1+isegm2.eq.-300).and.
     *            (nseg.ge.3))) then
               gradt(3) = 0.
               dd(3) = -1e20
            endif
         end if
      else
         dd(3) = -1e20
      end if

c     DETERMINE WHICH IS DESIRED (SMALLEST ABSOLUTE VALUE; ALL SHOULD 
c        BE NEGATIVE, BUT NEAR-ZERO MIGHT BE POSITIVE)
      if  (abs(dd(1)).le.abs(dd(2)).and.abs(dd(1)).le.abs(dd(3)))  then
         md = 1
         d = dd(1)
      else if (abs(dd(2)).le.abs(dd(3)))  then
         md = 2
         d = dd(2)
      else
         md = 3
         d = dd(3)
      end if

c     CONTINUE FINDING RAY
      length = length - d*sqrt(gradt(1)**2+gradt(2)**2+gradt(3)**2)
      x = x + d*gradt(1)
      y = y + d*gradt(2)
      z = z + d*gradt(3)
      if (nray.lt.1) write(*,*)'Weired thing happened in cover'
      if(ray(1,nray).ne.x.or.ray(2,nray).ne.y.or.ray(3,nray).ne.z) then
         nray=nray+1
         if(nray.gt.MAXSEG) then
            write(*,*)'nray is too large:',nray   
            write(*,*)'nx, ny, nz=',nx,ny,nz
            write(*,*)'xs, ys, zs=',xs,ys,zs
            write(*,*)'xr, yr, zr=',xr,yr,zr
            stop
         endif
         ray(1,nray)=x
         ray(2,nray)=y
         ray(3,nray)=z
      endif

      nseg = nseg + 1
      isegm2 = isegm1
      isegm1 = iseg

      if  (nseg.gt.MAXSEG)  then
c         write (*,*)'ray too long'
         isgood=0
         goto 800
      endif
      if  (md.eq.1)  then
         if  (gradt(1).ge.0) then
            i=i-1
            if  (i.lt.1)  goto 800
            iseg = -1
         else
            i=i+1
            if  (i.ge.nx) goto 800
            iseg = 1
         end if
      else if  (md.eq.2)  then
         if  (gradt(2).ge.0) then
            j=j-1
            if  (j.lt.1)  goto 800
            iseg = -20
         else
            j=j+1
            if  (j.ge.ny)  goto 800
            iseg = 20
         end if
      else
         if  (gradt(3).ge.0) then
            k=k-1
            if  (k.lt.1)  goto 800
            iseg = -300
         else
            k=k+1
            if  (k.ge.nz)  goto 800
            iseg = 300
         end if
      end if

      goto 100
c300  continue
 
800   continue
      end
