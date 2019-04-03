
      subroutine rayl2(isp,w,ni,delay)
      implicit real*8 (a-h, o-z)
      include 'RaySPDR.inc'
      real*8 la,lo
      dimension  w(3,msg+1)
      integer    in(8)
     

c     dinter - is the step for interpolation
c     finding the boundary of cell

      dinter=1.000
      rlength=0.0e10
      rdelay =0.0e10
      r2d = 90./asin(1.)
      dpi = asin(1.)/ 90.
* ibl is the flag that is writen in order to mark
* rays penetraiting lower 2900 km.
      ibl=0

      nrp1=ni
      pl=0.0
c     open(32,file='a.dat')
      call dist2(w(2,1),w(3,1),w(1,1),
     &                w(2,ni+1),w(3,ni+1),w(1,ni+1),dst)
      
c      write(*,*)'dst,theoretic time:',dst, dst/5.0

      dst=0.0
      pl =0.0
      do i=1,nrp1
         i1=i+1

         call dist2(w(2,i),w(3,i),w(1,i),
     &                w(2,i+1),w(3,i+1),w(1,i+1),dst)

         pl=pl+dst
c  decide on number of subsegments and compute length
         nms=nint(dst/dinter)+1

c        find dx, dy and dz for each step

         ar=w(2,i)*dpi
         as=w(2,i+1)*dpi
         br=w(3,i)*dpi
         bs=w(3,i+1)*dpi
         rr=w(1,i)
         rs=w(1,i+1)

         x1 = rr*sin(ar)*cos(br)
         y1 = rr*sin(ar)*sin(br)
         z1 = rr*cos(ar)
         x2 = rs*sin(as)*cos(bs)
         y2 = rs*sin(as)*sin(bs)
         z2 = rs*cos(as)

         dx = (x2-x1)/nms
         dy = (y2-y1)/nms
         dz = (z2-z1)/nms

         x  = x1 - 0.5*dx
         y  = y1 - 0.5*dy
         z  = z1 - 0.5*dz

         ssl=dst/nms

         do j=1,nms

c---------gradually increase distance toward the second point

          x = x + dx
          y = y + dy
          z = z + dz
          r = sqrt(x**2 + y**2 + z**2)

          acosa=z/r
          if(acosa.LT.-1.)acosa=-1.
          if(acosa.GT.1)acosa=1.
          la = acos(acosa)*r2d
          acosa=(x/r)/sin(la*dpi)
          if(acosa.LT.-1.)acosa=-1.
          if(acosa.GT.1.)acosa=1.
          lo = acos(acosa)*r2d
          if(y.LT.0.00000)lo=360.00000-lo

          call cellvel(isp,la,lo,r,Velc,nbl)
c         write(*,*)ip,jp,kp,Velc
c         write(*,*)ssl,Velc
          rdelay =rdelay+ssl/Velc
c         write(32,*)i, dst, ssl,Velc 

        enddo ! end of this segment
      enddo ! end of all ray segments

c     write(*,*)"travel time 3D and segment length ="
c     write(*,*)rdelay, pl

      delay = rdelay

c     close(32)
      return

      end

      subroutine dist2(rla1,rlo1,r1,rla2,rlo2,r2,dst)
      implicit real*8 (a-h, o-z)
      dpi = asin(1.)/ 90.

      as=rla1*dpi
      ar=rla2*dpi
      bs=rlo1*dpi
      br=rlo2*dpi
      rs=r1
      rr=r2
      x1 = rr*sin(ar)*cos(br)
      y1 = rr*sin(ar)*sin(br)
      z1 = rr*cos(ar)
      x2 = rs*sin(as)*cos(bs)
      y2 = rs*sin(as)*sin(bs)
      z2 = rs*cos(as)
      dx = (x2-x1) 
      dy = (y2-y1) 
      dz = (z2-z1) 
      dst=sqrt(dx**2+dy**2+dz**2)
      return
      end



      subroutine km2deg(ala,alo,adp,dx,dy,bla,blo,bdp)
      implicit real*8 (a-h, o-z)
c     This subroutine calculate position of new point
c     in polar coordinates basing on the coordinates
c     of main point in radians ( la is colatitude) and dx and dy in kilometers
      dpi = asin(1.)/ 90.
      dps=adp*SIN(ala)
      blo=alo+atan2(dx,dps)
      bla=ala+atan2(dy,adp)
        if(bla.gt.(180.*dpi))then
           bla=360.*dpi-bla
           blo=blo+180.*dpi
        endif
        if(bla.lt.0.)then
           bla=abs(bla)
           blo=blo+180.*dpi
        endif
        if(blo.lt.0.)blo=360.*dpi+blo
        if(blo.gt.(360.*dpi))blo=blo-(360.*dpi)
      bdp=sqrt(adp**2+dx**2+dy**2)
      return
      end

c                                                                       
c 3d pseudo-bending for the continuous, spherical earth
c   based on kazuki koketsu (eri, univ. tokyo) algorithm
c
c  ray is tracing from receiver to source
c  input coordinates are in degrees and km
c  latitude from -90 to 90 and longitude from -180 to 180
c  depth is radius from the earth's center

c 10/24/2007, output the calculated travel times

      subroutine pbr(isp,aas, bbs, hs, aar, bbr, hr, w, ni, tt)    

      implicit real*8 (a-h, o-z)                                        
      include 'RaySPDR.inc'

      dimension  r(msg+1), a(msg+1), b(msg+1)
      dimension  w(3,msg+1)
      integer    ni,i
      external   vel2,rtim,rlen
      real*8     RNULL
      common /coord/ shiftlo
c     data       ro / 6378.00 /
      data       ro / 6371.00 /
      data       RNULL /0.0e10/
      integer    true
      real       tt
      
c     write(*,*)aas, bbs, hs, aar, bbr, hr
c                                                                       
c       parameters for calculation                                      
c         xfac   = enhancement factor (see um and thurber, 1987)        
c         nloop  = number of bending iterations
c         n1, n2 = min & max of ray segments
c         mins  = min. length of segment (km)
c       initialization                                                  

c      flim   = 1.e-2/100. = 0.0001 sec
c      flim   = 1.e-2/10. = 0.001 sec


      ni     = msg+1
      xfac   = 1.5
      n1     = 2
      n2     = msg
      nloop  = 200
!     flim   = 1.e-5/10
!     mins   = 5.

!     flim   = 1.e-6/10
!     mins   = 2.

!      flim   = 1.e-5
!      mins   = 10.

!      flim   = 1.e-6
!      mins   = 5.

      flim   = 1.e-5
      mins   = 10. 

      dpi = asin(1.)/ 90.
      r2d = 90./asin(1.)

c     Check coordinates
     
      if(aas.LT.-90.OR.aas.GT.90.)then
        write(*,*)'Latitude of source is out of range in RaySPDR2.f'
        stop
      endif
      if(aar.LT.-90.OR.aar.GT.90.)then
        write(*,*)'Latitude of station is out of range in RaySPDR2.f'
        stop
      endif
      if(bbs.LT.-180.OR.bbs.GT.180.)then
        write(*,*)'Longitude of source is out of range in RaySPDR2.f'
        stop
      endif
      if(bbr.LT.-180.OR.bbr.GT.180.)then
        write(*,*)'Longitude of station is out of range in RaySPDR2.f'
        stop
      endif


c     Rotate coordinates in order to have
c     longitude and latitude range from 0 to 180. 
c     This program does not work with angles
c     greater than 180.       

c     Pass from latitude to colatitude

c      as = (90.00-aas) * dpi !colatitude
c      ar = (90.00-aar) * dpi !colatitude

c--- consider the ellipticity correction
      !aas2=atan(.99664719*tan(aas*dpi))
      aas2=atan(1.00000000*tan(aas*dpi))
      as  =asin(1.) - aas2

      !aar2=atan(.99664719*tan(aar*dpi))
      aar2=atan(1.00000000*tan(aar*dpi))
      ar  =asin(1.) - aar2


      if(bbr.LT.0.0)then
         bre=360.+bbr
      else
         bre=bbr
      endif

      if(bbs.LT.0.0)then
         bso=360.+bbs
      else
         bso=bbs
      endif
      dlo=abs(bso-bre)

      if(dlo.LT.180.)then


c     write(*,*)dlo

      shiftlo=0.0e10
      if(bso.LT.bre)then
         shiftlo=bso-(180.-dlo)/2.
         bbs=(180.-dlo)/2.
         bbr=bbs+dlo
      else
         shiftlo=bre-(180.-dlo)/2.
         bbr=(180.-dlo)/2.
         bbs=bbr+dlo
      endif

      else

      dlo=360.0000-dlo

      shiftlo=0.0e10
      if(bso.LT.bre)then
         shiftlo=bso-(dlo+(180.-dlo)/2.)
         bbs=(180.-dlo)/2.+dlo 
         bbr=bbs-dlo
      else    
         shiftlo=bre-(dlo+(180.-dlo)/2.)
         bbr=(180.-dlo)/2.+dlo
         bbs=bbr-dlo
      endif
 

      endif
c     write(*,*)'bbr,bbs,shift'
c     write(*,*)'aar,aas'
c     write(*,*)aar,aas
c     write(*,*)'bbr,bbs,shiftlo'
c     write(*,*)bbr,bbs,shiftlo

      bs = bbs * dpi
      br = bbr * dpi  

      ad = (as + ar) / 2.                                               

c--- change sign here
c      rs = ro - hs
c      rr = ro + hr

      rs = ro - hs
      rr = ro - hr
                                                      
c
c *** initial straight ray ***                                           
c       ni : number of ray segments
      ni = n1
      x1 = rr*sin(ar)*cos(br)                                       
      y1 = rr*sin(ar)*sin(br)                                       
      z1 = rr*cos(ar)                                               
      x2 = rs*sin(as)*cos(bs)                                       
      y2 = rs*sin(as)*sin(bs)                                       
      z2 = rs*cos(as)       
      dx = (x2-x1) / ni
      dy = (y2-y1) / ni
      dz = (z2-z1) / ni
c     write(*,*)(x2-x1),(y2-y1),(z2-z1)
c     write(*,*)sqrt((x2-x1)**2+(y2-y1)**2)
      do 155  j=1,ni+1
          x = x1 + dx*(j-1)
          y = y1 + dy*(j-1)
          z = z1 + dz*(j-1)                                             
          r(j) = sqrt(x**2 + y**2 + z**2)                         
          acosa=z/r(j)
          if(acosa.LT.-1.)acosa=-1.
          if(acosa.GT.1)acosa=1.
          a(j) = acos(acosa)                                   
          acosa=x/r(j)/sin(a(j))
          if(acosa.LT.-1.)acosa=-1.
          if(acosa.GT.1)acosa=1.
          b(j) = acos(acosa)
          if(y.LT.0.00000)b(j)=360.00000*dpi-b(j)

cc    write(*,*)a(j)*r2d,b(j)*r2d,r(j)
 155  continue
c
      to = rtim(isp,ni+1,r,a,b)
      tp = to
c      write(6,*) to
      do 112  i=1,ni+1
          w(1,i) = r(i)
          w(2,i) = a(i)
          w(3,i) = b(i)
 112  continue
c *** number of points loop ***
      loops = 0

c---- This is the old loop
c     do while(ni .le. n2)

cz--- change the convergence criterion    
      true=1
      do while(true.eq.1)  
          loops = 0
c                                                                       
c *** interation loop ***                                               
          do 3000  l=1,nloop                                                 
              loops = loops + 1
              do 2250  kk=2,ni
c               see um & thurber (1987) p.974.
c              write(*,*)'Tut 1'
                  if(mod(kk,2) .eq. 0) then
                      k = kk/2 + 1
                  else
                      k = ni+1 - (kk-1)/2
                  endif
                  r1 = r(k-1)                                 
                  a1 = a(k-1)                                 
                  b1 = b(k-1)                                 
                  x1 = r1*sin(a1)*cos(b1)                           
                  y1 = r1*sin(a1)*sin(b1)                           
                  z1 = r1*cos(a1)                                   
                  r3 = r(k+1)                                 
                  a3 = a(k+1)                                 
                  b3 = b(k+1)                                 
                  x3 = r3*sin(a3)*cos(b3)                           
                  y3 = r3*sin(a3)*sin(b3)                           
                  z3 = r3*cos(a3)                                   
                  dx = x3 - x1                                      
                  dy = y3 - y1                                      
                  dz = z3 - z1                                      
                  x2 = x1 + dx/2                                    
                  y2 = y1 + dy/2                                    
                  z2 = z1 + dz/2                                    
                  r2 = sqrt(x2**2 + y2**2 + z2**2)                  
                  acosa=z2/r2
                  if(acosa.LT.-1.)acosa=-1.
                  if(acosa.GT.1)acosa=1.
                  a2 = acos(acosa)                                  
                  sina = sin(a2)                                    
                  cosa = cos(a2)                                    
                  acosa=x2/r2/sina
                  if(acosa.LT.-1.)acosa=-1.
                  if(acosa.GT.1)acosa=1.
                  b2 = acos(acosa)                          
                  if(y.LT.0.00000)b2=360.00000*dpi-b2
c
c
                  dn = dx**2 + dy**2 + dz**2                        
                  ddn = sqrt(dn)                                    
                  dr = (r3-r1) / ddn                                
                  da = (a3-a1) / ddn                                
                  db = (b3-b1) / ddn

c  Begin find the gradients and velocities

c                 first find the length of segment
                  dseg=sqrt((dx/2)**2+(dy/2)**2+(dz/2)**2)

                  ddseg=dseg/2.
c                 write(*,*)'dseg ',dseg 

c                 Now ddseg will be a distance to find dV
c                 along the coordinates 

c                 Determine velocity at 3 points
c                 write(*,*)'before',ni
c                 write(*,*)
c           write(*,*)'Pered opredeleniem glubini ------------------'
c           do j=1,ni+1
c           write(*,*)j,ni,r(j),a(j)*r2d,b(j)*r2d
c           enddo
c                 write(*,*)'grohnulas pri',ni,' segmentov'
c                 write(*,*)'l,k ',l,k
c                 write(*,*)r1,a1*r2d,b1*r2d


                  v1 = vel2(isp,r1,a1,b1,1)                            

c                 write(*,*)'V1= ',v1
c                 write(*,*)r2,a2*r2d,b2*r2d

                  v2 = vel2(isp,r2,a2,b2,2)                            

c                 write(*,*)'V2= ',v2
c                 write(*,*)r3,a3*r2d,b3*r2d

                  v3 = vel2(isp,r3,a3,b3,3)                            

c                 write(*,*)'V3= ',v3

c                 Begin to determine coordinates
c                 of pints surroundibg point a2,b2,r2
c                 at the distance ddseg

                  upz = r2+ddseg
                  dwz = r2-ddseg

                   if(upz.gt.ro)then
                      upz=ro
                      dwz=upz-dseg
                   endif

                   if(dwz.le.0.)then
                      dwz=0.00000001
                      upz=r0
                   endif

C The following if-endif is just for P & S, thus comment out for SKS & PKP !!!
C This gives the lowermost mantle Vp in the outer core
                   if(dwz.le.3479.50)then
                      dwz=3479.500
                      upz=dwz+dseg
                   endif
C -------------------------------------------------


c                 find dV along depth
                  vr1 = vel2(isp,upz,a2,b2,4)
                  vr2 = vel2(isp,dwz,a2,b2,5)
                  vr=(vr1-vr2)/dseg

c                 write(*,*)'vr= ',vr

c                 find dV along longitude

                  call  km2deg(a2,b2,r2,ddseg,RNULL,adV,bdV,rdV) 
                  vb2 = vel2(isp,rdV,adV,bdV,6)
                  call  km2deg(a2,b2,r2,-1.*ddseg,RNULL,adV,bdV,rdV)
                  vb1 = vel2(isp,rdV,adV,bdV,7)
                  vb=-1.*(vb1-vb2)/dseg

c                 write(*,*)'vb= ',vb

c                 find dV along latitude

                  call  km2deg(a2,b2,r2,RNULL,ddseg,adV,bdV,rdV)
                  va2 = vel2(isp,rdV,adV,bdV,8)
                  call  km2deg(a2,b2,r2,RNULL,-1.*ddseg,adV,bdV,rdV)
                  va1 = vel2(isp,rdV,adV,bdV,9)
                  va=-1.*(va1-va2)/dseg

c                 write(*,*)'va= ',va
c                 write(*,*)'**************'
               
c                 write(*,*)'after',ni


c                 do j=1,ni+1
c                 write(*,*)j,ni,r(j),a(j)*r2d,b(j)*r2d
c                 enddo
      
			             
c               spherical
c                   velocity gradient
c                 va = va / r2
c                 vb = vb / r2 / sina
c                   (tangential vector) = (slowness vector) / s
c                 write(*,*)'perturbations !!!'
                  pr = dr
                  pa = r2 * da
                  pb = r2 * sina * db
c                 write(*,*)pr,pa,pb
                  vrd = pr*vr + pa*va + pb*vb
c                 write(*,*)vrd
                  rvr = vr - vrd*pr
                  rva = va - vrd*pa
                  rvb = vb - vrd*pb
c                 write(*,*)rvr,rva,rvb
                  rvs = sqrt(rvr*rvr + rva*rva + rvb*rvb)               
c                 write(*,*)rvs
c                 write(*,*)'posle rvs'
                  if(rvs .eq. 0.) then                              
                      r(k) = r2                                   
                      a(k) = a2                                   
                      b(k) = b2                                   
                  else                                              
                      rvr = rvr / rvs                                 
                      rva = rva / rvs                                 
                      rvb = rvb / rvs                                 
c                     write(*,*)rvr,rva,rvb
                      cc   = (1./v1+1./v3)/2.                          
c                     write(*,*)'cc ',cc
                      rcur = vr*rvr + va*rva + vb*rvb               
c                     write(*,*)'rcur ', rcur
                      if(rcur.LE.0.0)then
c                        write(*,*)'Negative'
                         rcur=abs(rcur)
                      endif
                      rcur = (cc*v2+1.) / (4.*cc*rcur)                
c                     write(*,*)'rcur ',rcur
                      rcur = -rcur + sqrt(rcur**2+dn/(8.*cc*v2))     
c                     write(*,*)'rcur ',rcur
                      rdr = rvr * rcur
                      rda = rva * rcur
                      rdb = rvb * rcur
c                     write(*,*)'rdr,rda,rdb ',rdr,rda,rdb
                      rrp  = r2 + rdr
                      ap  = a2 + rda/r2
                      bp  = b2 + rdb/(r2*sina)
c                     write(*,*)rrp,ap,bp
                      r(k) = (rrp-r(k))*xfac + r(k)
                      a(k) = (ap-a(k))*xfac + a(k)
                      b(k) = (bp-b(k))*xfac + b(k)
                  endif
c                 write(*,*)'after vector',ni
c                 write(*,*)r(k),a(k)*r2d,b(k)*r2d
c                 write(*,*)'++++++++++++++++++++++++++++++++++++++'
 2250             continue                                              
c                                                                       
              idstn=ni
              do 212  j=1,ni+1
                  w(1,j) = r(j)
                  w(2,j) = a(j)
                  w(3,j) = b(j)
c                 write(*,*)j,ni,w(1,j),w(2,j)*r2d,w(3,j)*r2d
 212          continue
              ni=idstn
c             here trace each iteration
c             write(*,*)'pered rtim',ni
              tk = rtim(isp,ni+1,r,a,b)
           
c              write(*,*)to,tk,ni
c              write(6,*) ni, l, (tk-atim)/atim*100.
              if(abs(to-tk) .le. to*flim)  go to 310                        
              to = tk                                                       
 3000     continue                                                          
 310      continue   

c          write(*,*)'loop=',loops, to, tk, abs(to-tk), to*flim
c          skip increasing of segment number if minimum length
c          of segment is exceed or maximum number of segments
c          was reached

           if(dseg.lt.mins.or.ni.ge.n2) go to 66666
c
c           double the number of points.
c           write(*,*)'double number',  to, ni
c
          to=tk
          ni = ni * 2
         
c         write(*,*)ni
          do i=1,ni/2+1
              r(i*2-1) = w(1,i)
              a(i*2-1) = w(2,i)
              b(i*2-1) = w(3,i)
          enddo
          do k=2,ni,2
              r1 = r(k-1)
              a1 = a(k-1)                                 
              b1 = b(k-1)                                 
              x1 = r1*sin(a1)*cos(b1)                           
              y1 = r1*sin(a1)*sin(b1)                           
              z1 = r1*cos(a1)                                   
              r3 = r(k+1)                                 
              a3 = a(k+1)                                 
              b3 = b(k+1)                                 
              x3 = r3*sin(a3)*cos(b3)                           
              y3 = r3*sin(a3)*sin(b3)                           
              z3 = r3*cos(a3)                                   
              dx = x3 - x1                                      
              dy = y3 - y1                                      
              dz = z3 - z1                                      
              x2 = x1 + dx/2                                    
              y2 = y1 + dy/2                                    
              z2 = z1 + dz/2                                    
              r2 = sqrt(x2**2 + y2**2 + z2**2)                  
              acosa=z2/r2
              if(acosa.LT.-1.)acosa=-1.
              if(acosa.GT.1)acosa=1.
              a2 = acos(acosa)
              sina = sin(a2)                                    
              acosa=x2/r2/sina
              if(acosa.LT.-1.)acosa=-1.
              if(acosa.GT.1)acosa=1.
              b2 = acos(acosa)
              if(y.LT.0.00000)b2=360.00000*dpi-b2
              r(k) = r2
              a(k) = a2
              b(k) = b2
          enddo
          tk = rtim(isp,ni+1,r,a,b)

66666     continue

          if(abs(to-tk) .le. to*flim)  go to 99999
          to = tk 
      enddo
c                                                                       
99999 continue

      tt = tk ! the final output of the 
      idstn=ni
      do i=1,ni+1
          w(1,i) = r(i)
          w(2,i) = a(i)
          w(3,i) = b(i)
      enddo
      ni=idstn
c Return coordinates to the origin
      idstn=ni
      do k=1,ni+1
          w(2,k) = w(2,k)*r2d
          w(3,k) = w(3,k)*r2d+shiftlo
          if(w(3,k).lt.0.)w(3,k)=360.+w(3,k)
c         write(*,*)k,ni,w(1,k),w(2,k),w(3,k)
      enddo 
      ni=idstn

c     write(*,*)'Calculated travel time for stop and seg:',tk,to,ni

c     do i=1,ni+1
c     write(*,*)w(1,i),w(2,i),w(3,i)
c     enddo
c     write(*,*)'Ray is done with ',ni,' segments'
      return
      end                                                               
c                                                                       
c
      function rlen(n, r, a, b)  
      implicit real*8 (a-h, o-z)                                        
      include 'RaySPDR.inc'
      integer n
      dimension  r(msg+1), a(msg+1), b(msg+1)
      rlen = 0.                                                         
      do 100  i=1,n-1                                                   
          j  = i + 1                                                    
          xi = r(i)*sin(a(i))*cos(b(i))                                 
          yi = r(i)*sin(a(i))*sin(b(i))                                 
          zi = r(i)*cos(a(i))                                           
          xj = r(j)*sin(a(j))*cos(b(j))                                 
          yj = r(j)*sin(a(j))*sin(b(j))                                 
          zj = r(j)*cos(a(j))                                           
          dl = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2                     
          rlen = rlen + sqrt(dl)                                        
  100 continue                                                          
      return                                                            
      end                                                               
c
c                                                                       
      function rtim(isp,m, r, a, b)             
      implicit real*8 (a-h, o-z)                                        
      include 'RaySPDR.inc'
      dimension  r(msg+1), a(msg+1), b(msg+1)
      integer m
      if(m.GT.(msg+1))write(*,*)'*'
      rtim = 0.
      rv1 = 1./vel2(isp,r(1),a(1),b(1),-1)
      do 100  j=1,m-1
          x1 = r(j)*sin(a(j))*cos(b(j))
          y1 = r(j)*sin(a(j))*sin(b(j))
          z1 = r(j)*cos(a(j))                                               
          x2 = r(j+1)*sin(a(j+1))*cos(b(j+1))
          y2 = r(j+1)*sin(a(j+1))*sin(b(j+1))
          z2 = r(j+1)*cos(a(j+1))
          dl = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
          rv2 = 1./vel2(isp,r(j+1),a(j+1),b(j+1),-2)
          sm = (rv1 + rv2) / 2.         
          rtim = rtim + sqrt(dl)*sm
          rv1 = rv2
  100 continue                                                          
      return                                                            
      end                                                               
c
c
      function vel2(isp,r, pa, ra, k)                              
      implicit real*8 (a-h, o-z)                                    
      common /coord/ shiftlo
      r2d = 90./asin(1.)                                                        

c     convert to degree and rotate coordinates 
      clat=pa*r2d
      clon=ra*r2d+shiftlo
      if(clon.lt.0.)clon=360.+clon

c     write(*,*)k
c     if(k.EQ.-1)write(*,*)clat,clon,r
      call cellvel(isp,clat,clon,r,V,nbl)
c     write(*,*)'exit vel'
      vel2=V

      return
      end


      subroutine cellvel(isp,la,lo,era,Vl,nbl)

      implicit real*8 (a-h, o-z)
      include 'RaySPDR.inc'

c     This subroutine determine the velocity in the
c     block using co Lat,Lon and R.
c     Lat and Lon in degrees
c     era - is the radius
c     dp is the depth

c     Note, the coordinates are colongitude and colatitude

      real*8 la,lo,era,dp
      integer  num
      dimension rak(npcom),vpak(npcom),vsak(npcom),dak(npcom)
      dimension xlayer(0:nlayer)
      real      x, y, z
      real      tlat, tlon, tdep
      real      x0,y0,z0,x1,y1,z1,x2,y2,z2
      real      v

      common /model/   xlayer,rak,vpak,vsak,dak

      parameter (nxi=72,nyi=36,nzi=16)
      parameter (dx=5.,dy=5.)

      dpi = asin(1.)/ 90.
      r2d = 90./asin(1.)

c     Firstly determine number of layer

      Re=6371.00
      dp=Re-era ! depth of the point
     
      nlr=0
      nla=0
      nlo=0

      if(dp.LT.0.00)then
        dp=0.0
      endif

      do while(lo.GE.360.)
         lo=lo-360.
      enddo

      do while(lo.LT.0.0e10)
         lo=lo+360.00
      enddo

c---- change co-latitude back to latitude
  
      la=90-la ! it should be from -90 to 90
c--- Now need to consider the ellipticity correction
      !atemp = 1.0/.99664719;
      atemp = 1.0;
      la=atan(atemp*tan(la*dpi))*r2d

c--- judging the layer fro AK135 model
      do i=2,np15
 
       if(dp.GE.dak(i-1).AND.dp.LT.dak(i))then
           ilr=i-1   ! calculate the layer number for the AK135 model
           go to 5
       endif    

      enddo 
5     continue

****  IF   *******************************************************
c     Where is the point - inside/outside of the region

      if(lo.GT.xn(1)+0.1.AND.lo.LT.xn(nx)-0.1.AND.la.GT.yn(1)+0.1
     &   .AND.la.LT.yn(ny)-0.1)then
      
c      inside of the study region

      if(dp.GE.zn(nz)-1)then
c        write(*,*)'Velocity determination
c     &   is deeper than regional nodes go to the global scale'
         go to  24
      endif


c--- decide the index number of the box
c  use Prothero's intmap here
      x = lo
      y = la
      z = dp

      !write(*,*) "in the VEL3:",x,y,z,xn(1),xn(nx),yn(1),yn(ny),zn(1),zn(nz)
      call vel3(isp,x,y,z,v)
      
      !write(*,*)'Random velocity:',vel(3,4,5),vel(7,8,9) 
      !write(*,*)"In cellvel",isp,x,y,z,v
      VL=v

      if(isp.eq.1) kp=kp+nz

      nbl=ip+(jp-1)*nx+(kp-1)*nx*ny

      else
c           off the region of study
c           now search for the larger block from the global model

24    continue

c     write(*,*)'outside the small region!'

       ip = 0
       jp = 0
       kp = 0
c--- just divide the Earth into many blocks

      do i=1,nzi

       if(dp.GE.xlayer(i-1).AND.dp.LT.xlayer(i))then
           nlr=i
           go to 1
       endif

      enddo
1     continue

      do i=1,nxi

       if(lo.GE.real(dx*real(i-1)).AND.
     &    lo.LT.real(dx*real(i)))then
           nlo=i
           go to 2
       endif

      enddo
2     continue

c     and then Row (lat) number

c--- note here la is latitude, not co-latitude
      do i=1,nyi

       if(la.GE.dy*real(i-1).AND.la.LT.dy*real(i))then
           nla=i
           go to 3
       endif

      enddo
3     continue

      num=(nla-1)*nxi+nlo+(nlr-1)*nxi*nyi+nx*ny*nz
      
      if(isp.eq.0) VL=vpak(ilr)
      if(isp.eq.1) Vl=vsak(ilr)

      nbl=num

      endif

      return
      end

      subroutine cderiv(isp,w,ni,dth2,dth3,dth4)
      implicit real*8 (a-h, o-z)
      parameter (msg = 512)
      dimension  w(3,msg+1)
      dpi = asin(1.)/ 90.
      r2d = 90./asin(1.)
      r0        = 6371.00
      epc       = 1.0e-5
      nrp       = ni+1
      pe        = w(2,nrp)*dpi
      re        = w(3,nrp)*dpi
      he        = w(1,nrp)
      nrp1      = nrp-1
      pe1       = w(2,nrp1)*dpi
      re1       = w(3,nrp1)*dpi
      he1       = w(1,nrp1)
      call conaz(re,pe,re1,pe1,des,az)
      des       = des*(he/r0)
      dhe       = he1-he
      adh       = abs(dhe)
      call cellvel(isp,pe*r2d,re*r2d,he,ve,nbl)
      if(des.lt.epc.and.adh.lt.epc) dhe = 1.0
      th        = atan2(des,dhe)
      dtddel    = (sin(th)*he/ve)
c     write(*,*)'rp= ',dtddel
      dtdr      = abs(cos(th)/ve)
c     write(*,*)'az= ',az*r2d
      if(az*r2d.GT.360.)az=az-360.*dpi
      dth3 =-dtddel*cos(az)
      dth4 = dtddel*sin(az)*sin(pe)
      dth2 =-sqrt((he/ve)**2-dtddel**2)/he
      return
      end

      subroutine conaz(re,pe,rs,ps,del,az)
      implicit real*8 (a-h, o-z)
      data r0,eps,pi15/6371.00,1.0e-7,4.712389/
      sps  = sin(ps)
      cps  = cos(ps)
      spe  = sin(pe)
      cpe  = cos(pe)
      ses  = sin(re-rs)
      ces  = cos(re-rs)
      x    = sps*ses
      y    = cpe*sps*ces - spe*cps
      s    = sqrt(x*x + y*y)
      c    = spe*sps*ces + cpe*cps
      del  = atan2(s,c)*r0
      az   = 0.0
      ax   = abs(x)
      ay   = abs(y)
      if(ax.lt.eps.and.ay.lt.eps) return
      az   = pi15-atan2(y,x)
      return
      end

cz- This file is taken from simul2000 algorithm by Clifford Thurber

cz--- now x is longitude, y is latitude and z is the depth

      subroutine vel3(isp,x,y,z,v)
c  This routine is Cliff Thurber's
c  common block variables:
      common/weight/ wv(8),ip,jp,kp,kpg
      include 'RaySPDR.inc'
c
c  use Prothero's intmap here
      call intmap(x,y,z,ip,jp,kp)
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
c     write(16,100)x,xl,ip,y,yl,jp,z,zl,kp
c      if(ip.eq.0 .or. jp.eq.0 .or. kp.eq.0) then
c         write(*,*)x,xl,ip,y,yl,jp,z,zl,kp
c      endif
c100	format(3(2f7.3,i3))
      xf=(x-xn(ip))/(xn(ip1)-xn(ip))
      yf=(y-yn(jp))/(yn(jp1)-yn(jp))
      zf=(z-zn(kp))/(zn(kp1)-zn(kp))
      xf1=1.0-xf
      yf1=1.0-yf
      zf1=1.0-zf
c
      wv(1)=xf1*yf1*zf1
      wv(2)=xf*yf1*zf1
      wv(3)=xf1*yf*zf1
      wv(4)=xf*yf*zf1
      wv(5)=xf1*yf1*zf
      wv(6)=xf*yf1*zf
      wv(7)=xf1*yf*zf
      wv(8)=xf*yf*zf
c     write(*,*)wv(1),wv(2),wv(3),wv(4),wv(5),wv(6),wv(7),wv(8)
c  calculate velocity
c  S-velocity is stored after P-velocity
c  (or V*Q if iuseq=1)
      kpg=kp
      if(isp.eq.1) kp=kp+nz
      kp1=kp+1
      v=wv(1)*vel(ip,jp,kp)+wv(2)*vel(ip1,jp,kp)
     2 +wv(3)*vel(ip,jp1,kp)+wv(4)*vel(ip1,jp1,kp)
     * +wv(5)*vel(ip,jp,kp1)+wv(6)*vel(ip1,jp,kp1)
     * +wv(7)*vel(ip,jp1,kp1)+wv(8)*vel(ip1,jp1,kp1)
       !write(*,*)wv(1),wv(2),wv(3),wv(4),wv(5),wv(6),wv(7),wv(8),v
      return
c***** end of subroutine vel3 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine intmap(x,y,z,ip,jp,kp)
c  Modified by W. Prothero so a single call can get the indices
c  common block variables:
      include 'RaySPDR.inc'
c
      ip=int((x+xl)/bld)
      ip=ixloc(ip)
      jp=int((yl+y)/bld)
      jp=iyloc(jp)
      kp=int((z+zl)/bld)
      kp=izloc(kp)
c  If an array element=0, the position is off the map.
      return
c***** end of subroutine intmap *****
      end


c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c

      subroutine ttmder(isp, dth1, stepl,wlat,wlon,theta,tt)
c
c  declaration statements:
      dimension in(8)     
c
c  common block variables:
      real dth1(4),stepl,wlat,wlon,theta 
      integer isp
      real xe,ye,ze,xc,yc,zc,tlat,tlon,tdep,v
      real dx, dy, dz, ds, uds, tt
      common/weight/ wv(8),ip,jp,kp,kpg
      include 'RaySPDR.inc'

c  resegment ray path
c  calculate travel time derivatives with respect to hypocentral
c  parameters - use coords of first two points on raypath
c  determine slowness at source

c--- now the ray path is stored in X-Y-Z 

      xc=rp(1,1)
      yc=rp(2,1)
      zc=rp(3,1)

c--- change it into lat-lon-dep to find velocity
      call car2sph_ft(xc,yc,zc,wlat,wlon,tlat,tlon,tdep,theta)

      !write(*,*)tlat,tlon,tdep,xl,yl,zl
      call vel3(isp,tlon,tlat,tdep,v)
             
      us=1.0/v

c**
c  determine cartesian derivatives from direction cosines

      dx=rp(1,2)-rp(1,1)
      dy=rp(2,2)-rp(2,1)
      dz=rp(3,2)-rp(3,1)

      ds=sqrt(dx*dx+dy*dy+dz*dz)
      uds=-us/ds
c  hypocentral derivatives
      dth1(2)=uds*dx
      dth1(3)=uds*dy
      dth1(4)=uds*dz
c  origin time derivative
      dth1(1)=1.0

c  skip next section if all velocity nodes are fixed (nparvi=0)
      if(nparvi.eq.0) goto 88
c  skip next section if only doing earthquake location
c     if(kopt.eq.0) go to 88
c
c  travel time and velocity partial derivatives
      tt=0.0
      half=0.5
cz  initialize the velocity derivative
      do i=1, npari
         dtm(i)=0.0
      enddo
c  loop over segments comprising the ray path
c     write(*,*)'Number of segments and stepl:',nrp, stepl
      nrp1=nrp-1
      pl=0.0
      do 50 i=1,nrp1
         i1=i+1
         rx=rp(1,i)
         ry=rp(2,i)
         rz=rp(3,i)
         dx=rp(1,i1)-rx
         dy=rp(2,i1)-ry
         dz=rp(3,i1)-rz
c  compute segment length
         sl=sqrt(dx*dx+dy*dy+dz*dz)
         pl=pl+sl

c  decide on number of subsegments and compute length
         nseg=nint(sl/stepl)+1
         fnsegi=1.0/float(nseg)
         ssl=sl*fnsegi
         dxs=dx*fnsegi
         dys=dy*fnsegi
         dzs=dz*fnsegi
c-----
         xp=rx-half*dxs
         yp=ry-half*dys
         zp=rz-half*dzs

c  loop over subsegments
         do 55 is=1,nseg
            xp=xp+dxs
            yp=yp+dys
            zp=zp+dzs
c now xp, yp, zp are X-Y-Z, change them into lat-lon-dep 
           
            call car2sph_ft(xp,yp,zp,wlat,wlon,tlat,tlon,tdep,theta)

c--- now it is possible that the ray may be outside the effective domain
c--- first determine if tlat, tlon and tdep are in the domain.

            if( tlon.le.(xn(1)+0.1) .or. tlon.ge.(xn(nx)-0.1) .or.
     &          tlat.le.(yn(1)+0.1) .or. tlat.ge.(yn(ny)-0.1) .or.
     &          tdep.le.(zn(1)+1.0) .or. tdep.ge.(zn(nz)-1.0) ) goto 1717

            call vel3(isp, tlon, tlat, tdep, v)
            dt=ssl/v
            tt=tt+dt
c
c
c  The next section is a change from 'block' to 'linear'
c   partial derivatives, by C.Thurber,may10,1983.
c  Nodes with non-zero weight
            in(1)=ip-1+nx2*(jp-2)+nxy2*(kp-2)-nxy2*(2*isp)
            in(2)=in(1)+1
            in(3)=in(1)+nx2
            in(4)=in(3)+1
            in(5)=in(1)+nxy2
            in(6)=in(5)+1
            in(7)=in(5)+nx2
            in(8)=in(7)+1
c
c  Assign zero weight to boundary nodes (these nodes are not 
c  included in the inversion, but are in the velocity array,
c  thus we want to avoid writing to negative or incorrect 
c  elements of the partial derivative matrix)

            if(ip.eq.1) then
c              write(16,1610) xp,yp,zp,v,ip,jp,kpg
               wv(1)=0.0
               wv(3)=0.0
               wv(5)=0.0
               wv(7)=0.0
            else
               if(ip.eq.nx1) then
c                 write(16,1610) xp,yp,zp,v,ip,jp,kpg
                  wv(2)=0.0
                  wv(4)=0.0
                  wv(6)=0.0
                  wv(8)=0.0
               end if
            endif

            if(jp.eq.1) then
c              write(16,1610) xp,yp,zp,v,ip,jp,kpg
               wv(1)=0.0
               wv(2)=0.0
               wv(5)=0.0
               wv(6)=0.0
            else
               if(jp.eq.ny1) then
c                 write(16,1610) xp,yp,zp,v,ip,jp,kpg
                  wv(3)=0.0
                  wv(4)=0.0
                  wv(7)=0.0
                  wv(8)=0.0
               endif
            endif

            if((kpg.eq.1).or.(kpg.eq.(nz1+1))) then
c              write(16,1610) xp,yp,zp,v,ip,jp,kpg
               do 30 izg=1,4
                  wv(izg)=0.0
   30          continue
            else
               if((kpg.eq.nz1).or.(kpg.eq.(2*nz1))) then
c                 write(16,1610) xp,yp,zp,v,ip,jp,kpg
                  do 35 izg=5,8
                     wv(izg)=0.0
   35             continue
               endif
            endif
 1610       format(' ASSIGNING ZERO WEIGHTS IN TTMDER',
     2         'xp=',f7.2,',yp=',f7.2,',zp=',
     3         f7.2,',v=',f5.3,',ip=',i2,',jp=',i2,',kpg=',i2)
c
c  Accumulate model partial derivatives
            do 48 kk=1,2
               kk1=kk-1
               do 47 jj=1,2
                  jj1=jj-1
                  do 46 ii=1,2
                     ii1=ii-1
                     ijk=ii+2*jj1+4*kk1
c skip boundary nodes
                     if(wv(ijk).lt.0.05) goto 46
c DEP
c write out DWS for all nodes (including fixed and linked) when nitmax=1
c (useful for planning fixed and linked)
c Include weight factor like for inversion
c Note that for shots on nit=0, combined weight, including residurpal
c wt'g is not yet calculated so use wtsht.
cz--- Keep it unused temperarily ( by H. Zhang)
cz             if(nitmax.eq.1) then
cz               if((ne.gt.(neqs+nbls)).and.(nit.eq.0)) then
cz                 hitall(in(ijk))=hitall(in(ijk))+wv(ijk)*wtsht
cz               else
cz                 hitall(in(ijk))=hitall(in(ijk))+wv(ijk)*wtcomb(no,ne)
cz               endif
cz             endif
c  skip fixed nodes
                     if(nfix(in(ijk)).eq.1) goto 46
                     ini=ndexfx(in(ijk))
c
c  start cht 1998
                     if (imerge(in(ijk)).eq.1) ini=ndexfx(jequal(in(ijk)))
c  end cht 1998
c
c check for writing to an array element that is outside of the inversion
c  solution array
                     if((ini.lt.1).or.(ini.gt.nparvi)) then
                        write(16,1606) ini,ijk
 1606                   format(' *** Error in TTMDER, accessing',
     2                     ' gridpoint outside of velocity inversion',
     3                     ' gridpoints, ini=',i5,', ijk=',i5,/,
     4                     22x,'Probably boundary gridpoints are',
     5                     ' too close (TTMDER tries to write DTM',
     6                     ' elements with wv >= 0.05)')
                        write(16,1603) ne,no,xp,yp,zp,v,ip,jp,kp,kpg
 1603                   format(' ne=',i5,', no=',i5,', xp=',f8.2,
     2                     ', yp=',f8.2,', zp=',f8.2,', v=',f8.3,/,
     3                     21x,'ip=',i6,',   jp=',i6,',   kp=',i6,
     4                     ',   kpg=',i6)
                        write(16,1607) (j,in(j),j,wv(j),j=1,8)
 1607                   format(' in(',i1,')=',i6,' wv(',i1,')=',e15.5)
                        write(16,1608)
 1608                   format(' * * * * STOP * * * * (to avoid',
     2                     ' writing outside of defined DTM array)')
                        stop
                     end if

                     inp=ini
cDEP Now include weight factor for hit(DWS)
cz                     if((ne.gt.(neqs+nbls)).and.(nit.eq.0)) then
cz                       hit(in(ijk))=hit(in(ijk))+wv(ijk)*wtsht
cz                     else
cz                       hit(in(ijk))=hit(in(ijk))+wv(ijk)*wtcomb(no,ne)
cz                     endif
cz                     hit(in(ijk))=hit(in(ijk))+wv(ijk)*TimeWeight

c     set up the model derivative matrix
                     dtm(inp)=dtm(inp)+wv(ijk)*ssl
                     
c***
 46               continue
 47            continue
 48         continue
c
 55      continue
 50   continue
1717  continue
      !write(*,*)'tt from ttmder:',tt
 88   continue
c**
c
      return
c***** end of subroutine ttmder *****
      end
c

c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c

      subroutine bldmap
c  common block variables:
      include 'RaySPDR.inc'
c
c     array size limits
c     ixkms=iykms=izkms=1500
c
c     write(6,400)
c     400 format(' subroutine bldmap')
      xl=bld-xn(1)
      ixmax=(xn(nx)+xl)/bld
      yl=bld-yn(1)
      iymax=(yn(ny)+yl)/bld
      zl=bld-zn(1)
      izmax=(zn(nz)+zl)/bld
c     write(6,402)ixmax,iymax,izmax
c 402 format(' array sizes: ',3i5)
c
c  Check for array size overflow
      if(ixmax.gt.ixkms.or.iymax.gt.iykms.or.izmax.gt.izkms)goto 330
      ix=1
      do 10 i=1,ixmax
c
         ix1=ix+1
c
         xnow=float(i)*bld-xl
         if (xnow.ge.xn(ix1)) ix=ix1
c
         ixloc(i)=ix
 10   continue
c  Fill remainder of array with zeroes.
      do 12 i=ixmax,ixkms
         ixloc(i)=0
 12   continue
c
c
      iy=1
      do 15 i=1,iymax
c
         iy1=iy+1
c
         ynow=float(i)*bld-yl
         if (ynow.ge.yn(iy1)) iy=iy1
c
         iyloc(i)=iy
 15   continue
c
c     Fill rest of array with zeroes.
      do 17 i=iymax,iykms
         iyloc(i)=0
 17   continue
c
      iz=1
      do 20 i=1,izmax
c
         iz1=iz+1
c
         znow=float(i)*bld-zl
         if (znow.ge.zn(iz1)) iz=iz1
c
         izloc(i)=iz
 20   continue
c
c     Fill remainder of array with zeroes.
      do 22 i=izmax,izkms
         izloc(i)=0
 22   continue
      return
 330  continue
      write(16,331)ixkms,iykms,izkms
 331  format(' ***** error in array size in common/locate/',/,
     *     ' maximum map dimensions (km)=',/,' x=',i5,' y=',i5,' z=',i5)
      write(16,332)ixmax,iymax,izmax
 332  format(' Actual map size (km): ',/,' x=',i5,' y=',i5,' z=',i5)
      stop
c*****end of subroutine bldmap *****
      end

      subroutine input_vel
c     this routine reads in the initial velocity model in the
c     form of velocity specified on a uniform but not
c     necessarily evenly spaced grid of points
c     (reads from file03 )
c
c  common block variables:
      include 'RaySPDR.inc'
c
c  declaration statements:
      integer ixf(maxpar),iyf(maxpar),izf(maxpar)
c
c     start cht 1998
      
c     end cht 1998
      character*1 vtype(2)
      parameter(zero=0.0,izero=0)
c
      write(16,*)'iuses=',iuses
      write(16,*)'invdel=',invdel
      write(16,*)'iuseq=',iuseq

      do n=1,maxpar
         cnode(n)='0'
      enddo

      ierror=0
      vtype(1)='P'
      vtype(2)='S'
c
c     for this version the gridpoints can be unevenly spaced
c     the origin of the coordinate system is at (x,y,z)=(0,0,0)
c     which will not in general correspond to the point
c     xn(1),yn(1),zn(1).
c     xn,yn,zn should be factors of bld (ie: a.0 for bld=1.0 or a.b for bld=0.1)
c
c     input the number of gridpoints in x, y and z directions
c     and bld factor (1.0 or 0.1 km) used to set up velocity interpolation grid
c     read(3,3002) bld,nx,ny,nz
      read(3,*) bld, nx, ny, nz
 3002 format(f4.1,3i3)
      if((bld.ne.1.0).and.(bld.ne.0.1).and.(bld.ne.0.01).and.(bld.ne.0.001)) then
         write(16,1625) bld
 1625    format(/, '******** STOP *********, bld must be 1.0 or 0.1,
     &   not ',f6.2)
      endif
      atemp=iuses*(nx-2)*(ny-2)*(nz-2)
      if(atemp.le.maxpar)goto 40
      write(16,42)
 42   format('0Too many nodes for program array sizes.')
      stop
 40   continue
c
c     start cht 1998
      do 123 k=1,atemp
         imerge(k)=0
         jequal(k)=0
 123  continue
c
c  end cht 1998
c
c  input the x grid, y grid, and z grid
cfh read in free format (makes life easier...)
c     read(3,3004) (xn(i),i=1,nx)
c     read(3,3004) (yn(i),i=1,ny)
c     read(3,3004) (zn(i),i=1,nz)
      read(3,*) (xn(i),i=1,nx) ! longitude
      read(3,*) (yn(i),i=1,ny) ! latitude
      read(3,*) (zn(i),i=1,nz) ! depth
 3003 format(3i3)
 3004 format(20f6.1)
c
      write(16,3005) bld,nx,ny,nz
 3005 format(//,' velocity grid size:',/,
     *     'bld =',f4.1,5x,' nx =',i3,5x,'ny =',i3,5x,'nz =',i3)
c
cfh give all these numbers the same format
      write(16,3006) (xn(i),i=1,nx)
cfh 3006 format(/,' xgrid',/,3x,12f7.1,8f6.1)
 3006 format(/,' xgrid',/,3x,20f7.1)
      write(16,3007) (yn(i),i=1,ny)
cfh 3007 format(/,' ygrid',/,3x,12f7.1,8f6.1)
 3007 format(/,' ygrid',/,3x,20f7.1)
      write(16,3008) (zn(i),i=1,nz)
cfh 3008 format(/,' zgrid',/,3x,8f6.1,12f7.1/)
 3008 format(/,' zgrid',/,3x,20f7.1/)
c
c  set all nodes to have fixed velocity - cht 2002
cz this set is removed by HZ
      inf=0


c
c  start cht 1998
c  lines moved followed by new code
c  compute total number of gridpoints (nodes)
      nodes=nx*ny*nz
      nxy=nx*ny
      nx2=nx-2                  ! number non-edge nodes in row
      nxy2=nx2*(ny-2)           ! number non-edge nodes in layer
      nz2=nz-2
      nodes2=nz2*nxy2
c  peripheral nodes
      nx1=nx-1
      ny1=ny-1
      nz1=nz-1
c
c
c  end cht 1998
c
c  now read in the velocity values
 65   write(16,3101)
c     do 38 kv=1,iuses
         kv=1
         do 37 k=1,nz
            k2=k + (kv-1)*nz
            write(16,3015) k,vtype(kv),zn(k)
            do 36 j=1,ny
cfh            read(3,3011) (vel(i,j,k2),i=1,nx)
               read(3,*) (vel(i,j,k2),i=1,nx)
               write(16,3013) (vel(i,j,k2),i=1,nx)
 36         continue
 37      continue
c38   continue
c CHANGE FOR VP/VS INVERSION
         if((iuses.eq.2).or.(iuseq.eq.1)) then
            do 100 k=1,nz
               if(iuseq .eq. 0) then
                  write(16,3016) k,zn(k)
               else
                  write(16,3017) k,zn(k)
               endif
               do 99 j=1,ny
                  if(iuseq .eq. 0) then
cfh                 read(3,3011) (vpvs(i,j,k),i=1,nx)
                     read(3,*) (vpvs(i,j,k),i=1,nx)
                     write(16,3013) (vpvs(i,j,k),i=1,nx)
                  else
                     read(3,3014) (qval(i,j,k),i=1,nx)
                     write(16,3014) (qval(i,j,k),i=1,nx)
                  endif
 99            continue
 100        continue
c  compute Vs from Vp and Vp/Vs or compute 1/tstar
            kv=2
            if(iuseq .eq. 0) then
               do 120 k=1,nz
                  ks=k+nz
                  write(16,3015) k,vtype(kv),zn(k)  
                  do 115 j=1,ny
                     do 110 i=1,nx
                        vel(i,j,ks)=vel(i,j,k)/vpvs(i,j,k)
 110                 continue
                     write(16,3013) (vel(i,j,ks),i=1,nx)
 115              continue
 120           continue
            else
               do 140 k=1,nz
                  ks=k+nz
                  write(16,3018) k,zn(k)
                  do 135 j=1,ny
                     do 130 i=1,nx
                        vel(i,j,ks)=vel(i,j,k)*qval(i,j,k)
 130                 continue
                     write(16,3014) (vel(i,j,ks),i=1,nx)
 135              continue
 140           continue
            endif
         endif
c
 3013    format(20f6.2)
 3014    format(20f7.1)
 3015    format(/,' layer',i3,5x,a1,' velocity',10x,'z =',f7.1)
 3016    format(/,' layer',i3,5x,'Vp/Vs',10x,'z =',f7.1)
 3017    format(/,' layer',i3,5x,'Q',10x,'z =',f7.1)
 3018    format(/,' layer',i3,5x,'Q * Vp',10x,'z =',f7.1)
 3011    format(20f5.2)
 3101    format(//,' velocity values on three-dimensional grid')
c  Number of medium parameters to invert for
         npar=nodes2*iuses
         nparv=npar
c     if(invdel.ne.0)npar=(npar+nsts*iuses)
c
c  Check to see whether medium parameters fit within array sizes
         if(nparv.gt.maxpar) goto 980
c
cfhdmep
c get number of Vp and Vp/Vs nodes that are free in the inversion
c    nodes reduced by fixednodes
         nparpi=nodes2
         nparsi=nodes2
c  fix specified nodes by setting nfix(k)=1, else=0
         if(inf.eq.0) goto 496
         do 70 i=1,inf
            iizf=izf(i)-2
c  if s velocity node
c     fhdmep           if(izf(i).gt.nz) iizf=izf(i)-4
            if(izf(i).gt.nz) then
               iizf=izf(i)-4
               nparsi=nparsi-1
            else
               nparpi=nparpi-1
            endif
c
            k=iizf*nxy2 + (iyf(i)-2)*nx2 + (ixf(i)-1)
            nfix(k)=1
            if(cnode(k).ne.'0') then
               write(16,1686) cnode(mnode)
 1686          format(' *-*-* ERROR velocity input.  This node has',
     2  ' already been ',/,' *-*-* assigned cnode= ',a1)
               ierror=1
               write(16,1681) ixf(i),iyf(i),izf(i),k
 1681          format('input3 fixed node',i5,' ixf,iyf,izf:',3i3,
     2       ' node number:',i5)
            endif
           cnode(k)='F'
 70     continue
c
        write(16,1611)
 1610   format(/,' velocity FIXED at the following nodes(1):')
 1611   format(/,' VELOCITY INVERSION GRID    0=free, F=Fixed',/,
     2  '   M=Master, C=Constant Pert. Link, ',
     3  'L=Linear Pert. Link')
 311    do 495 kv=1,iuses
           nz1=nz-1
           ny1=ny-1
           do 320 k=2,nz1
              if(iuseq.eq.0) then
                 if(kv.eq.1) write(16,1009) k,vtype(kv),zn(k)
 1009            format(/,' layer',i3,5x,a1,'-velocity nodes',
     2                10x,'z =',f7.1)
                 if(kv.eq.2) write(16,3016) k,zn(k)
              else
                 write(16,3017) k,zn(k)
              endif
              kk=k+(kv-1)*nz2

              do 310 j=2,ny1
                 n1=(kk-2)*nxy2+(j-2)*nx2+1
                 n2=n1+nx2-1
                 write(16,1006) (cnode(i),i=n1,n2)
c     write(16,1005) (nfix(i),i=n1,n2)
 310          continue
 320       continue
c 1005 format('    ',18i6)
 1006      format('  ',40(2x,a1))
 495    continue
 496    continue
c
c  ndexfx: index from full nodes to nodes reduced by fixed (invert nodes)
c  mdexfx: index from inversion solution nodes to full velocity nodes
        in=0
c
c  start cht 1998
        infl=inf+ilink
c
        write(16,6161) inf,ilink,infl
 6161 format(/,' number of fixed, linked, fixed+linked nodes: ',3i5)
c
      do 80 n=1,nparv
c
c  remove fixed and linked nodes from inversion solution node set
c
         if(nfix(n).eq.1) goto 80
c  start of imerge if-then-else
         if(imerge(n).eq.0) go to 888
c
c  calculate x-y-z indices of linked velocity grid
         k=(n-1)/nxy2+2
         j=2+(n-1+(2-k)*nxy2)/nx2
         i=1+n+nx2*(2-j)+nxy2*(2-k)
         if(k.ge.nz) k=k+2      ! if s velocity node
c
c  calculate x and z indices of "master" velocity grid
         nma=jequal(n)
c
         km=(nma-1)/nxy2+2
         jm=2+(nma-1+(2-km)*nxy2)/nx2
         im=1+nma+nx2*(2-jm)+nxy2*(2-km)
         if(km.ge.nz) km=km+2   ! if s velocity node
c
cDEP It seems better to have "constant" link for velocity derivatives
cDEP  to result in constant perturbation between master and linked,
cDEP  rather than just constant velocity.  So do not change initial
cDEP  model here. Of course a constant velocity initial model could
cDEP  be input if that was desired.
cDEPc  start of constant/linear link if-then-else
cDEP      if (ltype(nma).eq.1) then
cDEP      write(16,2626) i,j,k,n,im,jm,km,nma
cDEP 2626 format('  setting value for constant-type linked node: ',
cDEP     &  4i4,2x,4i4)
cDEPc
cDEP       if(iuseq.eq.0) then
cDEP           vel(i,j,k)=vel(im,jm,km)
cDEP           if (k.gt.nz)
cDEP     &     vpvs(i,j,k-nz)=vpvs(im,jm,km-nz)
cDEPc
cDEP       else
cDEP           vel(i,j,k+nz)=vel(im,jm,km+nz)
cDEP       endif
cDEPc
cDEP      else
cDEPc   continuing constant/linear link if-then-else
cDEP       il=i
cDEP       jl=j
cDEP       ktemp=k
cDEP       kmtemp=km
cDEP      if (iuseq.ne.0) then
cDEP       k=k+nz
cDEP       km=km+nz
cDEP      endif
cDEP       kl=k
cDEPc
cDEP       if (i.ne.im) then
cDEP       il=2*i-im
cDEP       vel(i,j,k)=vel(im,jm,km)+(vel(il,jl,kl)-vel(im,jm,km))*
cDEP     & (xn(i)-xn(im))/(xn(il)-xn(im))
cDEP       endif
cDEPc
cDEP       if (j.ne.jm) then
cDEP       jl=2*j-jm
cDEP       vel(i,j,k)=vel(im,jm,km)+(vel(il,jl,kl)-vel(im,jm,km))*
cDEP     & (yn(i)-yn(im))/(yn(il)-yn(im))
cDEP       endif
cDEPc
cDEP       if (k.ne.km) then
cDEP       kl=2*k-km
cDEP       vel(i,j,k)=vel(im,jm,km)+(vel(il,jl,kl)-vel(im,jm,km))*
cDEP     & (zn(i)-zn(im))/(zn(il)-zn(im))
cDEP       endif
cDEPc
cDEP      write(16,2627) i,j,k,n,im,jm,km,nma,il,jl,kl,vel(i,j,k)
cDEP 2627 format('  setting value for linear-type linked nodes: ',/,
cDEP     &  4i4,2x,4i4,2x,3i4,f7.2)
cDEPc
cDEP      if (k.gt.nz) then
cDEP       if (i.ne.im) then
cDEP       il=2*i-im
cDEP       vpvs(i,j,k)=vpvs(im,jm,km)+(vpvs(il,jl,kl)-vpvs(im,jm,km))*
cDEP     & (xn(i)-xn(im))/(xn(il)-xn(im))
cDEP       endif
cDEPc
cDEP       if (j.ne.jm) then
cDEP       jl=2*j-jm
cDEP       vpvs(i,j,k)=vpvs(im,jm,km)+(vpvs(il,jl,kl)-vpvs(im,jm,km))*
cDEP     & (yn(i)-yn(im))/(yn(il)-yn(im))
cDEP       endif
cDEPc
cDEP       if (k.ne.km) then
cDEP       kl=2*k-km
cDEP       vpvs(i,j,k)=vpvs(im,jm,km)+(vpvs(il,jl,kl)-vpvs(im,jm,km))*
cDEP     & (zn(i)-zn(im))/(zn(il)-zn(im))
cDEP       endif
cDEP      endif
cDEPc
cDEP      endif
c
         go to 80
c     
c  end of linked node section
c
 888     continue
c  add to index if not linked or fixed
c
         in=in+1
         ndexfx(n)=in
         mdexfx(in)=n
c         write(16,1888) in,n
c 1888    format(' adding inversion node number ',
c     &   i5,' corresponding to gridpoint number ',i5)
c
 80   continue
      inf2=nparv-in
      if(inf2.eq.infl) goto 85
c
c  end cht 1998
c
      write(16,1615) infl,nparv,in,inf2
 1615 format(/,' **** number of fixed and linked nodes input,',i4,
     2  ', does not equal velocity nodes,',i4,', minus invert',
     3  ' nodes,',i4,'.  Continue with inf=',i5,' ****',/)
      infl=inf2
 85   continue
      nparvi=nparv-infl
      npari=nparvi
c     if(invdel.ne.0) npari=nparvi+nstsi*iuses
      if(invdel.eq.0) goto 95
c  also set indices if station delays are included in inversion
      i1=nparv+1
      do 90 i=i1,npar
         is=i-nparv
c s-delay
cz         if(is.gt.nsts) is=is-nsts
         if(nfixst(is).eq.1) goto 90
         in=in+1
         ndexfx(i)=in
         mdexfx(in)=i
 90   continue
      npari=in
 95   continue
      write(16,1620) npar,nparv,npari,nparvi
 1620 format(' INPUT3:npar,nparv,npari,nparvi',4i6)
c  Check to see whether medium parameters fit within array sizes
      if(npari.gt.mxpari) goto 990
      if(npar.gt.maxpar) goto 995
c  Stop if there was an error reading in fixed and linked nodes
      if(ierror.ne.0) then
         write(16,1683)
         write(6,1683)
 1683    format(/,'STOP SIMUL2000! , Error in velocity input file ')
         stop
      endif
c
c  Set up an array which is used to point to node indices, for any x,y,z
c     write(*,*)'nparsi=',nparsi
      call bldmap
c     write(*,*)'nparsi=',nparsi
c
      return
c
 980  continue
      write(6,1683)
      write(16,1698) nparv,maxpar
 1698 format(/,'  ****** STOP ******',/,i8,' velocity nodes, program',
     2 ' arrays only allow',i6)
      stop
 990  continue
      write(6,1683)
      write(16,1699) npari,mxpari
 1699 format(/,'  ****** STOP ******',/,i8,' parameters to invert for',
     2     ', program arrays only allow',i6)
      stop
 995  continue
      write(6,1683)
      write(16,1695) npar,maxpar
 1695 format(/,'  ****** STOP ******',/,i8,' parameters, arrays are ',
     2 'for ',i8)
      stop
c
c***** end of subroutine input3 *****
      end
