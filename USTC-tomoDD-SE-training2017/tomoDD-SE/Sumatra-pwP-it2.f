      program reglbP1

* ----------------------------------------------------------------------------

c modified for iterative solution of bounce point, Feb 2008

* NOTE: scale is to scale perturbations (if scaling is not needed --> scale=1.0)
*       ALWAYS check scale = ... ???

c      parameter(nunkn=500000,scale=1.0)
      parameter(nunkn=500000)
      parameter(np1=6379,np15=449,npcom=100000,nlayer=16,nlayeri=19)
      implicit real*8 (a-h, o-z)
      parameter (msg=256)

* Note: msg is the max number of ray segments in a ray path !!!

* Global parameterization (5 deg by 5 deg with 16 layers):
      parameter (nx=72,ny=36,nz=16)
      parameter (dx=5.,dy=5.)

* Regional parameterization (0.5 deg by 0.5 deg with 19 layers):
c milo = min co-longitude; malo = max co-longitude !!!
c mila = min co-latitude ; mala = max co-latitude  !!!

* SUMATRA: 90E to 135E; 25N to 15S --> use colon and colat as input
      parameter (milo=90.,malo=135.,mila=65.,mala=105.)
      parameter (dxi=0.5,dyi=0.5)
* GO TO CHECK!!!

      dimension  w(3,msg+1),ws(3,msg+1)
      real*8 xor(nunkn)
      character*54 head1,head2*154
      character*8 phasej
 
      dimension r(npcom),v(npcom),d(npcom)
      dimension r1(npcom),v1(npcom),d1(npcom)
      dimension xlayer(0:nlayer),dz(nlayer),dzsum(nlayer)
      dimension xlayeri(0:nlayeri)
      common /value/ pi,con,rcon,deg,r,v,d,
     >r1,v1,d1,dz,dzsum

      common /bpstore/ xdep1, xlat1,xlon1,xdep2, xlat2,xlon2,
     & xdep3, xlat3,xlon3,xdep4, xlat4,xlon4

      common /model/ xlayer,xlayeri

      common /perturb/ xan(nunkn)

      pi=4.0*atan(1.0)
      deg=111.195
      con=pi/180.
      rcon=180./pi
      dpi = asin(1.)/ 90.
      r2d = 90./asin(1.)
      ieq=0
      rms=0.0e10

* The phase data file:
      open(2,file=
     >'TestpwP-GLB',
     >form='formatted')

* To write calculated time delays:
      open(3,file=
     >'residual1.dat',
     >form='formatted')

* To write relocation coefficients:
      open(11,file=
     >'reloc1.dat')
      
* Read the starting 1-D velocity model:
      open(10,file='/mnt/s24/sumatra/IlikUW/AK135/ak135.15.SKS')
      rewind(10)
 
* To save the calculated ray-segment lengths:
      open(14,file=
     >'raylength1.dat')

* Read the layer division for global parameterization:
      open(19,file=
     >'/mnt/s24/sumatra/IlikUW/P_Sumatra/RAY_3D/layer-16.dat')

* Read the layer division for regional parameterization:
      open(15,file=
     >'/mnt/s24/sumatra/IlikUW/P_Sumatra/RAY_3D/layer-19.dat')

* To provide anomalies from 3-D models: 
* In the first iteration, comment out the following 
* open 20 & 21 when we start using a 1-D velocity model)

* Local 3-D model (velocity perturbations in %):
c     open(20,file=
c    >'/home/sriwid/P_VRANCEA/MODEL/P_Vrancea-reg1x1.it1-w23')

* Global 3-D model (velocity perturbations in %):
c     open(21,file=
c    >'/home/sriwid/P_VRANCEA/MODEL/P_Vrancea-glb5x5.it1-w23')

* Save 'kounter':
      open(22,file=
     >'P-kount-it1')

      kount=0
      nmdat=0
      numP=0
      numpP=0
      numpwP=0

* NOTE AK135 contains: radius,depth,Vp,Vs,density
* For P-wave then v(i)=vp !!!
 
      do i=1,np15

      read(10,*)r(i),d(i),vp,vs,rho
      zzz=d(i)
      d(i)=float(nint(d(i)))
      r(i)=float(nint(r(i)))
      v(i)=vp

      enddo
 
***************************************************

* Read global layer division data:
      do i=0,nlayer
      read(19,*)xlayer(i)
      enddo
* Read regional layer division data:
      do i=0,nlayeri
      read(15,*)xlayeri(i)
      enddo
 
***************************************************

      nxi=(malo-milo)/dxi
      nyi=(mala-mila)/dyi
      nzi=nlayeri
      nmax1=nxi*nyi*nzi
      nmax2=nx*ny*nz+nmax1

* Read the target block number:
* comment out 20 and 21 for 1st iteration !!! 
* (as it starts from 1-D ak135 vel model)

* Regional model:
c     read(20,911)head1
c     read(20,912)head2
c     read(20,*)(xor(i),i=1,nmax1)
cc     xor(nmax1+1)=0.0

* Global model:
c     read(21,911)head1
c     read(21,912)head2
c     read(21,*)(xor(i),i=nmax1+1,nmax2)
cc     xor(nmax2+1)=0.0

* NOTE: For calculating tt in ak135, 'zero' the following xan(j):

      do j=1,nmax2+1

* Tracing at 3D model: (for next iterations to update 3-D model)
c        xan(j)=xor(j)*0.01*scale
* Tracing at 1D model: (for 1st itereration)
         xan(j)=0.0

      enddo

****************************************

10    read(2,150,end=200)
     >nblock,nst,res,iphj,jev,vlon,vlat,depth,rp,
     >slon,slat,backaz,azclst,
     >ihold,iflg2,dbot,ymp,dist,azsr,scor,
     >resisc,iphi,phasej,tdelta,ntel,prec,
     >obstt,prett,ecor,scor,elcor,
     >gblon,gblat,stadel,bdep,tbath,twater
 
****************************************

c  Add Corrections to obstt (JP Oct. 07)
      delay=obstt-ecor-elcor-tbath-twater
      nmdat=nmdat+1
 
* Residual-time computation for P waves:
        tt=0.0
        tta=0.0

c NOTE: All data filtering now done prior to ray-tracing  by 'assemble_ehb.f'
c [JP Jan 08]
         
*****************************************
        
        if(slon.GT.180.)then
           slon=slon-360.
        endif

c JP edit - zero initial values
      aas=0
      bbs=0
      hs=0
      aar=0
      bbr=0
      hr=0

c  CT edit - add pP - sept 2007
        if(phasej.eq.'P       '.or.phasej.eq.'Pn      ') then

        numP=numP+1

        hr=0.0e10 
        aar=slat
        bbr=slon
        hs=depth
        aas=vlat
        bbs=vlon
        dth2=0.0
        dth3=0.0
        dth4=0.0

        call pbr(aas, bbs, hs, aar, bbr, hr,w,ni)
        ieq=ieq+1

        call cderiv( w, ni, dth2, dth3, dth4)
        dth1=1.0
        write(11,120)nblock,dth1,dth2,dth3,dth4,jev,vlon,vlat,depth

        ichr=14
        ichd=3

c  CT edit - add pP - sept 2007
        iflag=0
        call rayl(ichr,ichd,w,ni,delay,kount,iflag)
        rms=rms+delay**2

c  CT edit - add pP - sept 2007
      else

c  JP edit - add kounter for pP and pwP - Jan. 08
         if (phasej.eq.'pP      ') then
            numpP=numpP+1
         else
            numpwP=numpwP+1
         endif

c  check that event is within study region [fix lat/colat (JP)]
      if ((90-vlat).lt.mila.or.(90-vlat).gt.mala.or.vlon.lt.milo.or
     &    .vlon.gt.malo) go to 10
c
      ibp=0

c  JP edit - save original values 
      vlato=vlat
      vlono=vlon
      slato=slat
      slono=slon
      deptho=depth

        ieq=ieq+1
 199  continue

c  first do event to bounce point
        hr=0.0e10 
        aar=gblat
        bbr=gblon
        hs=depth
        aas=vlat
        bbs=vlon
        dth2=0.0
        dth3=0.0
        dth4=0.0

        call pbr(aas, bbs, hs, aar, bbr, hr,w,ni)

c path into bounce point - store point 2 and check 1
      bpuplat=90-xlat2
      bpuplon=xlon2
      bpupdep=6371-xdep2
c
c verify that xlat4-xlon4-xdep4 equals original bounce point within tolerance
      if (abs(90-xlat1-gblat).gt.0.02.or.abs(xlon1-gblon).gt.0.02.or.
     & abs(6371-xdep1).gt.2) write(6,9876) 90-xlat1,xlon1,(6371-xdep1),
     & gblat,gblon,hr
 9876 format('  ** ERROR - p PATH END NOT AT BOUNCE POINT **',6f10.4)

        call cderiv( w, ni, dth2, dth3, dth4)
        dth1=1.0
c        write(11,120)nblock,dth1,dth2,dth3,dth4,jev,vlon,vlat,depth

c store these values and write out when pP path converges
      nblocks=nblock
      dth1s=dth1
      dth2s=dth2
      dth3s=dth3
      dth4s=dth4
      jevs=jev
      vlons=vlon
      vlat=vlats
      depth=depths
c
        ichr=14
        ichd=3

c  CT edit - add pP - sept 2007
c store these values and call rayl later after convergence
      ichrs=ichr
      ichds=ichd
      do 98 k=1,ni+1
      do 97 kk=1,3
      ws(kk,k)=w(kk,k)
   97 continue
   98 continue
      nis=ni
      delays=delay
      kounts=kount
      iflags=1
c
c       iflag=1        
c       call rayl(ichr,ichd,w,ni,delay,kount,iflag)
c        rms=rms+delay**2

c  now do bounce point to station
        hr=0.0 
        aar=slat
        bbr=slon
        hs=0.0
        aas=gblat
        bbs=gblon

        call pbr(aas, bbs, hs, aar, bbr, hr,w,ni)
c
c path away from bounce point - store point 3 and check 4
      bpdnlat=90-xlat3
      bpdnlon=xlon3
      bpdndep=6371-xdep3
c
c verify that xlat1-xlon1-xdep1 equals original bounce point within tolerance
      if (abs(90-xlat4-gblat).gt.0.02.or.abs(xlon4-gblon).gt.0.02.or.
     & abs(6371-xdep4).gt.2)write(6,9876) 90-xlat4,xlon4,(6371-xdep4),
     & gblat,gblon,hr
c
c        call cderiv( w, ni, dth2, dth3, dth4)
c        dth1=1.0
c        write(11,120)nblock,dth1,dth2,dth3,dth4,jev,vlon,vlat,depth
c
c check for too many bp iterations
      if (ibp.gt.10) go to 222
c
        ichr=14
        ichd=3
c
c scale bpdn lat and lon to same depth as up
      bpdnlats=gblat+(bpupdep/bpdndep)*(bpdnlat-gblat)
      bpdnlons=gblon+(bpupdep/bpdndep)*(bpdnlon-gblon)
c
c now check differences against threshold
      gblatpre=0.5*(bpuplat+bpdnlats)
      gblonpre=0.5*(bpuplon+bpdnlons)
c
      if (abs(gblat-gblatpre).gt.0.005.or.abs(gblon-gblonpre).gt.0.005)
     & then
      gblat=gblatpre
      gblon=gblonpre

c JP edit - reset to original values
      vlat=vlato
      vlon=vlono
      slat=slato
      slon=slono
      depth=deptho

      ibp=ibp+1

      go to 199
c
      endif
c
  222 continue
c
c output converged result
c
      write(11,120)nblocks,dth1s,dth2s,dth3s,dth4s,jevs,
     & vlons,vlats,depths
c
       call rayl(ichrs,ichds,ws,nis,delays,kounts,iflags)
c
c  CT edit - add pP - sept 2007
        iflag=2
        call rayl(ichr,ichd,w,ni,delay,kount,iflag)
        rms=rms+delay**2

      endif
      goto 10

200   continue

      write(22,*)kount,nmdat,numP,numpP,numpwP

* ---------------------------------------------------------------

120   format(i7,4f10.2,i6,f8.3,f8.3,f6.1)

150   format(2i6,f7.2,i4,i6,f8.3,f8.3,
     >f6.1,
     >f8.4,
     >f8.3,f8.3,2f9.3,
     >2i2,f7.2,f7.2,f7.2,f9.3,f7.2,f7.2,i4,a8,f8.3,i5,i3,
     >2f10.2,3f7.2,
     >3f8.3,f7.3,2f7.2)

555   format(i7,1x,f12.2)
666   format(f7.2)

911   format(a54)
912   format(a154)

      rms=sqrt(rms/(ieq-1))

      write(*,*)'Total used rays: ',ieq
      write(*,*)'Number of P, pP, pwP phases used: ',numP,numpP,numpwP
      write(*,*)'Total RMS: ',rms
 
      stop 'No worries mate ...'
      end
 
* ---------------------------------------------------------------

      subroutine rayl(ichr,ichd,w,ni,delay,kount,iflag)
      implicit real*8 (a-h, o-z)
      real*8 la,lo
      parameter (msg = 256)
      dimension  w(3,msg+1)

* dinter - is the step for interpolation
* finding the boundary of cell

      dinter=1.000
      rlength=0.0e10

c only zero rdelay for P or final segment for pP (JP)
      if(iflag.eq.0.or.iflag.eq.1)then
        rdelay =0.0e10
      endif

      r2d = 90./asin(1.)
      dpi = asin(1.)/ 90.

* ibl is a flag that is written in order to mark
* rays penetraiting lower 2900 km.
      ibl=0

* Begin here for each point:

      i=2
      do while(i.lt.ni+2)

* Find location of two points:

      call cellvel(w(2,i-1),w(3,i-1),w(1,i-1),Vel1,nbl1)
      call cellvel(w(2,i),w(3,i),w(1,i),Vel2,nbl2)

* is it in the same block? 
* if Yes - find the parameters
* and begin again for next two points (see go to 99)

         if(nbl1.EQ.nbl2)then

c  bug fix, CT, Oct 2007 - set kflag
            kflag=1
         
* Find length between points:

            call dist(w(2,i-1),w(3,i-1),w(1,i-1),
     &                w(2,i),w(3,i),w(1,i),dst)

            rlength=rlength+dst
            rdelay =rdelay+dst/((Vel1+Vel2)/2.) 
            go to 99

         endif
          
* if they are not in the same block then
* find the last point in the same block

         if(nbl1.NE.nbl2)then

c  bug fix, CT, Oct 2007 - unset kflag
            kflag=0
         
         Vel2=Vel1

* find last point within the block
* firstly find the distance between points

            call dist(w(2,i-1),w(3,i-1),w(1,i-1),
     &                w(2,i),w(3,i),w(1,i),dst)

* find number of steps for fixed length dinter:

            nms=int(dst/dinter)

* find dx, dy and dz for each step:

                  ar=w(2,i-1)*dpi
                  as=w(2,i)*dpi
                  br=w(3,i-1)*dpi
                  bs=w(3,i)*dpi
                  rr=w(1,i-1)
                  rs=w(1,i)
                 x1 = rr*sin(ar)*cos(br)
                 y1 = rr*sin(ar)*sin(br)
                 z1 = rr*cos(ar)
                 x2 = rs*sin(as)*cos(bs)
                 y2 = rs*sin(as)*sin(bs)
                 z2 = rs*cos(as)
                 dx = (x2-x1)/nms
                 dy = (y2-y1)/nms
                 dz = (z2-z1)/nms

* begin looking for the last point:

          do 10 j=1,nms

* gradually increase distance toward the second point:
          x = x1 + dx*j
          y = y1 + dy*j
          z = z1 + dz*j
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

c         write(*,*)'tut?',i

* where is it ?
            call cellvel(la,lo,r,Velc,nbl)      

* if it is in the next cell begin to find parameters
* for the last distance, redefine first point and
* begin for new two points ( see go to 66)
* if they still  are in the same block
* increase distance more and performe approximation again

            if(nbl.NE.nbl1)then 

               call dist(w(2,i-1),w(3,i-1),w(1,i-1),
     &                la,lo,r,dst)

                   rlength=rlength+dst
                   rdelay =rdelay+dst/((Vel1+Vel2)/2.)

* Write results into the file if ray is within grid net:

                   if(nbl1.GT.0)then
                     write(ichr,*)nbl1,rlength
                     kount=kount+1
                   else
                     ibl=1
                   endif

                   rlength=0.0e10
                   w(2,i-1)=la
                   w(3,i-1)=lo
                   w(1,i-1)=r
                              go to 66
            endif

          Vel2=Velc

10    continue 


         endif
99    continue
      i=i+1
66    continue

      enddo

c  bug fix, CT, Oct 2007 - check kflag
                   if(kflag.eq.1)then
                     write(ichr,*)nbl1,rlength
                     kount=kount+1
                     rlength=0.0e10
c                     write(ichr,*)'0 -1'
c                     kount=kount+1
c                     iflag=0
                   endif

* Ray tracing is finished:

c   CT edit - add pP - sept 2007
      if (iflag.eq.0.or.iflag.eq.2) then
        write(ichr,*)'0 -1'
        write(ichd,*)rdelay,delay,(delay-rdelay),' ',ibl
        delay=delay-rdelay
        kount=kount+1
      endif
      
      return
      end

* ----------------------------------------------------------------------

      subroutine dist(rla1,rlo1,r1,rla2,rlo2,r2,dst)
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

* ----------------------------------------------------------------------

      subroutine km2deg(ala,alo,adp,dx,dy,bla,blo,bdp)
      implicit real*8 (a-h, o-z)

* This subroutine calculates the position of new point
* in polar coordinates based on the coordinates
* of main point in radians ( la is colatitude) and dx and dy in kilometers

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

* ----------------------------------------------------------            

c  3D pseudo-bending for the continuous, spherical earth
c  based on Kazuki Koketsu (ERI, Univ. Tokyo, 1998) algorithm
c
c  ray is tracing from receiver to source
c  input coordinates are in degrees and km
c  latitude from -90 to 90 and longitude from -180 to 180
c  depth is radius from the earth's center

      subroutine pbr(aas, bbs, hs, aar, bbr, hr,w,ni)    
      implicit real*8 (a-h, o-z)                                        
      parameter (msg = 256)
      dimension  r(msg+1), a(msg+1), b(msg+1)
      dimension  w(3,msg+1)
      integer ni,i
      external   vel,rtim,rlen
      real*8 RNULL
      common /coord/ shiftlo
      common /bpstore/ xdep1, xlat1,xlon1,xdep2, xlat2,xlon2,
     & xdep3, xlat3,xlon3,xdep4, xlat4,xlon4
      data      ro / 6371.00 /
      data    RNULL /0.0e10/

c  JP edit - store original values   
      aaso=aas
      bbso=bbs
      hso=hs
      aaro=aar
      bbro=bbr
      hro=hr
                                                        
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
      n2     = 256
      nloop  = 200
      flim   = 1.e-2/100.
      mins   = 20.


      dpi = asin(1.)/ 90.
      r2d = 90./asin(1.)

* Check coordinates
c      write(*,*)aas,bbs,aar,bbr     
      if(aas.LT.-90.OR.aas.GT.90.)then
        write(*,*)'Latitude of source is out of range'
        stop
      endif
      if(aar.LT.-90.OR.aar.GT.90.)then
        write(*,*)'Latitude of station is out of range'
        stop
      endif
      if(bbs.LT.-180.OR.bbs.GT.180.)then
        write(*,*)'Longitude of source is out of range'
        stop
      endif
      if(bbr.LT.-180.OR.bbr.GT.180.)then
        write(*,*)'Longitude of station is out of range'
        stop
      endif

* Rotate coordinates in order to have
* longitude and latitude range from 0 to 180. 
* This program does not work with angles
* greater than 180.       

* Pass from latitude to colatitude:

      as = (90.00-aas) * dpi
      ar = (90.00-aar) * dpi

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

c        write(*,*)'bbr,bbs,shift'
c        write(*,*)'aar,aas'
c        write(*,*)aar,aas
c        write(*,*)'bbr,bbs,shiftlo'
c        write(*,*)bbr,bbs,shiftlo

         bs = bbs * dpi
         br = bbr * dpi  

      ad = (as + ar) / 2.                                               
      rs = ro - hs
      rr = ro + hr                                                      
 
c *** initial straight ray ***                                           
c     ni : number of ray segments

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
      to = rtim(ni+1,r,a,b)
      tp = to
c      write(6,*) to
      do 112  i=1,ni+1
          w(1,i) = r(i)
          w(2,i) = a(i)
          w(3,i) = b(i)
 112  continue

c *** number of points loop ***
      loops = 0
      do while(ni .le. n2)
                                                                        
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
 
                  dn = dx**2 + dy**2 + dz**2                        
                  ddn = sqrt(dn)                                    
                  dr = (r3-r1) / ddn                                
                  da = (a3-a1) / ddn                                
                  db = (b3-b1) / ddn

* Begin to find the gradients and velocities:
* first find the length of segment

                  dseg=sqrt((dx/2)**2+(dy/2)**2+(dz/2)**2)

                  ddseg=dseg/2.
c                 write(*,*)'dseg ',dseg 

* Now ddseg will be a distance to find dV
* along the coordinates 

* Determine velocity at 3 points:

c                 write(*,*)'before',ni
c                 write(*,*)
c           write(*,*)'Pered opredeleniem glubini ------------------'
c           do j=1,ni+1
c           write(*,*)j,ni,r(j),a(j)*r2d,b(j)*r2d
c           enddo
c                 write(*,*)'grohnulas pri',ni,' segmentov'
c                 write(*,*)'l,k ',l,k
c                 write(*,*)r1,a1*r2d,b1*r2d


                  v1 = vel(r1,a1,b1,1)                            

c                 write(*,*)'V1= ',v1
c                 write(*,*)r2,a2*r2d,b2*r2d

                  v2 = vel(r2,a2,b2,2)                            

c                 write(*,*)'V2= ',v2
c                 write(*,*)r3,a3*r2d,b3*r2d

                  v3 = vel(r3,a3,b3,3)                            

c                 write(*,*)'V3= ',v3

* Begin to determine coordinates
* of pints surroundibg point a2,b2,r2
* at the distance ddseg

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

* The following if-endif is just for P & S, thus comment out for SKS & PKP !!!
* This gives the lowermost mantle Vp in the outer core
                   if(dwz.le.3479.50)then
                      dwz=3479.500
                      upz=dwz+dseg
                   endif
* ----------------------------------------------------------------------------

* Find dV along depth:

                  vr1 = vel(upz,a2,b2,4)
                  vr2 = vel(dwz,a2,b2,5)
                  vr=(vr1-vr2)/dseg

c                 write(*,*)'vr= ',vr

* Find dV along longitude:

                  call  km2deg(a2,b2,r2,ddseg,RNULL,adV,bdV,rdV) 
                  vb2 = vel(rdV,adV,bdV,6)
                  call  km2deg(a2,b2,r2,-1.*ddseg,RNULL,adV,bdV,rdV)
                  vb1 = vel(rdV,adV,bdV,7)
                  vb=-1.*(vb1-vb2)/dseg

c                 write(*,*)'vb= ',vb

* Find dV along latitude:

                  call  km2deg(a2,b2,r2,RNULL,ddseg,adV,bdV,rdV)
                  va2 = vel(rdV,adV,bdV,8)
                  call  km2deg(a2,b2,r2,RNULL,-1.*ddseg,adV,bdV,rdV)
                  va1 = vel(rdV,adV,bdV,9)
                  va=-1.*(va1-va2)/dseg

c                 write(*,*)'va= ',va
c                 write(*,*)'**************'
               
c               write(*,*)'after',ni


c           do j=1,ni+1
c           write(*,*)j,ni,r(j),a(j)*r2d,b(j)*r2d
c           enddo
      
* Spherical velocity gradient:

c                 va = va / r2
c                 vb = vb / r2 / sina

c                 (tangential vector) = (slowness vector) / s
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

c  "Tut esli rcur < 0.0 proishodit hernia
c   poetomu postavlen abs. Ne yasno mozhno li eto delat
c   ili net no rabotaet. Obichno oshibka poyavliaetsia
c   ochen redko v nekotorih tochkah
c   v etom sluchae abs prosto ne daet oshibki y posledniaya iteraciya
c   uzhe ne imeet rcur negativnim y podgoniaet normalno reshenie
c   (mozhet bit)"

c                     write(*,*)'rcur ', rcur
                      if(rcur.LE.0.0)then
                      write(*,*)'Negative'
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
                      rp  = r2 + rdr
                      ap  = a2 + rda/r2
                      bp  = b2 + rdb/(r2*sina)
c                     write(*,*)rp,ap,bp
                      r(k) = (rp-r(k))*xfac + r(k)
                      a(k) = (ap-a(k))*xfac + a(k)
                      b(k) = (bp-b(k))*xfac + b(k)
                  endif
c             write(*,*)'after vector',ni
c             write(*,*)r(k),a(k)*r2d,b(k)*r2d
c             write(*,*)'++++++++++++++++++++++++++++++++++++++'
 2250         continue                                              
c                                                                       
c             write(*,*)'pered ciclom', ni
              idstn=ni
              do 212  j=1,ni+1
                  w(1,j) = r(j)
                  w(2,j) = a(j)
                  w(3,j) = b(j)
c             write(*,*)j,ni,w(1,j),w(2,j)*r2d,w(3,j)*r2d
 212          continue
              ni=idstn
c              here trace each iteration
c             write(*,*)'pered rtim',ni
              tk = rtim(ni+1,r,a,b)
           
c             write(*,*)to,tk,ni
c             write(6,*) ni, l, (tk-atim)/atim*100.
c               tk > to $b$n;~$bdd;_$5$;$k$?$a(b abs $b$oiu$1$j$$!%(b
              if(abs(to-tk) .le. to*flim)  go to 310                        
              to = tk                                                       
 3000     continue                                                          
 310      continue   
              to=tk

* skip increasing of segment number if minimum length
* of segment is exceed or maximum number of segments
* was reached

           if(dseg.lt.mins.or.ni.ge.n2) go to 66666
 
* Double the number of points:
c           write(*,*)'double number'

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
          tk = rtim(ni+1,r,a,b)
c           tk > tp $b$n;~$bdd;_$5$;$k$?$a(b abs $b$oiu$1$j$$!%(b
c here i change tp and put to
          if(abs(to-tk) .le. to*flim)  go to 99999
          to = tk 
      enddo
c                                                                       
99999 continue
c     write(*,*)'poslednii',ni
      idstn=ni
      do i=1,ni+1
          w(1,i) = r(i)
          w(2,i) = a(i)
          w(3,i) = b(i)
      enddo
      ni=idstn
66666 continue

* Return coordinates to the origin:

      idstn=ni
      do k=1,ni+1
          w(2,k) = w(2,k)*r2d
          w(3,k) = w(3,k)*r2d+shiftlo
       if(w(3,k).lt.0.)w(3,k)=360.+w(3,k)
c     write(*,*)k,ni,w(1,k),w(2,k),w(3,k)
      enddo 
c
      xdep1=w(1,1)
      xlat1=w(2,1)
      xlon1=w(3,1)
      xdep2=w(1,2)
      xlat2=w(2,2)
      xlon2=w(3,2)
      xdep3=w(1,ni)
      xlat3=w(2,ni)
      xlon3=w(3,ni)
      xdep4=w(1,ni+1)
      xlat4=w(2,ni+1)
      xlon4=w(3,ni+1)
c
      ni=idstn
c     do i=1,ni+1
c     write(*,*)w(1,i),w(2,i),w(3,i)
c     enddo
c     write(*,*)'Ray is done with ',ni,' segments'
c     write(*,*)aas,bbs,hs,aar,bbr,hr

c  JP edit - reset to original values    
      aas=aaso
      bbs=bbso
      hs=hso
      aar=aaro
      bbr=bbro
      hr=hro

      return
      end                                                               
 
* --------------------------------------------------------------

      function rlen(n, r, a, b)                                         
      implicit real*8 (a-h, o-z)                                        
      parameter (msg = 256)
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
 
* --------------------------------------------------------------------

      function rtim(m, r, a, b)                      
      implicit real*8 (a-h, o-z)                                        
      parameter (msg = 256)
      dimension  r(msg+1), a(msg+1), b(msg+1)
      integer m
      if(m.GT.(msg+1))write(*,*)'*'
      rtim = 0.
      rv1 = 1./vel(r(1),a(1),b(1),-1)
      do 100  j=1,m-1
          x1 = r(j)*sin(a(j))*cos(b(j))
          y1 = r(j)*sin(a(j))*sin(b(j))
          z1 = r(j)*cos(a(j))                                               
          x2 = r(j+1)*sin(a(j+1))*cos(b(j+1))
          y2 = r(j+1)*sin(a(j+1))*sin(b(j+1))
          z2 = r(j+1)*cos(a(j+1))
          dl = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
          rv2 = 1./vel(r(j+1),a(j+1),b(j+1),-2)
          sm = (rv1 + rv2) / 2.         
          rtim = rtim + sqrt(dl)*sm
          rv1 = rv2
  100 continue                                                          
      return                                                            
      end                                                               
 
* --------------------------------------------------------------------

      function vel(r, pa, ra, k)                              
      implicit real*8 (a-h, o-z)                                    
      common /coord/ shiftlo
      r2d = 90./asin(1.)                                                        

* Convert to degree and rotate coordinates :

      clat=pa*r2d
      clon=ra*r2d+shiftlo
      if(clon.lt.0.)clon=360.+clon

c     write(*,*)k
c     if(k.EQ.-1)write(*,*)clat,clon,r
      call cellvel(clat,clon,r,V,nbl)
c     write(*,*)'exit vel'
      vel=V

      return
      end

* --------------------------------------------------------------------

      subroutine cellvel(la,lo,era,Vl,nbl)
      implicit real*8 (a-h, o-z)

* This subroutine determines the velocity in the
* block using co-Lat, Lon and R.
* Lat and Lon in degrees
* era - is the radius
* dp is the depth

      parameter (nx=72,ny=36,nz=16)
      parameter (dx=5.,dy=5.,nunkn=500000)
      parameter(np1=6379,np15=449,npcom=100000,nlayer=16)

c CHECK!!!
c Note, the coordinates are colongitude and colatitude !!!

c INA: 90 E to 155 E; 15N to 20S
*     parameter (milo=90.,malo=155.,mila=75.,mala=110.)
*     parameter (dxi=0.5,dyi=0.5,nlayeri=19)
                                                         
c VRANCEA: 0E to 65E; 60N to 25N
*     parameter (milo=0.,malo=65.,mila=30.,mala=65.)
*     parameter (dxi=0.5,dyi=0.5,nlayeri=19)

c INDONESIA : 90E to 135E; 25N to 15S --> use colon and colat as input
      parameter (milo=90.,malo=135.,mila=65.,mala=105.)
      parameter (dxi=0.5,dyi=0.5,nlayeri=19)
      real*8 la,lo,era,dp
      integer  num
      dimension r(npcom),v(npcom),d(npcom)
      dimension r1(npcom),v1(npcom),d1(npcom)
      dimension xlayer(0:nlayer),dz(nlayer),dzsum(nlayer)
      dimension xlayeri(0:nlayeri)

      common /value/ pi,con,rcon,deg,r,v,d,
     >r1,v1,d1,dz,dzsum

      common/model/xlayer,xlayeri
      common /perturb/ xan(nunkn)

* First determine the number of layers:

      Re=6371.00
      dp=Re-era
     
      nlr=0
      nla=0
      nlo=0

* Additional parameters for the study region:

      nxi=(malo-milo)/dxi
      nyi=(mala-mila)/dyi
      nzi=nlayeri

      if(dp.LT.0.00)then
c       dp=Re+dp
        dp=0.0
      endif

      do while(lo.GE.360.)
         lo=lo-360.
      enddo

      do while(lo.LT.0.0e10)
         lo=lo+360.00
      enddo

      if(la.LT.0.000)then
         la=abs(la)
      endif

      if(la.GT.180.)then
         la=360.-la
      endif

      if(la.EQ.180.)then
         nla=nla-0.0001000
      endif

      do i=2,np15
 
       if(dp.GE.d(i-1).AND.dp.LT.d(i))then
           ilr=i-1
           go to 5
       endif    

      enddo 
5     continue

****  IF   *******************************************************

* Where is the point - inside or outside the study region?

      if(lo.GE.milo.AND.lo.LE.malo.AND.la.GE.mila
     &   .AND.la.LE.mala)then

* inside of the study region:

      if(dp.GT.xlayeri(nzi))then
c     write(*,*)'Velocity determination
c    & is deeper than regional nodes go to the global scale'
      go to  24
      endif

* number of the layers in the study region:

      do i=1,nzi

       if(dp.GE.xlayeri(i-1).AND.dp.LT.xlayeri(i))then
           nlr=i
           go to 21
       endif

      enddo
21     continue

* Then column (lon) number:

      do i=1,nxi

       if(lo.GE.(milo+real(dxi*real(i-1))).AND.
     &    lo.LT.(milo+real(dxi*real(i))))then
           nlo=i
           go to 22
       endif

      enddo
22    continue

c     and then Row (lat) number

      do i=1,nyi

       if(la.GE.(mila+dyi*real(i-1)).AND.
     &    la.LT.(mila+dyi*real(i)))then
           nla=i
           go to 23
       endif

      enddo
23     continue



      num=(nla-1)*nxi+nlo+(nlr-1)*nxi*nyi

** ELSE *************************************** ELSE

      else

* outside the region of study:

24     continue
 
      if(dp.GT.xlayer(nz))then

c     write(*,*)'Velocity determination
c    &                               is deeper than nodes'

      go to 4
      endif

      do i=1,nz

       if(dp.GE.xlayer(i-1).AND.dp.LT.xlayer(i))then
           nlr=i
           go to 1
       endif

      enddo 
1     continue

* Fixed boundaries:

      if(dp.ge.xlayer(nz))nlr=nz
      if(dp.le.0)nlr=1

* Then column (lon) number:

      do i=1,nx
       
       if(lo.GE.real(dx*real(i-1)).AND.
     &    lo.LT.real(dx*real(i)))then 
           nlo=i 
           go to 2 
       endif 

      enddo
2     continue

* and then row (lat) number:

      do i=1,ny

       if(la.GE.dy*real(i-1).AND.la.LT.dy*real(i))then
           nla=i 
           go to 3 
       endif

      enddo
3     continue

      if(nlr.LE.0.OR.nlr.GT.nz)then
         write(*,*)dp,'NLR: ',nlr
         stop 'error of depth'
      endif
      if(nlo.LE.0.OR.nlo.GT.nx)then
         write(*,*)lo,'NLO: ',nlo
         stop 'error of LON'
      endif
      if(nla.LE.0.OR.nla.GT.ny)then
         write(*,*)la,'NLA: ',nla
         stop 'error of LAT'
      endif

* The number of cell in global parameterization will be

      num=(nla-1)*nx+nlo+(nlr-1)*nx*ny+nxi*nyi*nzi

      endif

******************* ENDIF ******************

* and velocity in this block
c      write(*,*)num,nla,nlo,nlr
      Vl=v(ilr)*(1.+xan(num))
      nbl=num
      return
4     continue
      Vl=v(ilr)
      nbl=-1
      return
      end

* ---------------------------------------------------------------------

      subroutine cderiv(w,ni,dth2,dth3,dth4)
      implicit real*8 (a-h, o-z)
      parameter (msg = 256)
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
      call cellvel(pe*r2d,re*r2d,he,ve,nbl)
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

* ---------------------------------------------------------------------

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

