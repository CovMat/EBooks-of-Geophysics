cz- This file is taken from simul2000 algorithm by Clifford Thurber

      subroutine vel3(isp,x,y,z,v)
c  This routine is Cliff Thurber's
c  common block variables:
      common/weight/ wv(8),ip,jp,kp,kpg
      include 'RayFDD.inc'
c
c  use Prothero's intmap here
      call intmap(x,y,z,ip,jp,kp)
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
c	write(16,100)x,xl,ip,y,yl,jp,z,zl,kp
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
      return
c***** end of subroutine vel3 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine intmap(x,y,z,ip,jp,kp)
c  Modified by W. Prothero so a single call can get the indices
c  common block variables:
      include 'RayFDD.inc'
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

      subroutine ttmder(isp, dth1, stepl)
c
c  declaration statements:
      dimension in(8)     
c
c  common block variables:
      real dth1(4),  stepl
      integer isp
      common/weight/ wv(8),ip,jp,kp,kpg
      include 'RayFDD.inc'
c
c
c   This subroutine has been modified for Q inversion 13-March-1998.
c   Following Andreas Rietbrock's version.
c   Then (when iuseq=1), the second part of the model array
c   is used to store Vp*Qp (not Vp/Vs).
c   The solution is only for Q, not for Vp though.  
c   Thus the solution arrays that would otherwise be related
c   to the first part of the model array (Vp) are used.
c   Therefore in ttmder, vel3 is called with "1" to retrieve
c   the V*Q value; and vel3 is called with "0" to get the 
c   indices for the solution arrays.
c
c
c  resegment ray path
c  calculate travel time derivatives with respect to hypocentral
c  parameters - use coords of first two points on raypath
c  determine slowness at source
      xe=rp(1,1)
      ye=rp(2,1)
      ze=rp(3,1)              
      call vel3(isp,xe,ye,ze,v)       
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
c
         xp=rx-half*dxs
         yp=ry-half*dys
         zp=rz-half*dzs
c  loop over subsegments
         do 55 is=1,nseg
            xp=xp+dxs
            yp=yp+dys
            zp=zp+dzs
c
          call vel3(isp, xp, yp, zp, v)
          dt=ssl/v
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
c
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
c
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
c
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
cz                   inp=ini+(no-1)*npari
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
   46             continue
   47          continue
   48       continue
c
   55    continue
   50 continue
  88    continue
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
      include 'RayFDD.inc'
c
c     array size limits
c     ixkms=iykms=izkms=1500
c
c     write(6,400)
c 400 format(' subroutine bldmap')
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
   10 continue
c  Fill remainder of array with zeroes.
      do 12 i=ixmax,ixkms
         ixloc(i)=0
   12 continue
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
   15 continue
c
c  Fill rest of array with zeroes.
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
   20 continue
c
c  Fill remainder of array with zeroes.
      do 22 i=izmax,izkms
         izloc(i)=0
  22  continue
      return
 330   continue
      write(16,331)ixkms,iykms,izkms
 331  format(' ***** error in array size in common/locate/',/,
     *' maximum map dimensions (km)=',/,' x=',i5,' y=',i5,' z=',i5)
      write(16,332)ixmax,iymax,izmax
  332 format(' Actual map size (km): ',/,' x=',i5,' y=',i5,' z=',i5)
      stop
c***** end of subroutine bldmap *****
      end



      subroutine input_vel
c  this routine reads in the initial velocity model in the
c  form of velocity specified on a uniform but not
c  necessarily evenly spaced grid of points
c  (reads from file03 )
c
c  common block variables:
      include 'RayFDD.inc'
c
c  declaration statements:
      integer ixf(maxpar),iyf(maxpar),izf(maxpar)
c
c  start cht 1998

c  end cht 1998
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
c  for this version the gridpoints can be unevenly spaced
c  the origin of the coordinate system is at (x,y,z)=(0,0,0)
c  which will not in general correspond to the point
c  xn(1),yn(1),zn(1).
c  xn,yn,zn should be factors of bld (ie: a.0 for bld=1.0 or a.b for bld=0.1)
c
c input the number of gridpoints in x, y and z directions
c  and bld factor (1.0 or 0.1 km) used to set up velocity interpolation grid
      read(3,3002) bld,nx,ny,nz
 3002 format(f4.1,3i3)
      if((bld.ne.1.0).and.(bld.ne.0.1)) then
        write(16,1625) bld
 1625   format(/, '******** STOP *********, bld must be 1.0 or 0.1,
     2   not ',f6.2)
      endif
        atemp=iuses*(nx-2)*(ny-2)*(nz-2)
        if(atemp.le.maxpar)goto 40
        write(16,42)
 42     format('0Too many nodes for program array sizes.')
        stop
 40     continue
c
c  start cht 1998
      do 123 k=1,atemp
      imerge(k)=0
      jequal(k)=0
 123  continue
c
c  end cht 1998
c
c  input the x grid, y grid, and z grid
cfh read in free format (makes life easier...)
c        read(3,3004) (xn(i),i=1,nx)
c        read(3,3004) (yn(i),i=1,ny)
c        read(3,3004) (zn(i),i=1,nz)
        read(3,*) (xn(i),i=1,nx)
        read(3,*) (yn(i),i=1,ny)
        read(3,*) (zn(i),i=1,nz)
 3003 format(3i3)
 3004 format(20f6.1)
c
      write(16,3005) bld,nx,ny,nz
 3005 format(//,' velocity grid size:',/,
     * 'bld =',f4.1,5x,' nx =',i3,5x,'ny =',i3,5x,'nz =',i3)
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
      nx2=nx-2             ! number non-edge nodes in row
      nxy2=nx2*(ny-2)      ! number non-edge nodes in layer
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
   65 write(16,3101)
c     do 38 kv=1,iuses
         kv=1
         do 37 k=1,nz
            k2=k + (kv-1)*nz
            write(16,3015) k,vtype(kv),zn(k)
            do 36 j=1,ny
cfh               read(3,3011) (vel(i,j,k2),i=1,nx)
               read(3,*) (vel(i,j,k2),i=1,nx)
               write(16,3013) (vel(i,j,k2),i=1,nx)
   36       continue
   37    continue
c  38 continue
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
   99     continue
  100  continue
c  compute Vs from Vp and Vp/Vs or compute 1/tstar
        kv=2
       if(iuseq .eq. 0) then
           do 120 k=1,nz
              ks=k+nz
              write(16,3015) k,vtype(kv),zn(k)  
              do 115 j=1,ny
                 do 110 i=1,nx
                    vel(i,j,ks)=vel(i,j,k)/vpvs(i,j,k)
  110           continue
                 write(16,3013) (vel(i,j,ks),i=1,nx)
  115        continue
  120     continue
       else
         do 140 k=1,nz
            ks=k+nz
            write(16,3018) k,zn(k)
            do 135 j=1,ny
               do 130 i=1,nx
                  vel(i,j,ks)=vel(i,j,k)*qval(i,j,k)
  130           continue
                 write(16,3014) (vel(i,j,ks),i=1,nx)
  135        continue
  140     continue
       endif
      endif
c
 3013 format(20f6.2)
 3014 format(20f7.1)
 3015 format(/,' layer',i3,5x,a1,' velocity',10x,'z =',f7.1)
 3016 format(/,' layer',i3,5x,'Vp/Vs',10x,'z =',f7.1)
 3017 format(/,' layer',i3,5x,'Q',10x,'z =',f7.1)
 3018 format(/,' layer',i3,5x,'Q * Vp',10x,'z =',f7.1)
 3011 format(20f5.2)
 3101 format(//,' velocity values on three-dimensional grid')
c  Number of medium parameters to invert for
      npar=nodes2*iuses
      nparv=npar
cz      if(invdel.ne.0)npar=(npar+nsts*iuses)
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
cfhdmep           if(izf(i).gt.nz) iizf=izf(i)-4
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
 1686   format(' *-*-* ERROR velocity input.  This node has',
     2  ' already been ',/,' *-*-* assigned cnode= ',a1)
             ierror=1
             write(16,1681) ixf(i),iyf(i),izf(i),k
 1681        format('input3 fixed node',i5,' ixf,iyf,izf:',3i3,
     2       ' node number:',i5)
           endif
           cnode(k)='F'
   70   continue
c
       write(16,1611)
 1610 format(/,' velocity FIXED at the following nodes(1):')
 1611 format(/,' VELOCITY INVERSION GRID    0=free, F=Fixed',/,
     2  '   M=Master, C=Constant Pert. Link, ',
     3  'L=Linear Pert. Link')
  311 do 495 kv=1,iuses
         nz1=nz-1
         ny1=ny-1
         do 320 k=2,nz1
          if(iuseq.eq.0) then
            if(kv.eq.1) write(16,1009) k,vtype(kv),zn(k)
 1009       format(/,' layer',i3,5x,a1,'-velocity nodes',
     2         10x,'z =',f7.1)
            if(kv.eq.2) write(16,3016) k,zn(k)
          else
            write(16,3017) k,zn(k)
          endif
            kk=k+(kv-1)*nz2

            do 310 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1006) (cnode(i),i=n1,n2)
c               write(16,1005) (nfix(i),i=n1,n2)
  310       continue
  320    continue
c 1005 format('    ',18i6)
 1006 format('  ',40(2x,a1))
  495 continue
  496 continue
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
         if(km.ge.nz) km=km+2      ! if s velocity node
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
  888 continue
c  add to index if not linked or fixed
c
         in=in+1
         ndexfx(n)=in
         mdexfx(in)=n
c         write(16,1888) in,n
c 1888    format(' adding inversion node number ',
c     &   i5,' corresponding to gridpoint number ',i5)
c
   80 continue
      inf2=nparv-in
      if(inf2.eq.infl) goto 85
c
c  end cht 1998
c
        write(16,1615) infl,nparv,in,inf2
 1615   format(/,' **** number of fixed and linked nodes input,',i4,
     2  ', does not equal velocity nodes,',i4,', minus invert',
     3  ' nodes,',i4,'.  Continue with inf=',i5,' ****',/)
        infl=inf2
   85 continue
      nparvi=nparv-infl
      npari=nparvi
cz      if(invdel.ne.0) npari=nparvi+nstsi*iuses
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
   90 continue
      npari=in
   95 continue
      write(16,1620) npar,nparv,npari,nparvi
 1620 format(' INPUT3:npar,nparv,npari,nparvi',4i6)
c  Check to see whether medium parameters fit within array sizes
      if(npari.gt.mxpari) goto 990
      if(npar.gt.maxpar) goto 995
c  Stop if there was an error reading in fixed and linked nodes
      if(ierror.ne.0) then
        write(16,1683)
        write(6,1683)
 1683   format(/,'STOP SIMUL2000! , Error in velocity input file ')
        stop
      endif
c
c  Set up an array which is used to point to node indices, for any x,y,z
c      write(*,*)'nparsi=',nparsi
      call bldmap
c      write(*,*)'nparsi=',nparsi
c
      return
c
  980 continue
      write(6,1683)
      write(16,1698) nparv,maxpar
 1698 format(/,'  ****** STOP ******',/,i8,' velocity nodes, program',
     2 ' arrays only allow',i6)
      stop
  990 continue
      write(6,1683)
      write(16,1699) npari,mxpari
 1699 format(/,'  ****** STOP ******',/,i8,' parameters to invert for',
     2 ', program arrays only allow',i6)
      stop
  995 continue
      write(6,1683)
      write(16,1695) npar,maxpar
 1695 format(/,'  ****** STOP ******',/,i8,' parameters, arrays are ',
     2 'for ',i8)
      stop
c
c***** end of subroutine input3 *****
      end
