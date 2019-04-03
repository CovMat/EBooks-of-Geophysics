c--- All rights reserved by (c) Haijiang Zhang and Clifford Thurber, 2005-2015.
c--- This code cannot be used for commercial purposes without permission from authors.

c tomoDD-SE: Double-difference tomography method for regional scale tomography
c This version uses the spherical pseudo-bending method to calculate travel times and ray paths between events and stations
c
c
c References:
c For a detailed description of the algorithm see:
c    Waldhauser, F. and W.L. Ellsworth, A double-difference earthquake
c    location algorithm: Method and application to the northern Hayward
c    fault, Bull. Seismol. Soc. Am., 90, 1353-1368, 2000.
c
c    Zhang, H., and C. H. Thurber (2003), Double-Difference Tomography: The Method and Its Application to the Hayward Fault, California, c    Bull. Seism. Soc. Am, 93, 1875-1889. 
c
c    Zhang, H., and C. Thurber (2006), Development and applications of double-difference tomography, Pure and Applied Geophys., 163, c    373-403, doi:10.1007/s00024-005-0021-y. 
c
c The code is continuously being updated and improved, so feel
c free to contact Haijiang Zhang at zhang11@ustc.edu.cn

c------------------------------------------------------------------------------------------------------------------
c       SOFTWARE LICENSE AGREEMENT

c     “DOUBLE-DIFFERENCE SEISMIC TOMOGRAPHY”

c      FOR ACADEMIC RESEARCH USE ONLY 

c YOU MUST ACCEPT THE FOLLOWING TERMS TO RECEIVE PERMISSION TO USE THIS SOFTWARE.
c TO ACCEPT THE FOLLOWING TERMS, SELECT "YES."  IF YOU DO NOT ACCEPT THESE TERMS, SELECT "NO," 
c AND DO NOT USE THE SOFTWARE.

c This Software License Agreement (“Agreement”) is by and between you, an individual (“you”) 
c and the Wisconsin Alumni Research Foundation ("WARF"), a nonstock, nonprofit Wisconsin corporation.
c Under this Agreement, WARF is granting you a limited non-commercial, non-transferable right to 
c install, execute and use the tomoDD Software (the “Software”) under United States Copyright 
c No. 6,813,639 (the “Copyright”), solely as necessary to install, execute and use the Software 
c in accordance with this Agreement, in the fields of academic research or non-profit scholarship 
c or both (the “Licensed Field”).

c BY SELECTING “YES,” YOU HEREBY REPRESENT THAT YOU ARE AN INDIVIDUAL EMPLOYED BY AN ACADEMIC OR 
c OTHER NON-PROFIT RESEARCH INSTITUTE (“Institution”) AND AGREE TO BE BOUND BY THE TERMS OF THIS 
c AGREEMENT, INCLUDING THE WARRANTY DISCLAIMERS, LIMITATIONS OF LIABILITY AND TERMINATION PROVISIONS
c LOCATED DIRECTLY BELOW.

c WARRANTY DISCLAIMER.  WARF warrants to you that it owns the Software or otherwise has the right
c to grant the license granted to you in this Agreement.  WARF IS PROVIDING THE SOFTWARE TO YOU ON 
c AN “AS IS” BASIS.  WARF MAKES NO WARRANTY AS TO THE PERFORMANCE OF THE SOFTWARE OR ANY RESULTS YOU
c MAY OBTAIN BY USING THE SOFTWARE.  WARF MAKES NO REPRESENTATION AND EXTENDS NO WARRANTIES OF ANY 
c KIND, EITHER EXPRESS OR IMPLIED, AS TO TITLE OR NON-INFRINGEMENT OF THIRD-PARTY RIGHTS, 
c MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE.

c Limitation of Liability.  YOU ASSUME ALL RISK IN THE USE OF THE SOFTWARE, AND WARF SHALL INCUR 
c NO LIABILITY FOR YOUR USE OF THE SOFTWARE OR FOR ANY ACT OR EVENT RELATED TO OR RESULTING FROM 
c YOUR USE OF THE SOFTWARE.  TO THE FULLEST EXTENT PERMITTED BY LAW, IN NO EVENT SHALL WARF BE 
c LIABLE OR OBLIGATED UNDER CONTRACT, NEGLIGENCE, STRICT LIABILITY, CONTRIBUTION, BREACH OF 
c WARRANTY, OR OTHER LEGAL OR EQUITABLE THEORY, FOR ANY DAMAGES, COSTS OR EXPENSES OF ANY KIND 
c WHATSOEVER, HOWSOEVER CAUSED, OR ANY LOST PROFITS, OR LOST DATA, COST OF PROCUREMENT OF SUBSTITUTE
c GOODS, TECHNOLOGY, SERVICES, OR ANY CLAIMS BY THIRD PARTIES (INCLUDING BUT NOT LIMITED TO 
c ANY DEFENSE THEREOF), OR OTHER SIMILAR COSTS.  YOU ACKNOWLEDGE THAT THE FOREGOING LIMITATIONS 
c ARE REASONABLE AND AN ESSENTIAL PART OF THIS AGREEMENT.

c Termination.  Without prejudice to any other rights, WARF may terminate this Agreement if you 
c fail to comply with the terms and conditions of this Agreement.  Upon termination of this 
c Agreement for any reason, you immediately must destroy all copies of the Software in your 
c possession or control and must send WARF a written notice confirming destruction is complete.

c License Grant.  By clicking “OK,” WARF grants to you a non-exclusive, royalty-free, 
c non-sublicensable, non-transferable license, limited to the Licensed Field, to use the 
c Software and to internally display (e.g. post in a classroom) the Software solely for 
c academic research or non-profit scholarly purposes.  WARF expressly reserves all other rights.

c Intellectual Property.  WARF retains all right, title and interest in and to the Software, 
c including, without limitation, all intellectual property rights that may exist in the Software, 
c and all copies of the Software, regardless of form or medium, including the Copyright.  
c Nothing in this Agreement grants to you any additional rights to the Software or any updates, 
c new releases, commercial versions, or any other software program.

c Limitations on Use.  You may not take any of the following actions:  (1) modify, translate 
c or otherwise create any derivative works of the Software; (2) publicly display (e.g. post 
c on the Internet) or publicly perform (e.g. present at a conference) the Software; 
c (3) sell, lease, license, rent or lend the Software; and (4) make any other use of the Software 
c beyond the scope of the rights granted to you herein, including its use for any commercial 
c purposes whatsoever.  You must fully reproduce and not obscure, alter or remove any of WARF’s,
c or the inventors, copyright notices or other proprietary notices that appear on the Software.

c Submissions.  If you decide to transmit, by any means or media, any materials or other 
c information (including without limitation, any ideas, concepts or techniques for improvements 
c to the Software) to WARF, whether as feedback, data, questions, comments, suggestions or the like,
c you are providing such submissions to WARF on an unrestricted and non-confidential basis, and 
c you hereby grant WARF a non-exclusive, royalty-free, worldwide, perpetual, assignable, 
c irrevocable license, with the right to sublicense, to use, copy, transmit, distribute, create 
c derivative works of, publicly display and publicly perform such submissions in the Licensed Field.
c  You may contact WARF for any purpose under this Agreement or otherwise at:

c      Wisconsin Alumni Research Foundation

c      P.O. Box 7365, Madison, Wisconsin 53726

c      Telephone:  (608) 263-2500  Fax:  (608) 263-1064

c Law; Electronic Agreement.  The terms of this Agreement are governed by the laws of the 
c United States of America, and the State of Wisconsin, and you hereby submit to the nonexclusive
c jurisdiction of those courts for purposes of any dispute arising under this Agreement.  
c You acknowledge and agree to be bound by this Agreement by executing it in electronic form.  
c By selecting “YES,” you affirmatively agree to the terms of this Agreement.  
c If you select “NO,”you will not be permitted to install or use the Software.

c By selecting “OK,” you represent and warrant to WARF that you have read, understand and agree 
c to the terms of this Agreement and that you have full authority to enter into this Agreement 
c without obtaining any other approvals.  Please print a copy of this Agreement for your records.
 

c ©2004-2009 Wisconsin Alumni Research Foundation.   All rights reserved.

c---------------------------------------------------------------------------


	include 'tomoFDD.inc'

        use tomoFDD

        implicit none
        include 'RaySPDR.inc'

	integer		absolute_use
	real		acond
	real		adamp(20)
	real		adep
        real            air_dep
	integer		aiter(0:20)
	integer		ajoint(20)
	real		alat
	real		alon
	real		amaxdcc(20)
	real		amaxdct(20)
	real		amaxres_cross(20)
	real		amaxres_net(20)
	integer		amcusp(1000)
	real		awt_ccp(20)
	real		awt_ccs(20)
	real		awt_ctd(20)
	real		awt_ctp(20)
	real		awt_cts(20)
        integer         CC_format
	integer		clust(MAXCL,MAXEVE)
	real		cohav
	real		damp
	character	dattim*25
	real		dep_grid(0:MAXGRID)
	real		dtav
	integer		dt_c1(MAXDATA)
	integer		dt_c2(MAXDATA)
	real		dt_cal(MAXDATA)
	real		dt_dt(MAXDATA)
	integer		dt_ic1(MAXDATA)
	integer		dt_ic2(MAXDATA)
	integer		dt_idx(MAXDATA)
	integer		dt_ista(MAXDATA)
	real		dt_offs(MAXDATA)
	real		dt_qual(MAXDATA)
	real		dt_res(MAXDATA)
	character	dt_sta(MAXDATA)*7
	real		dv(MXPARI)
	real		dt_wt(MAXDATA)
	real		dxav
	real		dyav
	real		dzav
        real            dmav
        real            dthP(4), dthS(4)
	real		etav
        real            emav
	integer		ev_cusp(MAXEVE)
	integer		ev_date(MAXEVE)
	real		ev_dep(MAXEVE)
	real		ev_herr(MAXEVE)
	real		ev_lat(MAXEVE)
	real		ev_lon(MAXEVE)
	real		ev_mag(MAXEVE)
	real		ev_res(MAXEVE)
	integer		ev_time(MAXEVE)
        integer         ev_type(MAXEVE)
	real		ev_x(MAXEVE)
	real		ev_y(MAXEVE)
	real		ev_zerr(MAXEVE)
	real		ev_z(MAXEVE)
	integer         evID, evID1, evID2
        integer         staProj(MAXSTA)
	logical		ex
	real		exav
	real		eyav
	real		ezav
	integer		fd_nx
	integer		fd_ny
	integer		fd_nz
	character       fn_abs*80
	character	fn_cc*80
	character	fn_ct*80
	character	fn_eve*80
	character	fn_inp*80
	character	fn_loc*80
	character	fn_reloc*80
	character	fn_res*80
	character	fn_srcpar*80
	character	fn_sta*80
	character	fn_stares*80
	character	fn_vel*80
        character       fn_vp*80
        character       fn_vs*80
	real		finc
	integer		fu0
	integer		fu1
	integer		fu3
	integer		i, i2
	integer		iargc
	integer		ibeg
	integer		iclust
	integer		icusp(MAXEVE)
	integer		idata
	integer		idy
	integer		iend
        integer         ifindi
	integer		ihr
        integer         iicusp(MAXEVE)  ! [1..nev] Index table into ev_cusp[]
	integer		imn
	integer		imo
	integer	        in
	integer		ineg
	integer		iphase
	integer		isolv
        integer         isp
	integer		istart
	integer		iter
	integer		itf
	integer		iunit, iunit0, iunit1
	integer		iyr
	integer		j
	integer		jiter
	integer		joint
	integer		juliam
	integer		k
	integer		kiter
	integer		l
	doubleprecision	lat
	integer		log
	doubleprecision	lon
        integer         m
	real		maxdcc
	real		maxdct
	real 		maxdist
	integer		maxiter
        real            maxoff
	real		maxres_cross
	real		maxres_net
	integer		mbad
	integer		minobs_cc
	integer		minobs_ct
	real		minwght
	integer		n
	integer		narguments
	integer		ncc
	integer		nccold
	integer		nccp
	integer		nccs
	integer		nclust
	integer		nct
	integer		nctold
	integer		nctp
	integer		ncts
	integer		ncusp
	integer		ndt
        integer         ndtold
	integer		nev
	integer		nevold
	integer		niter
        integer         nn
	integer		noclust(MAXEVE)
	real		noisef_dt
        real            norm_threshold
	integer		nsrc
	integer		nsta
	real		picav
        integer         RayTracing ! choose which ray tracing method to use
	real		resvar1
	real		rms_cc
	real		rms_cc0
	real		rms_cc0old
	real		rms_ccold
	real		rms_ct
	real		rms_ct0
	real		rms_ct0old
	real		rms_ctold
        real            rota
	real		sc
	real		sdc0_dep
	real		sdc0_lat
	real		sdc0_lon
        integer         pgood(MAXOBS,MAXEVE)
        integer         sgood(MAXOBS,MAXEVE)
        real            slow
        real            dvel
	integer		src_cusp(MAXEVE)
	real		src_dep(MAXEVE)
	real		src_dt(MAXEVE)
	real		src_dx(MAXEVE)
	real		src_dy(MAXEVE)
	real		src_dz(MAXEVE)
	real		src_et(MAXEVE)
	real		src_ex(MAXEVE)
	real		src_ey(MAXEVE)
	real		src_ez(MAXEVE)
	real		src_lat0(MAXEVE)
	doubleprecision	src_lat(MAXEVE)
	real		src_lon0(MAXEVE)
	doubleprecision	src_lon(MAXEVE)
	integer		src_nnp(MAXEVE)
	integer		src_nns(MAXEVE)
	integer		src_np(MAXEVE)
	integer		src_ns(MAXEVE)
	real		src_rmsc(MAXEVE)
	real		src_rmsn(MAXEVE)
        integer         src_type(MAXEVE)
	real		src_t0(MAXEVE)
	real		src_t(MAXEVE)
	real		src_x0(MAXEVE)
	real		src_x(MAXEVE)
	real		src_y0(MAXEVE)
	real		src_y(MAXEVE)
	real		src_z0(MAXEVE)
	real		src_z(MAXEVE)
	real		sta_az(MAXSTA)
        integer         eve_sta(MAXEVE,MAXOBS+1)
	real		sta_dist(MAXSTA)
        integer         sta_itmp(MAXSTA)
	character	sta_lab(MAXSTA)*7
	real		sta_lat(MAXSTA)
	real		sta_lon(MAXSTA)
        real            sta_elv(MAXSTA)
	integer		sta_nnp(MAXSTA)
	integer		sta_nns(MAXSTA)
	integer		sta_np(MAXSTA)
	integer		sta_ns(MAXSTA)
	real		sta_rmsc(MAXSTA)
	real		sta_rmsn(MAXSTA)
	integer         staID
	real		stepl
	character	str1*60
	character	str80*80
	character	str3*3
	real		tav
        real            theta
        real            threshold(20)
	real		tdep, dep0, dep1
	real		tlat
	real		tlon
	real		tmpr1
	real		tmpr2
	real		tmp_ttp(MAXOBS,MAXEVE)
	real		tmp_tts(MAXOBS,MAXEVE)
	real		tmp_xp(MAXOBS,MAXEVE)
	real		tmp_yp(MAXOBS,MAXEVE)
	real		tmp_zp(MAXOBS,MAXEVE)
	real		tmp_xs(MAXOBS,MAXEVE)
	real		tmp_ys(MAXOBS,MAXEVE)
	real		tmp_zs(MAXOBS,MAXEVE)
	integer		trimlen
        integer         unknowns
        real            v
	real		weight1
	real		weight2
	real		weight3
	real		wlat
	real		wlon
	real		wt_ccp
	real		wt_ccs
	real		wt_ctp
	real		wt_cts
	real		wtdd
	real		xav
	real		xc
        real            xp
        real            xr
        real            xs
        real            yp
        real            yr
        real            ys
	real		yav
	real		yc
	real		zav
        real            zp
        real            zr
        real            zs
        real            zc
        real            vp(500), vs(500)

c--- some variables for lsfitHV_FDD_lsqr_shot
        real            src_em(MXPARI)
        integer         col(MAXND*MAXDATA)
        integer         total ! the number of nonzero velocity derivatives
        real            anorm
        real            arnorm
        real            atol
        real            btol
        real            conlim
        real            d(MAXDATA+4+MXPARI*6)   ! Data vector
        real            dtavold
        real            dxavold
        real            dyavold
        real            dzavold
        real            dmavold
        real            etavold
        real            exavold
        real            eyavold
        real            ezavold
        real            emavold
        real            factor
        integer         istop
        integer         itnlim
c       integer         iw(2*(8*MAXDATA+4*MAXEVE+MAXND*MAXDATA)+1 ) ! lsqr index array
        integer,allocatable :: iw(:)
        integer         k2
        integer         leniw
        integer         lenrw
        integer         nar
        integer         nndt
        real            norm(MAXEVE*4+MXPARI)
        real            norm_abs(MAXEVE*4+MXPARI)
        real            norm1(MAXEVE*4+MXPARI)
        real            norm_test(MAXEVE*4+MXPARI)
        real            rnorm, rnorm_wt
c       real            rw(8*MAXDATA+4*MAXEVE+MAXND*MAXDATA)
        real,allocatable :: rw(:)
        real            se(MAXEVE*4+MXPARI)     ! Solution error
        real            w1(MAXEVE*4+MXPARI)     ! Work space
        real            w2(MAXEVE*4+MXPARI)     ! Work space
        real            wtinv(MAXDATA+4+6*MXPARI)! +4 = mean shift constr
        real            wt(MAXDATA+4+6*MXPARI)
        real            x(MAXEVE*4+MXPARI)      ! Solution vector
        real            xnorm, xnorm_vel
        integer         old_index(4*MAXEVE+MXPARI)
        integer         new_index(4*MAXEVE+MXPARI)
        real            temp1(MXPARI), temp2(MXPARI)
        real            dx, dy, dz
        integer         m1, m2, m3, m4, m5
        integer         diff_node(MXPARI)
        integer         fix(4*MAXEVE)
        real            max_norm
        real            norm_node(MAXNX,MAXNY,2*MAXNZ)
        integer         new_npari
        integer         new_loc
        real            atemp
        real            lat0,lat1,lon0,lon1
        real            dlat,dlon,ddep

        real*8          aar, bbr, hr, aar0, bbr0, hr0
        real*8          aas, bbs, hs, aas0, bbs0, hs0        
        real*8          rp_geo(3,msg+1)
        real*8          dtemp, dtemp2        
        real*8          vp2,vs2,rho        
        real*8          dpi, r2d
        integer         nseg2
        real            rand,r0,tt
        integer         i1,j1,k1
        integer*8       iw_size, rw_size

        real*8 rak(npcom),vpak(npcom),vsak(npcom),dak(npcom)
        real*8 xlayer(0:nlayer)

       
        common /model/ xlayer,rak,vpak,vsak,dak


        write(*,*)"Regional-Scale Double-Difference Seismic Tomography (tomoDD-SE)" 
        write(*,*)"©2007-2015 Wisconsin Alumni Research Foundation."
        write(*,*)"All rights reserved."
        write(*,*)"By running this program, you represent and warrant to WARF that you have read,"
        write(*,*)"understand and agree to the terms of the Agreement enclosed in copyright.txt" 
        write(*,*)"and that you have full authority to enter into this Agreement without obtaining"
        write(*,*)"any other approvals."

	call alloc_FDD()
        iw_size = 2*(8*MAXDATA+4*MAXEVE+MAXND*MAXDATA)+1
        allocate( iw(iw_size) )
        if (.not.allocated(iw)) print *, "Memory allocation ERROR for iw"

        rw_size = 8*MAXDATA+4*MAXEVE+MAXND*MAXDATA
        allocate( rw(rw_size) )
        if (.not.allocated(rw)) print *, "Memory allocation ERROR for rw"

        print *, "Size of IW()=",iw_size
        print *, "Size of RW()=",rw_size

        r2d = 90./asin(1.)
        dpi = asin(1.)/ 90.
        r0  = 6371.0
 
        minwght= 0.00001
        rms_ccold= 0
        rms_ctold= 0
        rms_cc0old= 0
        rms_ct0old= 0
c--- open log file:
        call freeunit(log)
        open(log,file='tomoFDD.log',status='unknown')
        str1= 'starting tomoFDD (v1.0 - 01/2003)...'
        call datetime(dattim)
        write(6,'(a45,a)') str1, dattim
        write(log,'(a45,a)') str1, dattim

c--- get input parameter file name:
        narguments = iargc()
        if(narguments.lt.1) then
           write(*,'(/,a)') 'PARAMETER INPUT FILE [tomoFDD.inp <ret>]:'
           read(5,'(a)') fn_inp
           if(trimlen(fn_inp).le.1) then
              fn_inp= 'tomoFDD.inp' !default input file name
           else
              fn_inp= fn_inp(1:trimlen(fn_inp))
           endif
        else
           call getarg(1,fn_inp)
        endif
        inquire(FILE= fn_inp,exist=ex)
        if(.not. ex) stop' >>> ERROR OPENING INPUT PARAMETER FILE.'

c--- get input parameters:
c--- This function should change 
c--- define coordinate center and the size of the box (xd,yd)
        call getinpSPDR(MAXEVE,log,fn_inp,
     &       fn_cc,fn_ct,fn_sta,fn_eve, fn_abs, fn_vel,
     &       fn_loc,fn_reloc,fn_res,fn_stares,fn_srcpar,
     &       fn_vp,fn_vs,idata,iphase,
     &       minobs_cc,minobs_ct,
     &       amaxres_cross,amaxres_net,amaxdcc,amaxdct,
     &       noisef_dt,maxdist,
     &       awt_ccp,awt_ccs,awt_ctp,awt_cts,adamp, awt_ctd,
     &       istart,maxiter,isolv,niter,aiter,
     &       iclust,ncusp,icusp,wlat, wlon, 
     &       stepl, !RayTracing,
     &       rota, ajoint, CC_format, threshold,
     &       weight1, weight2, weight3, air_dep,
     &       lat0,lat1,lon0,lon1,dep0,dep1,
     &       dlat,dlon,ddep)

        write(*,*)'iuses:',iuses
        write(log,*)"Rotation angle:",rota
        !write(*,*)'lat0,lat1,lon0,lon1,dep0,dep1',lat0,lat1,lon0,lon1,dep0,dep1
        theta=rota*3.1415926/180.0
c---  The function of rota has not been added yet, will be updated later

c--- write computational parameters out
        write(log,'("Coordinate center---",/, 
     &        "Lat:",f9.4,"Lon:",f9.4)')wlat, wlon
        write(log,'("Computational step size: ", f9.4,/,
     &       "Ray-tracing step size: ",f9.4)') finc, stepl
        write(log,'("Model derivatives calculation parameters---",/,
     &       "iuses: ", i3, "invdel: ",i3, "iuseq=", i3)') iuses, invdel, iuseq

        open(3,file='MOD',status='old')
        open(16,file=fn_vel)
        open(15,file='grid.dat')

        open(10,file='ak135.15.SKS')
        rewind(10)

c* Read the layer division from an array outside:
        open(19,file=
     >      'layer-16.dat')

c NOTE AK135 contains: radius,depth,Vp,Vs,density

        do i=1,np15

          read(10,*)rak(i),dak(i),vp2,vs2,rho
          dak(i)=float(nint(dak(i)))
          rak(i)=float(nint(rak(i)))
          vpak(i)=vp2
          vsak(i)=vs2
        enddo
        close(10)

c***************************************************
c     * Read global layer division data:
        do i=0,nlayer
           read(19,*)xlayer(i)
        enddo
        close(19)

c---  read model definition parameters, modified from simul2000, subroutine 
c---  input3; Now MOD stores longitude, latitude and depth information!
c---  NOT X-Y-Z as in the previous version.
        call input_vel

        write(log,*)'# of inversion nodes:'
        write(log,*)'npari: ',npari,' nparvi: ',nparvi
        write(log,*)'nparpi:',nparpi,' nparsi: ',nparsi
        
c---  clustering:
        if((idata.eq.1.and.minobs_cc.eq.0).or.
     &       (idata.eq.2.and.minobs_ct.eq.0).or.
     &       (idata.eq.3.and.minobs_ct+minobs_cc.eq.0)) then
           absolute_use=1       ! reading absolute data if available
c---  get data:
           call getdata(
     &          log, fn_cc, fn_ct, fn_sta, fn_eve, fn_srcpar, fn_abs,
     &          idata, iphase, ncusp, icusp,
     &          maxdist,amaxdct(1),amaxdcc(1),CC_format,
     &          ev_date, ev_time, ev_cusp, ev_lat, ev_lon, ev_dep,
     &          ev_mag, ev_herr, ev_zerr, ev_res, ev_type,
     &          sta_lab, sta_lat, sta_lon, sta_elv,
     &          dt_sta, dt_dt, dt_qual, dt_c1, dt_c2, dt_idx,
     &          dt_ista, dt_ic1, dt_ic2,dt_offs,
     &          nev, nsta, ndt, nccp, nccs, nctp, ncts,
     &          absolute_use)

           nclust= 1
           clust(1,1)= nev
           do i=1,nev
              clust(1,i+1)= ev_cusp(i)
           enddo
           write(*,'(/,"no clustering performed.")')
           write(log,'(/,"no clustering performed.")')


c  we do not need clustering analysis for double-difference tomography          
        else
           absolute_use=0       ! do not use absolute data when clustering
c---  get data:
           write(*,*)wlat, wlon
           write(*,*)fd_nx, fd_ny, fd_nz
           call getdata(
     &          log, fn_cc, fn_ct, fn_sta, fn_eve, fn_srcpar, fn_abs,
     &          idata, iphase, ncusp, icusp,
     &          maxdist,amaxdct(1),amaxdcc(1), CC_format,
     &          ev_date, ev_time, ev_cusp, ev_lat, ev_lon, ev_dep,
     &          ev_mag, ev_herr, ev_zerr, ev_res, ev_type,
     &          sta_lab, sta_lat, sta_lon, sta_elv,
     &          dt_sta, dt_dt, dt_qual, dt_c1, dt_c2, dt_idx,
     &          dt_ista, dt_ic1, dt_ic2,dt_offs,
     &          nev, nsta, ndt, nccp, nccs, nctp, ncts,
     &          absolute_use)

           call clusterFDD(log, nev, ndt,
     &          idata, minobs_cc, minobs_ct,
     &          dt_c1, dt_c2, ev_cusp,
     &          clust, noclust, nclust)
           
        endif

c---  open files
        call freeunit(fu0)
        open(fu0,file=fn_loc,status='unknown')
        call freeunit(fu1)
        open(fu1,file=fn_reloc,status='unknown')
        if(trimlen(fn_stares).gt.1) then
           call freeunit(fu3)
           open(fu3,file=fn_stares,status='unknown')
        endif
        
        jiter = 0               ! counter for iter with no updating (air quakes)
c---  big loop over clusters starts here:
        if(iclust.ne.0) then
           if (iclust.lt.0 .or. iclust.gt.nclust) then
              write(*,*) 'error: invalid cluster number ',iclust
              write(*,*) 'must be between 1 and nclust (',nclust,')'
              stop
           endif
           ibeg= iclust
           iend= iclust
        else
           ibeg= 1
           iend= nclust
        endif
        do iclust= ibeg,iend
           call datetime(dattim)
           write(log,'(/,"RELOCATION OF CLUSTER:",i2,5x,a,/,
     &          "----------------------")')iclust,dattim
           write(*,'(/,"RELOCATION OF CLUSTER:",i2,5x,a,/,
     &          "----------------------")')iclust,dattim
           
c---  get data for each cluster if clustering was invoked:
           if((nclust.eq.1.and.clust(iclust,1).eq.nev).or.
     &          (idata.eq.1.and.minobs_cc.eq.0).or.
     &          (idata.eq.2.and.minobs_ct.eq.0).or.
     &          (idata.eq.3.and.minobs_ct+minobs_cc.eq.0)) goto 50
           
           ncusp= clust(iclust,1)
           do i=1,ncusp
              icusp(i)= clust(iclust,i+1)
           enddo
           
           absolute_use=1       ! read absolute data after clustering
           if(idata.ne.0) call getdata(
     &          log, fn_cc, fn_ct, fn_sta, fn_eve, fn_srcpar, fn_abs,
     &          idata, iphase, ncusp, icusp,
     &          maxdist,amaxdct(1),amaxdcc(1), CC_format,
     &          ev_date, ev_time, ev_cusp, ev_lat, ev_lon, ev_dep,
     &          ev_mag, ev_herr, ev_zerr, ev_res, ev_type,
     &          sta_lab, sta_lat, sta_lon, sta_elv,
     &          dt_sta, dt_dt, dt_qual, dt_c1, dt_c2, dt_idx,
     &          dt_ista, dt_ic1, dt_ic2,dt_offs,
     &          nev, nsta, ndt, nccp, nccs, nctp, ncts,
     &          absolute_use)
           
 50        continue
           nccold= nccp+nccs
           nctold= nctp+ncts
           ncc= nccp+nccs
           nct= nctp+ncts
           nevold= nev
           
c---  get cluster centroid:
           sdc0_lat= 0
           sdc0_lon= 0
           sdc0_dep= 0
           do i=1,nev
              sdc0_lat= sdc0_lat + ev_lat(i)
              sdc0_lon= sdc0_lon + ev_lon(i)
              sdc0_dep= sdc0_dep + ev_dep(i)
           enddo
           sdc0_lat= sdc0_lat/nev
           sdc0_lon= sdc0_lon/nev
           sdc0_dep= sdc0_dep/nev
           
           write(log,'("Cluster centroid at:",1x,f10.6,2x,f11.6,2x,f9.6)')
     &          sdc0_lat,sdc0_lon,sdc0_dep
           
c---  get cartesian coordinates for epicenters
           do i=1,nev
              tlat= ev_lat(i)
              tlon= ev_lon(i)
              tdep= ev_dep(i)
              call sph2car_ft(tlat,tlon,tdep,wlat,wlon,xc,yc,zc,theta)
              ev_x(i)= xc *1000
              ev_y(i)= yc *1000
              ev_z(i)= zc *1000
           enddo
           
           write(log,'("# events:",i5)')nev
           
c---  write output (mdat.loc):
           write(fu0,'(i9,1x,f10.6,1x,f11.6,1x,f9.3,1x,f10.1,1x,f10.1,
     &          1x,f10.1,
     &          1x,f8.1,1x,f8.1,1x,f8.1,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2,
     &          1x,f3.1,1x,i3)')
     &          (ev_cusp(i),ev_lat(i),ev_lon(i),ev_dep(i),ev_x(i),ev_y(i),
     &          ev_z(i),ev_herr(i)*1000,ev_herr(i)*1000,ev_zerr(i)*1000,
     &          int(ev_date(i)/10000),
     &          int(mod(ev_date(i),10000)/100),mod(ev_date(i),100),
     &          int(ev_time(i)/1000000),int(mod(ev_time(i),1000000)/10000),
     &          mod(real(ev_time(i)),10000.0)/100,ev_mag(i),iclust,i=1,nev)
           
c---  get initial trial sources:
           call trialsrc_FDD_shot(istart,sdc0_lat,sdc0_lon,sdc0_dep,
     &          nev,ev_cusp,ev_lat,ev_lon,ev_dep,ev_type,
     &          nsrc,src_cusp,src_lat0,src_lon0,
     &          src_x0,src_y0,src_z0,src_t0,
     &          src_lat,src_lon,src_dep, src_type,
     &          src_x,src_y,src_z,src_t,wlat,wlon,theta)
           
           write(*,'("Initial trial sources =",i6)')nsrc
           write(log,'("# initial trial sources:",i6)')nsrc

c--- determine first if all of sources are in the computation domain
           do i=1, nsrc
              if(src_lon(i).le.xn(1) .or. src_lon(i).ge.xn(nx) .or.
     &           src_lat(i).le.yn(1) .or. src_lat(i).ge.yn(ny) .or.
     &           src_dep(i).le.zn(1) .or. src_dep(i).ge.zn(nz) ) then
                 write(log,*)'Warning! Events outside the computation domain'
                 write(*,*)'Warning! Events outside the computation domain'
                 write(*,*)  src_cusp(i),src_lat(i),src_lon(i),src_dep(i)
                 write(log,*)src_cusp(i),src_lat(i),src_lon(i),src_dep(i)
              endif
           enddo

c--- output station cordinates
           do i=1,nsta
              tlat=sta_lat(i)
              tlon=sta_lon(i)
              tdep=sta_elv(i)              
              call  sph2car_ft(tlat,tlon,tdep,wlat,wlon,xc,yc,zc, theta)
c              write(*,*)i,xc,yc,zc
           enddo

c--- set up the array eve_sta that stores the station information per event

        do i=1, MAXEVE
           eve_sta(i,1)=0
        enddo

        do i=1, ndt
           staID=dt_ista(i)
           evID1=dt_ic1(i)
           evID2=dt_ic2(i)
           call add_sta(eve_sta,evID1,staID)
           call add_sta(eve_sta,evID2,staID)
        enddo


           
c---  loop over iterations starts here:
c     define iteration step at which re-weighting starts: this is dynam. since
c     it depends on the number of neg depths runs before.
           
c     first reset aiter() and maxiter
           do i=1,niter
              aiter(i) = aiter(i) - jiter
           enddo
           maxiter = maxiter - jiter
           
           kiter= 0		! counter for iter with data skipping
           jiter= 0		! counter for iter with no updating (air quakes)
           mbad= 0		! counter for air quakes
           
           iter= 1
 55        call datetime(dattim)
           write(log,'(/,"===ITERATION ",i3," (",i3,") ",a)')
     &          iter-jiter, iter, dattim
           
c---  get weighting parameters for this iteration:
           do i=1,niter
              if(iter.le.aiter(i)) goto 75
           enddo
 75        maxres_cross= amaxres_cross(i)
           maxres_net= amaxres_net(i)
           maxdcc= amaxdcc(i)
           maxdct= amaxdct(i)
           wt_ccp= awt_ccp(i)
           wt_ccs= awt_ccs(i)
           wt_ctp= awt_ctp(i)
           wt_cts= awt_cts(i)
           damp  = adamp(i)
           joint = ajoint(i)
           wtdd  = awt_ctd(i)
           norm_threshold = threshold(i)
           
           write(log, '(/,"Weighting parameters for this iteration:",/,
     &          "  wt_ccp= ",f7.4,2X,"wt_ccs= ",f7.4,2X,
     &          "maxr_cc= ",f7.4,2X,"maxd_cc= ",f7.2,2X,/,
     &          "  wt_ctp= ",f7.4,2x,"wt_cts= ",f7.4,2x,"maxr_ct= ",f7.4,2x,
     &          "maxd_ct= ",f7.2,/,"  damp= ",f7.1, 2x,"joint= ", i3,2x,
     &          "wtdd= ",f7.4,2x,"THRES= ",f9.5)')
     &          wt_ccp,wt_ccs,maxres_cross,
     &          maxdcc,wt_ctp,wt_cts,maxres_net,
     &          maxdct,damp, joint, wtdd, norm_threshold
           
c---  calculate travel times  and slowness vectors:
           write(log,'(/,"~ getting partials for ",i5,
     &          " stations and ",i5," source(s) ...")') nsta,nsrc

           if(mbad.gt.0) goto 1111

c---- To save memory usage, directly insert the codes from subroutine partials
           
c       Calculating hypocenter derivatives for P and S wave.
c       The rays are now traced using spherical pseudo-bending algorithm 

	iunit = 0
	if (trimlen(fn_srcpar).gt.1) then
c       Open source-parameter file
	   call freeunit(iunit)
	   open(iunit,file=fn_srcpar,status='unknown')
	endif
	if (iuses.eq.2) then	
	   write(16,*)' S wave will also be ray traced'
	endif

	do i=1,nsta
	   !write(*,*)i, sta_lab(i),sta_lat(i),sta_lon(i)
c       Calculate travel-time grid for station i, i.e. treating 
c       station i as a "source" first

           ! limit the range of staion's location ! Bug fixed by Hao Guo
           if (abs(sta_lat(i)).gt.90) sta_lat(i) = sta_lat(i)-sta_lat(i)/abs(sta_lat(i))*180
           if (abs(sta_lon(i)).gt.180) sta_lon(i) = sta_lon(i)-sta_lon(i)/abs(sta_lon(i))*360

	   tlat=sta_lat(i)
	   tlon=sta_lon(i)
	   tdep=sta_elv(i)
	  
           aar=tlat ! these must be double-precisions
           bbr=tlon
           hr =tdep

           aar0=tlat ! these must be double-precisions
           bbr0=tlon
           hr0 =tdep
 
	   do j=1, nsrc

cz---- change the program here, do not calculate the partial derivatives
cz---- for every station and event pair. Instead only calculates the 
cz---- partial derivatives between actual event and station pair. 
              
cz---- check if this event is recorded on this station
              staID=i
              evID=j
              call find_id(eve_sta,staID,evID,k)
              if(k.eq.0) goto 1199   ! do not calculate this station-event pair
              !write(*,*)src_cusp(i),src_lat(i),src_lon(i),src_dep(i)
c       convert source coordinates into km (HZ)
	     
              ! limit the range of event's location ! Bug fixed by Hao Guo
              if (abs(src_lat(j)).gt.90) src_lat(j) = src_lat(j)-src_lat(j)/abs(src_lat(j))*180
              if (abs(src_lon(j)).gt.180) src_lon(j) = src_lon(j)-src_lon(j)/abs(src_lon(j))*360 

	      tlat=src_lat(j)
	      tlon=src_lon(j)
	      tdep=src_dep(j)

              aas=tlat ! double-precisions
              bbs=tlon
              hs =tdep

              aas0=tlat ! double-precisions
              bbs0=tlon
              hs0 =tdep
	   
              aar=aar0
              bbr=bbr0
              hr =hr0
              !write(*,*)sta_lab(i),src_cusp(j),aar,bbr,hr,aas,bbs,hs

c---- Call the Spherical Pseudo-bending method
c---  ni is the number of segments in the ray path and w stores the
c---  coordinates of the segments. It starts from receive to the source
              isp=0

              call pbr(isp,aar, bbr, hr, aas, bbs, hs, rp_geo, nseg2, tmp_ttp(k,j))
              !write(16,*)"P ",sta_lab(i),src_cusp(j), tmp_ttp(k,j)
              !write(*,*)'calculated SPD time',tmp_ttp(k,j)
c             write(*,*)"nseg2:",nseg2
              !call rayl2(isp,rp_geo,nseg2,dtemp2)

              do i1=1, nseg2+1
                 tlat  = 90-rp_geo(2,i1)
                 !dtemp = 1.0/.99664719;
                 dtemp = 1.0;
                 tlat  = atan(dtemp*tan(tlat*dpi))*r2d
                 tlon  = rp_geo(3,i1)
                 tdep  = r0-rp_geo(1,i1)
                 call  sph2car_ft(tlat,tlon,tdep,wlat,wlon,xc,yc,zc,theta)
                 rp(1,i1)=xc
                 rp(2,i1)=yc
                 rp(3,i1)=zc
                 !write(*,*)i1,tlat,tlon,tdep,xc,yc,zc
              enddo
              nrp = nseg2+1
              pgood(k,j) = 1

c       determine the hypocenter derivatives for P-wave
	      call ttmder(0, dthP, stepl, wlat, wlon, theta,tt)

	      if(src_type(evID).ne.0) then !shot or Blast data
c		 write(*,*)'Shot data....'
		 tmp_xp(k,j)=0
		 tmp_yp(k,j)=0
		 tmp_zp(k,j)=0
	      else ! earthquake data
		 tmp_xp(k,j)=-dthP(2)
		 tmp_yp(k,j)=-dthP(3)
		 tmp_zp(k,j)=-dthP(4)
	      endif
              
              !write(*,*)tmp_ttp(k,j),tt,tmp_xp(k,j),tmp_yp(k,j),tmp_zp(k,j)
	      n=0
	      do m=1,npari  
		 if (abs(dtm(m)).ge.0.05) then
		    n=n+1
		    tmp_vp_index(k,j,n+1)=m
		    tmp_vp(k,j,n)=-dtm(m)	
		 endif
	      enddo
	      tmp_vp_index(k,j,1)=n ! No. of nonzero model derivatives
	      
c       calculate the S-wave travel time by 3D ray tracing
	      
	      if (iuses.eq.2) then
		
                 aas=aas0 ! double-precisions
                 bbs=bbs0
                 hs =hs0

                 aar=aar0
                 bbr=bbr0
                 hr =hr0
c                write(*,*)aar,bbr,hr,aas,bbs,hs

c---- Call the Spherical Pseudo-bending method
c---  ni is the number of segments in the ray path and w stores the
c---  coordinates of the segments. It starts from receive to the source
                 isp=1
                 call pbr(isp,aar, bbr, hr, aas, bbs, hs, rp_geo, nseg2, tmp_tts(k,j))
                 !write(16,*)'S ',sta_lab(i), src_cusp(j), tmp_tts(k,j)
c                write(*,*)"nseg2:",nseg2
                 !call rayl2(isp,rp_geo,nseg2,dtemp2)

                 do i1=1, nseg2+1
                   tlat  = 90-rp_geo(2,i1)
                   !dtemp = 1.0/.99664719;
                   dtemp = 1.0
                   tlat  = atan(dtemp*tan(tlat*dpi))*r2d
                   tlon  = rp_geo(3,i1)
                   tdep  = r0-rp_geo(1,i1)
                   call sph2car_ft(tlat,tlon,tdep,wlat,wlon,xc,yc,zc,theta)
                   rp(1,i1)=xc
                   rp(2,i1)=yc
                   rp(3,i1)=zc
                 enddo
                 nrp = nseg2+1
                 sgood(k,j) = 1

		 call ttmder(1, dthS, stepl, wlat, wlon, theta)

		 if(src_type(evID).ne.0) then !shot or blast data
		    tmp_xs(k,j)=0
		    tmp_ys(k,j)=0
		    tmp_zs(k,j)=0 
		 else ! Earthquake data
		    tmp_xs(k,j)=-dthS(2)
		    tmp_ys(k,j)=-dthS(3)
		    tmp_zs(k,j)=-dthS(4) 
		 endif
	
		 n=0
		 do m=1,npari
		    if (abs(dtm(m)).ge.0.05) then
		       n=n+1
		       tmp_vs_index(k,j,n+1)=m
		       tmp_vs(k,j,n)=-dtm(m)
		    endif
		 enddo
		 tmp_vs_index(k,j,1)=n ! save the number of non model derivatives
		 
	      endif
	      
10	      continue

c           Write to source-parameter file
	      if (iunit .ne. 0)
     &         write(iunit,'(i9,2x,f9.4,2x,f9.4,2x,a7,2x,f9.4)')
     &         src_cusp(j), src_lat(j), src_lon(j), sta_lab(i)
     
1199          continue
	   enddo
	enddo
c       write(*,*)'Successful in ending partials_FDD'

	if (iunit .ne. 0) close(iunit) ! Source-parameter file

1111    continue ! in the case of airquake, directly jump to here
	
c---  get double difference vector:
           call dtres_FDD(log,ndt,MAXSTA,nsrc,nsta,
     &          dt_dt,dt_idx,
     &          dt_ista,dt_ic1,dt_ic2,
     &          src_cusp,src_t,tmp_ttp,tmp_tts,
     &          dt_cal,dt_res, eve_sta)
           
c---  get a priori weights and reweight residuals
           call weighting_FDD(log,ndt,mbad,amcusp,idata,kiter,ineg,
     &          maxres_cross,maxres_net,maxdcc,maxdct,minwght,
     &          wt_ccp,wt_ccs,wt_ctp,wt_cts, wtdd,
     &          dt_c1,dt_c2,dt_idx,dt_qual,dt_res,dt_offs,
     &          dt_wt)
           
c---  skip outliers and/or air quakes:
           if(ineg.gt.0) then

c--- To save memory usage, directly insert skipping code here

      write(log,'("skipping data...")')

c     Skip data with large resiudals
      if (kiter.eq.1) then
          ndtold = ndt
          nccold = 0
          nctold = 0
      endif
      ncc = 0
      nct = 0
      j = 1
      do i=1,ndt
         if (kiter.eq.1) then
            if (dt_idx(i).le.2) then
               nccold = nccold+1
            else
               nctold = nctold+1
            endif
         endif
         if (dt_wt(i).ge.minwght) then
            dt_sta(j) = dt_sta(i)
            dt_c1(j) = dt_c1(i)
            dt_c2(j) = dt_c2(i)
            dt_idx(j) = dt_idx(i)
            dt_qual(j) = dt_qual(i)
            dt_dt(j) = dt_dt(i)
            dt_cal(j) = dt_cal(i)
            dt_res(j) = dt_res(i)
            dt_wt(j) = dt_wt(i)
            dt_offs(j) = dt_offs(i)
            if (dt_idx(i).le.2) then
                ncc = ncc+1
            else
                nct = nct+1
            endif
            j = j+1
         endif
      enddo
      ndt = j-1
      write(log,'("# obs = ",i9," (",f5.1,"%)")')
     &ndt, (ndt*100.0/ndtold)
      if (nccold.gt.0.and.nctold.gt.0) then
         write(log,'("# obs cc = ",i9," (",f5.1,"%)")')
     &   ncc, (ncc*100.0/nccold)
         write(log,'("# obs ct = ",i9," (",f5.1,"%)")')
     &   nct, (nct*100.0/nctold)
      endif

c     Skip events
      do i=1,ndt
         dt_ic1(i) = dt_c1(i)   !dt_ic1 is just a workspace array here!
         dt_ic2(i) = dt_c2(i)   !dt_ic2 is just a workspace array here!
      enddo
      call sorti(ndt, dt_ic1)
      call sorti(ndt, dt_ic2)
      k = 1
      do i=1,nev
         if (ifindi(ndt, dt_ic1, ev_cusp(i)).gt.0 .or.
     &       ifindi(ndt, dt_ic2, ev_cusp(i)).gt.0) then
            ev_date(k) = ev_date(i)
            ev_time(k) = ev_time(i)
            ev_cusp(k) = ev_cusp(i)
            ev_lat(k) = ev_lat(i)
            ev_lon(k) = ev_lon(i)
            ev_dep(k) = ev_dep(i)
            ev_mag(k) = ev_mag(i)
            ev_herr(k) = ev_herr(i)
            ev_zerr(k) = ev_zerr(i)
            ev_res(k) = ev_res(i)
	    ev_type(k)=ev_type(i)
            ev_x(k) = ev_x(i)
            ev_y(k) = ev_y(i)
            ev_z(k) = ev_z(i)
            k = k+1
         endif
      enddo
      nev = k-1
      write(log,'("# events = ",i9)') nev

c     Skip sources
c     Uses sorted dt_ic[12] arrays from above
      if (nsrc.ne.1) then
         k = 1
         do i=1,nsrc
            if (ifindi(ndt, dt_ic1, src_cusp(i)).gt.0 .or.
     &          ifindi(ndt, dt_ic2, src_cusp(i)).gt.0) then
               src_cusp(k) = src_cusp(i)
               src_lat(k) = src_lat(i)
               src_lon(k) = src_lon(i)
               src_lat0(k) = src_lat0(i)
               src_lon0(k) = src_lon0(i)
               src_dep(k) = src_dep(i)
               src_x(k) = src_x(i)
               src_y(k) = src_y(i)
               src_z(k) = src_z(i)
               src_t(k) = src_t(i)
               src_x0(k) = src_x0(i)
               src_y0(k) = src_y0(i)
               src_z0(k) = src_z0(i)
               src_t0(k) = src_t0(i)
	       src_type(k)=src_type(i)
               do j=1,eve_sta(i,1)
                  tmp_ttp(j,k) = tmp_ttp(j,i)
                  tmp_tts(j,k) = tmp_tts(j,i)
                  tmp_xp(j,k) = tmp_xp(j,i)
                  tmp_yp(j,k) = tmp_yp(j,i)
                  tmp_zp(j,k) = tmp_zp(j,i)
                  tmp_xs(j,k) = tmp_xs(j,i)
                  tmp_ys(j,k) = tmp_ys(j,i)
                  tmp_zs(j,k) = tmp_zs(j,i)
                  pgood(j,k)  = pgood(j,i)
                  do m=1, tmp_vp_index(j,i,1)
                     tmp_vp(j,k,m)=tmp_vp(j,i,m)
                     tmp_vp_index(j,k,m+1)=tmp_vp_index(j,i,m+1)
                  enddo
                  tmp_vp_index(j,k,1)=tmp_vp_index(j,i,1)
                  if (iuses.eq.2) then
                     sgood(j,k)=sgood(j,i)
                     do m=1, tmp_vs_index(j,i,1)
                        tmp_vs(j,k,m)=tmp_vs(j,i,m)
                        tmp_vs_index(j,k,m+1)=tmp_vs_index(j,i,m+1)
                     enddo
                     tmp_vs_index(j,k,1)=tmp_vs_index(j,i,1)
                  endif
                  eve_sta(k,j+1) = eve_sta(i,j+1)
               enddo
               eve_sta(k,1) = eve_sta(i,1)
               k = k+1
            endif
         enddo
         nsrc = k-1
      endif

c    Clean stations
      do i=1,nsta
         sta_itmp(i) = 0
      enddo
      do j=1,ndt
         do i=1,nsta
            if (dt_sta(j).eq.sta_lab(i)) then
               sta_itmp(i) = 1
               goto 200	! break
            endif
         enddo
200      continue
      enddo

c--- initialize the station projector
      do i=1,nsta
         staProj(i)=0 ! this indicates the station is removed
      enddo

      k = 1
      do i=1,nsta
         if (sta_itmp(i).eq.1) then
            sta_lab(k) = sta_lab(i)
            sta_lat(k) = sta_lat(i)
            sta_lon(k) = sta_lon(i)
            sta_dist(k) = sta_dist(i)
            sta_az(k) = sta_az(i)
            sta_np(k) = sta_np(i)
            sta_ns(k) = sta_ns(i)
            sta_nnp(k) = sta_nnp(i)
            sta_nns(k) = sta_nns(i)
            sta_rmsc(k) = sta_rmsc(i)
            sta_rmsn(k) = sta_rmsn(i)
c---- set up the station projector
            staProj(i) = k
            k = k+1
         endif
      enddo
      nsta = k-1
      write(log,'("# stations = ",i9)') nsta

cz--- re-organize several arrays
       do j=1,nsrc
          k=1
          do i=1,eve_sta(j,1)
             n=eve_sta(j,i+1)
             if(staProj(n).gt.0) then
               tmp_ttp(k,j) = tmp_ttp(i,j)
               tmp_tts(k,j) = tmp_tts(i,j)
               tmp_xp(k,j) = tmp_xp(i,j)
               tmp_yp(k,j) = tmp_yp(i,j)
               tmp_zp(k,j) = tmp_zp(i,j)
               tmp_xs(k,j) = tmp_xs(i,j)
               tmp_ys(k,j) = tmp_ys(i,j)
               tmp_zs(k,j) = tmp_zs(i,j)
               pgood(k,j)  = pgood(i,j)
               do m=1, tmp_vp_index(i,j,1)
                  tmp_vp(k,j,m)=tmp_vp(i,j,m)
                  tmp_vp_index(k,j,m+1)=tmp_vp_index(i,j,m+1)
               enddo
               tmp_vp_index(k,j,1)=tmp_vp_index(i,j,1)
               if (iuses.eq.2) then
                  sgood(k,j)=sgood(i,j)
                  do m=1, tmp_vs_index(i,j,1)
                     tmp_vs(k,j,m)=tmp_vs(i,j,m)
                     tmp_vs_index(k,j,m+1)=tmp_vs_index(i,j,m+1)
                  enddo
                  tmp_vs_index(k,j,1)=tmp_vs_index(i,j,1)
               endif
               eve_sta(j,k+1)=staProj(eve_sta(j,i+1))
               k=k+1
             endif
           enddo
           eve_sta(j,1) = k-1
       enddo



c     Index station labels and cuspids
      call indexxi(nev, ev_cusp, iicusp)
      do i=1,nev
         icusp(i) = ev_cusp(iicusp(i)) !icusp is just a workspace array here!
      enddo
      do i=1,ndt
         do j=1,nsta
            if (dt_sta(i).eq.sta_lab(j)) then
               dt_ista(i) = j
               dt_ic1(i) = iicusp(ifindi(nev, icusp, dt_c1(i)))
               dt_ic2(i) = iicusp(ifindi(nev, icusp, dt_c2(i)))
               goto 300	! continue 2
            endif
         enddo
         write(*,'("FATAL ERROR (indexing).")')
         stop   
300      continue
      enddo

c              write(*,*)'After skiping'
c              do i=1,nsrc
c                 do j=1,eve_sta(i,1)
c                    write(16,*)sta_lab(eve_sta(i,j+1)),src_cusp(i)
c                    write(16,*)(tmp_vp_index(j,i,k+1),k=1,tmp_vp_index(j,i,1))
c                 enddo
c              enddo
c--   Dont mess anymore with this cluster if we have wiped out all events
              if(nev.lt.2) then 
                 write(log,*)' Cluster has less than 2 events.'
                 write(*,*)' Cluster has less than 2 events.'
                 goto 778
              endif
           else
              write(log,'("no data skipped.")')
              
           endif
           
           
c---  least square fitting:
           
           if(isolv.eq.2) then
              if(joint .eq. 1)  then ! joint inversion
                 write(*,*) 'Preparing joint inversion....'
                 
                 if(iter.eq.1) then
                    resvar1= -999
                    unknowns=4*nev+npari
                    call resstat_FDD(log,idata,ndt,unknowns,dt_res,dt_wt,dt_idx,
     &                   rms_cc,rms_ct,rms_cc0,rms_ct0,
     &                   rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
     &                   resvar1)
                 endif


C       New version, which can deal well with column scaling associate with
c       model derivatives
c       also add smooth constraint into the inversion
c       add function to deal with the case when the absolute data are not
c       used any more

      write(log,'(/,"~ setting up G matrix.. ")')

c     If mean shift not contstrained

      m = 0   ! indicates the nonzero model derivatives
c     Prepare sparse data and design vector

c--- Absolute_use parameter can be controled by wtdd
c----if wtdd=0 and reweighting is applied, then absolute data are removed from inversion
c----Use Norm_threshold to control if all the nodes are inverted
	
c----Initialize weigh matrix
      do i=1,MAXDATA+4+6*MXPARI
	 wt(i)=0.0
	 wtinv(i)=0.0
      enddo

      do i=1,ndt
c        Weight data first
         wt(i) = dt_wt(i)

         if (wt(i).ne.0) then
            wtinv(i) = 1.0/wt(i)
         else
            wtinv(i) = 1.0
         endif
c         if (abs(wt(i)).gt.1) write(log,*)'Warning! Weight=',wt(i)
         !write(*,*) "dt_res(i)",dt_res(i),dt_ic1(i),dt_ic2(i)
         d(i) = dt_res(i)*1000.0 * wt(i)
         iw(1+i) = i
         iw(1+ndt+i) = i
         iw(1+2*ndt+i) = i
         iw(1+3*ndt+i) = i
         iw(1+4*ndt+i) = i
         iw(1+5*ndt+i) = i
         iw(1+6*ndt+i) = i
         iw(1+7*ndt+i) = i

c        Set up non-zero G matrix elements and apply weights
         if (nsrc.eq.1) then
            k1 = 1
            k2 = 1
         else
            call find_id(eve_sta,dt_ista(i),dt_ic1(i),k1)
            if(k1.eq.0) write(*,*)sta_lab(dt_ista(i)),src_cusp(dt_ic1(i))
            call find_id(eve_sta,dt_ista(i),dt_ic2(i),k2)
            if(k2.eq.0) write(*,*)sta_lab(dt_ista(i)),src_cusp(dt_ic2(i))
         endif

         if (dt_idx(i).eq.2 .or. dt_idx(i).eq.4) then ! S data
c--- do not use the event-station pair when it fails to calculate the
c--- ray path
            if(sgood(k1,dt_ic1(i)).eq.0.or.sgood(k2,dt_ic2(i)).eq.0) then
               wt(i)=0
               d(i)=0
            endif

            rw(i)       = tmp_xs(k1,dt_ic1(i)) * wt(i) 
            rw(ndt+i)   = tmp_ys(k1,dt_ic1(i)) * wt(i) 
            rw(2*ndt+i) = tmp_zs(k1,dt_ic1(i)) * wt(i) 

	    if(src_type(dt_ic1(i)).eq.1) then ! shot data
	       rw(3*ndt+i)= 0.0
	    else  !earthquake (0) or Blast (2)
	       rw(3*ndt+i) = wt(i)
	    endif
	    		  
	    if(dt_ic1(i).ne.dt_ic2(i)) then ! difference data            
                rw(4*ndt+i) = -tmp_xs(k2,dt_ic2(i)) * wt(i) 
                rw(5*ndt+i) = -tmp_ys(k2,dt_ic2(i)) * wt(i) 
                rw(6*ndt+i) = -tmp_zs(k2,dt_ic2(i)) * wt(i) 
		if(src_type(dt_ic2(i)).eq.1) then ! shot data
		   rw(7*ndt+i) = 0.0
		else ! earthquake (0) or blast data (2)
		   rw(7*ndt+i) = -wt(i)	   
		endif
	        do k=1, npari
	           temp1(k)=0.0
	           temp2(k)=0.0
	        enddo

C           number of non-zero derivatives in Vs 
		do k=1, tmp_vs_index(k1,dt_ic1(i),1) 
		   j=tmp_vs_index(k1,dt_ic1(i),k+1)
		   temp1(j)=tmp_vs(k1,dt_ic1(i),k) ! model derivatives
		enddo

		do k=1, tmp_vs_index(k2,dt_ic2(i),1)
		   j=tmp_vs_index(k2,dt_ic2(i),k+1)
		   temp2(j)=tmp_vs(k2,dt_ic2(i),k)
		enddo
		do k=1, npari
		   if(abs(temp1(k)-temp2(k)).gt.0.01) then ! 0.01 could be changed
		      m=m+1
		      rw(8*ndt+m)=(temp1(k)-temp2(k))*wt(i)
		      iw(1+8*ndt+m)=i
		      col(m)=k
		   endif
		   
		enddo
		
	     else		! absolute travel time
		   rw(4*ndt+i) = 0.0
		   rw(5*ndt+i) = 0.0
		   rw(6*ndt+i) = 0.0
		   rw(7*ndt+i) = 0.0	       
		   do k=1,tmp_vs_index(k1,dt_ic1(i),1)
		      if(abs(tmp_vs(k1,dt_ic1(i),k)).gt.0.01) then
			 m=m+1
			 rw(8*ndt+m)=tmp_vs(k1,dt_ic1(i),k)*wt(i)
			 iw(1+8*ndt+m)=i
			 col(m)=tmp_vs_index(k1,dt_ic1(i),k+1)	
		      endif
		   enddo
	
	     endif		! end of S wave
	     
	  else			! if it is P-wave
	     
             if(pgood(k1,dt_ic1(i)).eq.0.or.pgood(k2,dt_ic2(i)).eq.0) then
                wt(i)=0
                d(i)=0
             endif

	     rw(i)       = tmp_xp(k1,dt_ic1(i)) * wt(i)
	     rw(ndt+i)   = tmp_yp(k1,dt_ic1(i)) * wt(i)
	     rw(2*ndt+i) = tmp_zp(k1,dt_ic1(i)) * wt(i)

	     if(src_type(dt_ic1(i)).eq.1) then ! shot data
		rw(3*ndt+i) = 0.0
	     else ! earthquake (0) or blast data (2)
		rw(3*ndt+i) = wt(i)
	     endif
	     
	     if(dt_ic1(i).ne.dt_ic2(i)) then	! double difference data

		rw(4*ndt+i) = -tmp_xp(k2,dt_ic2(i)) * wt(i)
		rw(5*ndt+i) = -tmp_yp(k2,dt_ic2(i)) * wt(i)
		rw(6*ndt+i) = -tmp_zp(k2,dt_ic2(i)) * wt(i)

		if(src_type(dt_ic2(i)).eq.1) then !shot data
		   rw(7*ndt+i) = 0
		else ! earthquake (0) or Blast data (2)
		   rw(7*ndt+i) = -wt(i)
		endif

		do k=1, npari
		   temp1(k)=0.0
		   temp2(k)=0.0
		enddo
		do k=1, tmp_vp_index(k1,dt_ic1(i),1)
		   j=tmp_vp_index(k1,dt_ic1(i),k+1)
		   temp1(j)=tmp_vp(k1,dt_ic1(i),k)
		enddo
		do k=1, tmp_vp_index(k2,dt_ic2(i),1)
		   j=tmp_vp_index(k2,dt_ic2(i),k+1)
		   temp2(j)=tmp_vp(k2,dt_ic2(i),k)
		enddo
		do k=1, npari
		   if(abs(temp1(k)-temp2(k)).gt.0.01) then
		      m=m+1
		      rw(8*ndt+m)=(temp1(k)-temp2(k))*wt(i)
		      iw(1+8*ndt+m)=i
		      col(m)=k
		   endif
		   
		enddo
		
	     else		! absolute travel time
		   rw(4*ndt+i) = 0.0
		   rw(5*ndt+i) = 0.0
		   rw(6*ndt+i) = 0.0
		   rw(7*ndt+i) = 0.0	       
		   
		   do k=1,tmp_vp_index(k1,dt_ic1(i),1)
		      if (abs(tmp_vp(k1,dt_ic1(i),k)) .gt.0.01) then
			 m=m+1
			 rw(8*ndt+m)=tmp_vp(k1,dt_ic1(i),k)*wt(i)
			 iw(1+8*ndt+m)=i
			 col(m)=tmp_vp_index(k1,dt_ic1(i),k+1)	
		      endif
		   enddo
		
	     endif		! end of P-wave
	     
	  endif
	enddo
	
	write(*,*)'# of nonzero model derivatives before smoothing', m
        if(m.gt.MAXND*MAXDATA) then
           write(*,*)'Error: M>MAXND*MAXDATA'
           stop
        endif
	

cz The following section will be overlaped after the smooting constarint
cz--------------------------------------------------------------
	nar=m+8*ndt

	do i=1, ndt
c       Set up column indexes with non-zero elements
	   iw(1+nar+      i) = 4*dt_ic1(i) - 3
	   iw(1+nar+  ndt+i) = 4*dt_ic1(i) - 2
	   iw(1+nar+2*ndt+i) = 4*dt_ic1(i) - 1
	   iw(1+nar+3*ndt+i) = 4*dt_ic1(i)
	   iw(1+nar+4*ndt+i) = 4*dt_ic2(i) - 3
	   iw(1+nar+5*ndt+i) = 4*dt_ic2(i) - 2
	   iw(1+nar+6*ndt+i) = 4*dt_ic2(i) - 1

	   iw(1+nar+7*ndt+i) = 4*dt_ic2(i) 
	enddo

	do i=1,m
	   iw(1+nar+8*ndt+i)=col(i)+nev*4
	enddo
cz------------------------------------------------------------
	
cz    determine which model derivatives are not zeroes.
cz    later only invert these models ! 

	do i=1,4*nev+npari
	   norm(i) = 0.0
	   norm_abs(i)=0.0
	enddo

	do i=1,nar
	   norm(iw(1+nar+i)) = norm(iw(1+nar+i)) + rw(i)**2
	   norm_abs(iw(1+nar+i)) = norm_abs(iw(1+nar+i)) + abs(rw(i))
	enddo

cz Add on Feb 25, 2003
cz find the maximum absolute value of derivative weight sum
	max_norm=-999999.0
	
c	write(16,*)(norm(i),i=1,4*nev)
	do i=nev*4+1,nev*4+npari	   
	   if( norm_abs(i).gt.max_norm ) max_norm=norm_abs(i)
	enddo

	write(*,*)'max_norm=',max_norm

cz output norm(i) into unit 16 to know its distributions
	write(16,*)'DWS for P-wave.'

cz initialize norm_node
      do k=1,2*nz
	 do j=1,ny
	    do i=1,nx
	       norm_node(i,j,k)=0.0
	    enddo
	 enddo
      enddo

c--- check the values here
      !write(*,*)'nxy2=',nxy2
      !write(*,*)'nx2=',nx2
      
      do n=1, npari
         nn=mdexfx(n) 
         !write(16,*)'n,nn:',n,nn        
c        calculate x and z indices of grid
         k=(nn-1)/nxy2+2
         j=2+(nn-1+(2-k)*nxy2)/nx2
         i=1+nn+nx2*(2-j)+nxy2*(2-k)   
         
         if (k.ge.nz) k=k+2    ! if s node        
         norm_node(i,j,k)=norm_abs(n+4*nev)
	enddo

      do k=1,nz
         do j=1,ny
            write(16, '(100f10.1)')(norm_node(i,j,k),i=1,nx)
         enddo
      enddo

      if(iuses.eq.2) then
	 write(16,*)'DWS for S-wave.'
	 do k=nz+1,nz*2
	    do j=1,ny
	       write(16, '(100f10.1)')(norm_node(i,j,k),i=1,nx)
	    enddo
	 enddo
      endif

      write(16,*)'End of the norm node'
	
cz       judge if the norm(i) is larger than the norm_threshold of the maximum value.
cz       Then keep those nodes. Consider only model derivatives.
cz       A suggestion is that let norm_threshold equal to zero whe Abs. 
cz       data are more weighted.

c--- Check if there is any zero for location and origin time
c---- derivatives. This will be the case if shot and blast data 
c---- are included.

	j=0
	do i=1,4*nev
	   fix(i)=0 ! fix location or origin time
	   if(norm_abs(i).gt.0) then
	      j=j+1
	      old_index(j)=i
	      new_index(i)=j
	      fix(i)=1
	   endif
	enddo
	new_loc=j

	write(*,*)'4*nev=',4*nev
	write(*,*)'# of nonzero location and origin derivatives:',new_loc

	do i=nev*4+1,nev*4+npari
            if(i.le.nev*4+nparpi) then ! P node
	      if (norm_abs(i).ge.norm_threshold) then
		 diff_node(i-4*nev)=1
		 j=j+1
		 old_index(j)=i
		 new_index(i)=j
	      else
		 diff_node(i-4*nev)=0 ! fix the node
	      endif
            else ! S node
              if (norm_abs(i).ge.norm_threshold/3.0) then
                 diff_node(i-4*nev)=1
                 j=j+1
                 old_index(j)=i
                 new_index(i)=j
              else
                 diff_node(i-4*nev)=0 ! fix the node
              endif
            endif
	enddo

	new_npari=j-new_loc

	write(*,*)'new_loc=', new_loc, '   new_npari=',new_npari

c       add smoothing contraint into the inversion
c       only apply the smoothing constraint to the models with effective
c       derivatives
	n=0
	if(weight1.eq.0.and.weight2.eq.0.and. weight3.eq.0.0) then
	   write(*,*)'No smoothing constraint applied'
	else
	   write(*,*)'Applying smoothing constraint!'
	   do k=2, nz-2		! P velocity node
              !if(zn(k).ge.80) then
                    !weight1=2.0*weight1
                    !weight2=2.0*weight2
                    !weight3=1.2*weight3
              !endif

              !if(zn(k).ge.60) then
              !      weight1=1.5*weight1
              !      weight2=1.5*weight2
              !      weight3=1.2*weight3
              !endif


	      do j=2, ny-2
		 do i=2, nx-2
		    m1=(k-2)*nxy2+(j-2)*nx2+i-1
		    m2=(k-2)*nxy2+(j-2)*nx2+i+1-1
		    m3=(k-2)*nxy2+(j+1-2)*nx2+i-1
		    m4=(k+1-2)*nxy2+(j-2)*nx2+i-1
		    m1=ndexfx(m1)
		    m2=ndexfx(m2)
		    m3=ndexfx(m3)
		    m4=ndexfx(m4)
		    dx=abs(xn(i)-xn(i+1))
		    dy=abs(yn(j)-yn(j+1))
		    dz=abs(zn(k)-zn(k+1))	
		    
		    if(diff_node(m1).eq.1.and.diff_node(m2).eq.1) then
		       n=n+1
		       m=m+1
		       rw(8*ndt+m)=weight1 ! in X direction
		       iw(1+8*ndt+m)=ndt+n ! row index
		       col(m)=m1
		    
		       m=m+1
		       rw(8*ndt+m)=weight1*(-1.0) ! in X direction
		       iw(1+8*ndt+m)=ndt+n ! row index
		       col(m)=m2
		       d(ndt+n)=0.0
		    endif
		    
		    if(diff_node(m1).eq.1.and.diff_node(m3).eq.1) then
		       n=n+1
		       m=m+1
		       rw(8*ndt+m)=weight2*1.0 ! in Y direction
		       iw(1+8*ndt+m)=ndt+n ! row index
		       col(m)=m1
		    
		       m=m+1
		       rw(8*ndt+m)=weight2*(-1.0) ! in Y direction
		       iw(1+8*ndt+m)=ndt+n ! row index
		       col(m)=m3
		       d(ndt+n)=0.0
		    endif

		    if(diff_node(m1).eq.1.and.diff_node(m4).eq.1) then
		       m=m+1
		       n=n+1
		       rw(8*ndt+m)=weight3*1.0 ! in Z direction
		       iw(1+8*ndt+m)=ndt+n ! row index
		       col(m)=m1
		       
		       m=m+1
		       rw(8*ndt+m)=weight3*(-1.0) ! in Z direction
		       iw(1+8*ndt+m)=ndt+n ! row index
		       col(m)=m4
		       d(ndt+n)=0.0
		    endif

		    if (i.eq.nx-2) then ! At the right boundary
		       m5=(k-2)*nxy2+(j+1-2)*nx2+(nx-1)-1
		       m5=ndexfx(m5)

		       if(diff_node(m2).eq.1.and.diff_node(m5).eq.1) then
			  n=n+1
			  m=m+1
			  rw(8*ndt+m)=weight2*1.0 ! in Y direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m2
			  d(ndt+n)=0.0
		       
			  m=m+1
			  rw(8*ndt+m)=weight2*(-1.0) ! in Y direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m5
			  d(ndt+n)=0.0
		       endif
		    endif
		    
		    if (k.eq.nz-2) then ! At the bottom boundary
		       m1=(nz-1-2)*nxy2+(j-2)*nx2+i-1
		       m2=(nz-1-2)*nxy2+(j-2)*nx2+i+1-1
		       m3=(nz-1-2)*nxy2+(j+1-2)*nx2+i-1
		       m1=ndexfx(m1)
		       m2=ndexfx(m2)
		       m3=ndexfx(m3)
		       
		       if(diff_node(m1).eq.1.and.diff_node(m2).eq.1) then
			  n=n+1
			  m=m+1
			  rw(8*ndt+m)=weight1*1.0 ! in X direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m1
		       
			  m=m+1
			  rw(8*ndt+m)=weight1*(-1.0) ! in X direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m2
			  d(ndt+n)=0.0
		       endif
		       
		       if(diff_node(m1).eq.1.and.diff_node(m3).eq.1) then
			  n=n+1
			  m=m+1
			  rw(8*ndt+m)=weight2*1.0 ! in Y direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m1
		       
			  m=m+1
			  rw(8*ndt+m)=weight2*(-1.0) ! in Y direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m3
			  d(ndt+n)=0.0
		       endif
		       
		       if (i.eq.nx-2) then
			  m5=(nz-1-2)*nxy2+(j+1-2)*nx2+(nx-1)-1
			  m5=ndexfx(m5)
			
			  if(diff_node(m2).eq.1.and.diff_node(m5).eq.1) then
			     n=n+1
			     m=m+1
			     rw(8*ndt+m)=weight2*1.0 ! in Y direction
			     iw(1+8*ndt+m)=ndt+n ! row index
			     col(m)=m2
			     d(ndt+n)=0.0
			  
			     m=m+1
			     rw(8*ndt+m)=weight2*(-1.0) ! in Y direction
			     iw(1+8*ndt+m)=ndt+n ! row index
			     col(m)=m5
			     d(ndt+n)=0.0
			  endif
		       endif 
		    endif                                
		 enddo
	      enddo
	   enddo
	   
	   if (iuses.eq.2) then
	      do k=nz+2, 2*nz-2	! S velocity node
		 do j=2, ny-2
		    do i=2, nx-2
		       
		       m1=(k-4)*nxy2+(j-2)*nx2+i-1
		       m2=(k-4)*nxy2+(j-2)*nx2+i+1-1
		       m3=(k-4)*nxy2+(j+1-2)*nx2+i-1
		       m4=(k+1-4)*nxy2+(j-2)*nx2+i-1
		       m1=ndexfx(m1)
		       m2=ndexfx(m2)
		       m3=ndexfx(m3)
		       m4=ndexfx(m4)
		       dx=abs(xn(i)-xn(i+1))
		       dy=abs(yn(j)-yn(j+1))
		       dz=abs(zn(k-nz)-zn(k+1-nz))	
		       
		       if(diff_node(m1).eq.1.and.diff_node(m2).eq.1) then
			  m=m+1
			  n=n+1
			  rw(8*ndt+m)=weight1*1.0 ! in X direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m1
		       
			  m=m+1
			  rw(8*ndt+m)=weight1*(-1.0) ! in X direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m2
			  d(ndt+n)=0.0
		       endif
		       
		       if(diff_node(m1).eq.1.and.diff_node(m3).eq.1) then
			  m=m+1
			  n=n+1
			  rw(8*ndt+m)=weight2*1.0 ! in Y direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m1
	       
			  m=m+1
			  rw(8*ndt+m)=weight2*(-1.0) ! in Y direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m3
			  d(ndt+n)=0.0
		       endif
		       
		       if(diff_node(m1).eq.1.and.diff_node(m4).eq.1) then
			  m=m+1
			  n=n+1
			  rw(8*ndt+m)=weight3*1.0 ! in Z direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m1
			  
			  m=m+1
			  rw(8*ndt+m)=weight3*(-1.0) ! in Z direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m4
			  d(ndt+n)=0.0
		       endif

		       if (i.eq.nx-2) then ! At the right boundary
			  m5=(k-4)*nxy2+(j+1-2)*nx2+(nx-1)-1
			  m5=ndexfx(m5)
			  
			  n=n+1
			  m=m+1
			  rw(8*ndt+m)=weight2*1.0 ! in Y direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m2
			  d(ndt+n)=0.0
			  
			  m=m+1
			  rw(8*ndt+m)=weight2*(-1.0) ! in Y direction
			  iw(1+8*ndt+m)=ndt+n ! row index
			  col(m)=m5
			  d(ndt+n)=0.0
		       endif
		       
		       if (k.eq.2*nz-2) then ! At the bottom boundary
			  m1=(2*nz-1-4)*nxy2+(j-2)*nx2+i-1
			  m2=(2*nz-1-4)*nxy2+(j-2)*nx2+i+1-1
			  m3=(2*nz-1-4)*nxy2+(j+1-2)*nx2+i-1
			  m1=ndexfx(m1)
			  m2=ndexfx(m2)
			  m3=ndexfx(m3)
			  
			  if(diff_node(m1).eq.1.and.diff_node(m2).eq.1) then
			     n=n+1
			     m=m+1
			     rw(8*ndt+m)=weight1*1.0 ! in X direction
			     iw(1+8*ndt+m)=ndt+n ! row index
			     col(m)=m1
			  
			     m=m+1
			     rw(8*ndt+m)=weight1*(-1.0) ! in X direction
			     iw(1+8*ndt+m)=ndt+n ! row index
			     col(m)=m2
			     d(ndt+n)=0.0
			  endif
			  
			  if(diff_node(m1).eq.1.and.diff_node(m3).eq.1) then
			     n=n+1
			     m=m+1
			     rw(8*ndt+m)=weight2*1.0 ! in Y direction
			     iw(1+8*ndt+m)=ndt+n ! row index
			     col(m)=m1
			  
			     m=m+1
			     rw(8*ndt+m)=weight2*(-1.0) ! in Y direction
			     iw(1+8*ndt+m)=ndt+n ! row index
			     col(m)=m3
			     d(ndt+n)=0.0
			  endif
			  if (i.eq.nx-2) then
			     m5=(2*nz-1-4)*nxy2+(j+1-2)*nx2+(nx-1)-1
			     m5=ndexfx(m5)
			     
			     if(diff_node(m2).eq.1.and.diff_node(m5).eq.1) then
				n=n+1
				m=m+1
				rw(8*ndt+m)=weight2*1.0 ! in Y direction
				iw(1+8*ndt+m)=ndt+n ! row index
				col(m)=m2
				d(ndt+n)=0.0

				m=m+1
				rw(8*ndt+m)=weight2*(-1.0) ! in Y direction
				iw(1+8*ndt+m)=ndt+n ! row index
				col(m)=m5
				d(ndt+n)=0.0
			     endif
			  endif 
		       endif                               
		    enddo
		 enddo
	      enddo
	   endif 
	endif

	nndt=ndt+n		! # of rows
	total=m			! total number of nonzero model derivatives in matrix G
	write(*,*)'Total number of nonzero model derivatives after smoothing', total

	nar=8*ndt+total
	iw(1)=nar
c       write(*,*)'setting up the column indexes'
	do i=1, ndt
c       Set up column indexes with non-zero elements
	   iw(1+nar+      i) = 4*dt_ic1(i) - 3
	   iw(1+nar+  ndt+i) = 4*dt_ic1(i) - 2
	   iw(1+nar+2*ndt+i) = 4*dt_ic1(i) - 1
	   iw(1+nar+3*ndt+i) = 4*dt_ic1(i)
	   iw(1+nar+4*ndt+i) = 4*dt_ic2(i) - 3
	   iw(1+nar+5*ndt+i) = 4*dt_ic2(i) - 2
	   iw(1+nar+6*ndt+i) = 4*dt_ic2(i) - 1
	   iw(1+nar+7*ndt+i) = 4*dt_ic2(i)
	enddo

c	write(*,*)'Ending/....'
	do m=1,total
	   iw(1+nar+8*ndt+m)=col(m)+nev*4
	enddo
	
c       Scale G matrix so the L2 norm of each column is 1.
	write(log,'("~ scaling G columns ... ")')
c       write(*,*)'Scaling...'
c       G array scaling
	do i=1,4*nev+npari
	   norm(i) = 0.0
	enddo
	do i=1,nar
	   norm(iw(1+nar+i)) = norm(iw(1+nar+i)) + rw(i)**2
	enddo
c       write(16,*)'output the square of G columns'
	do i=1,nev*4+npari
c       write(16,*)i,norm(i)
	   norm(i) = sqrt(norm(i)/nndt)
	enddo
	
cz keep norm values for those effective nodes

	do i=1,new_loc+new_npari
	   j=old_index(i)
	   norm1(i)=norm(j)
	enddo

	do i=1,new_npari+new_loc
	   norm(i)=norm1(i)
	enddo

	do i=1, nar
	   i1=iw(1+nar+i)
	   if(i1.le.4*nev) then
	      if(fix(i1).eq.1) then ! invert location or origin time
		 iw(1+nar+i)=new_index(i1)
	      else
		 iw(1+nar+i)=new_npari+new_loc
		 rw(i)=0.0
	      endif
	   else	      
	      if(diff_node(i1-4*nev).eq.1) then
		 iw(1+nar+i)=new_index(i1)
	      else
c          ! set index larger than the effective range?
		 rw(i)= 0.0
		 iw(1+nar+i)=new_npari+new_loc ! 
	      endif      
	   endif
	enddo
	
 1	write(16,*)'new_loc+new_npari:',new_loc+new_npari
	do i=1,nar
	   if( iw(1+nar+i).ge.1.and.iw(1+nar+i).le.new_loc+new_npari) then
	      rw(i) = rw(i) / norm(iw(1+nar+i))
	   else
	      write(16,*)iw(1+nar+i), 'Column index is not in the range'
	   endif
	enddo
	
c       Testing...
	do i=1,new_loc+new_npari
	   norm_test(i) = 0.0
	enddo
	do i=1,nar
	   if(iw(1+nar+i).ge.1.and.iw(1+nar+i).le.new_loc+new_npari)
     &   	   norm_test(iw(1+nar+i)) = norm_test(iw(1+nar+i)) + rw(i)**2
	enddo
	do i=1,new_loc+new_npari
	   norm_test(i) = sqrt(norm_test(i)/nndt)
	   if (abs(norm_test(i)-1).gt.0.001) then
	      write(*,'("FATAL ERROR (lsqr: G scaling).")')
	      stop
	   endif
	enddo

c       Least square fitting using the algorithm of
c       Paige and Saunders, acm-trans. math. software, vol.8, no. 2,
c       jun., 1982, p. 195.

c       Set up input parameter first
	m = nndt
	n = new_loc+new_npari
	leniw = 2*nar+1
	lenrw = nar
	do i= 1,n
	   w1(i) = 0.0
	   w2(i) = 0.0
	   x(i) = 0.0
	   se(i) = 0.0
	enddo
	atol = 0.000001
	btol = 0.000001
c       conlim = 100000.0
	conlim=1200.0
c       itnlim = 100*n
	itnlim = 4*n
	istop = 0
	anorm = 0.0
	acond = 0.0
	rnorm = 0.0
	arnorm = 0.0
	xnorm = 0.0
	
	call datetime(dattim)
	write(log,'("~ lsqr ...    ", a)') dattim
	write(*,'("~ lsqr ...    ", a)') dattim

c d= data vector; w1,w2 = workspace; x= solution vector; se=solution error
      call lsqr(m, n, damp, 
     & leniw, lenrw, iw, rw, 
     & d, w1, w2, x, se, 
     & atol, btol, conlim, itnlim, 
     & istop, anorm, acond, rnorm, arnorm, xnorm)

      write(log,'("  istop = ",i1,"; acond (CND)=",f8.1,"; anorm =",f8.1,
     & "; arnorm =",f8.1,"; xnorm =",f8.1)')
     & istop, acond, anorm, arnorm, xnorm

      if (nsrc.eq.1) nsrc = nev

c     Rescale model vector
      do i=1,new_loc+new_npari
         x(i) = x(i) / norm(i)
         se(i) = se(i) / norm(i)
      enddo

      xnorm=0.0
      do i=1,new_loc+new_npari
         xnorm=xnorm+x(i)*x(i)
      enddo
      xnorm_vel=0.0
      do i=new_loc+1,new_loc+new_npari
         xnorm_vel=xnorm_vel+x(i)*x(i) ! hguo
      enddo
c     Unweight and rescale G matrix

c      do i=1,ndt
c         rw(i)       = rw(i)       * wtinv(i) * norm(iw(1+nar+i))
c         rw(ndt+i)   = rw(ndt+i)   * wtinv(i) * norm(iw(1+nar+ndt+i))
c         rw(2*ndt+i) = rw(2*ndt+i) * wtinv(i) * norm(iw(1+nar+2*ndt+i))
c         rw(3*ndt+i) = rw(3*ndt+i) * wtinv(i) * norm(iw(1+nar+3*ndt+i))
c         rw(4*ndt+i) = rw(4*ndt+i) * wtinv(i) * norm(iw(1+nar+4*ndt+i))
c         rw(5*ndt+i) = rw(5*ndt+i) * wtinv(i) * norm(iw(1+nar+5*ndt+i))
c         rw(6*ndt+i) = rw(6*ndt+i) * wtinv(i) * norm(iw(1+nar+6*ndt+i))
c       rw(7*ndt+i) = rw(7*ndt+i) * wtinv(i) * norm(iw(1+nar+7*ndt+i))	  
c	enddo

CZ    Unweight and rescale the model derivative matrix
      do i=1,nar
	 i1=iw(1+i)
	 i2=iw(1+nar+i)
         rw(i) = rw(i)* wtinv(i1) * norm(i2)
      enddo	

c     Compute residuals from d = G*x
      do i=1,ndt
         d(i) = -dt_res(i)*1000.0
      enddo
      call aprod(1, m, n, x, d, leniw, lenrw, iw, rw)
      do i=1,ndt
         dt_res(i) = -d(i)/1000.0
      enddo

       rnorm=0.0
       rnorm_wt=0.0
       do i=1,ndt
          rnorm=rnorm+d(i)*d(i)
          rnorm_wt=rnorm_wt+d(i)*d(i)*wt(i)*wt(i) ! hguo
       enddo
       write(*,*)'smooth, damp, xnorm, xnorm_vel, rnorm, rnorm_wt:',weight1,damp,sqrt(xnorm),sqrt(xnorm_vel),sqrt(rnorm),sqrt(rnorm_wt)
       write(log,*)'smooth, damp, xnorm, xnorm_vel, rnorm, rnorm_wt:',weight1,damp,sqrt(xnorm),sqrt(xnorm_vel),sqrt(rnorm),sqrt(rnorm_wt)

c     Get residual statistics (avrg, rms, var..)
      unknowns=new_loc+new_npari
      call resstat_FDD(log, idata, ndt, unknowns, dt_res, wt, dt_idx, 
     & rms_cc, rms_ct, rms_cc0, rms_ct0, 
     & rms_ccold, rms_ctold, rms_cc0old, rms_ct0old, 
     &             resvar1)


c     Scale errors
c The standard error estimates returned by LSQR increase monotonically
c with the iterations.  If LSQR shuts down early because of loose tolerances,
c or because the rhs-vector is special, the estimates will be too small.
c (I think they are most likely to be accurate if the rhs is random.)
c
c Remember that se(j) is covariance(j) / (m - n)
c where m - n = 1000000.  I've never quite understood why we
c divide by that number.

c Errors for the 95% confidence level,
c thus multiply the standard errors by 2.7955
      factor = 2.7955

cz Initialize solution 
      do i=1,nev
	 src_dt(i)=0.0
	 src_dx(i)=0.0
	 src_dy(i)=0.0
	 src_dz(i)=0.0
	 src_et(i)=0.0
	 src_ex(i)=0.0
	 src_ey(i)=0.0
	 src_ez(i)=0.0
      enddo

c     Store solution and errors
      do j=1,new_loc
	 i=old_index(j)
	 i1=int(i/4)
	 i2=mod(i,4)
	 if(i2.eq.0) then !origin time
	    src_dt(i1) = -x(j)
	    src_et(i1) = sqrt(se(j)) * sqrt(resvar1) *factor
	 elseif(i2.eq.1) then ! x-coordinate
	    src_dx(i1+1) = -x(j)
	    src_ex(i1+1) = sqrt(se(j)) * sqrt(resvar1) *factor
	 elseif(i2.eq.2) then ! y-coordinate
	    src_dy(i1+1) = -x(j)
	    src_ey(i1+1) = sqrt(se(j)) * sqrt(resvar1) *factor
	 else
	    src_dz(i1+1) = -x(j) ! z-coordinate
	    src_ez(i1+1) = sqrt(se(j)) * sqrt(resvar1) *factor
	 endif
      enddo
	   
      do i=1,nev
         src_cusp(i) = ev_cusp(i)
      enddo

c     initialize two vectors
      do i=1,npari
	 dv(i)=0.0
	 src_em(i)=0.0
      enddo
      do i=new_loc+1,new_npari+new_loc
	 dv(old_index(i)-4*nev) = -x(i)
c        Take weighted variance
	 atemp=sqrt(se(i))*sqrt(resvar1)*factor
	 src_em(old_index(i)-4*nev)=atemp
c         write(*,*) old_index(i), src_em(old_index(i))
      enddo

c     Get average errors and vector changes
      exavold = exav
      eyavold = eyav
      ezavold = ezav
      etavold = etav
      emavold = emav
      dxavold = dxav
      dyavold = dyav
      dzavold = dzav
      dtavold = dtav
      dmavold = dmav
       
      exav = 0.0
      eyav = 0.0
      ezav = 0.0
      etav = 0.0
      emav = 0.0
      dxav = 0.0
      dyav = 0.0
      dzav = 0.0
      dmav = 0.0
      do i=1,nev
         exav = exav + src_ex(i)
         eyav = eyav + src_ey(i)
         ezav = ezav + src_ez(i)
         etav = etav + src_et(i)
         dxav = dxav + abs(src_dx(i))
         dyav = dyav + abs(src_dy(i))
         dzav = dzav + abs(src_dz(i))
         dtav = dtav + abs(src_dt(i))
      enddo
      exav = exav/nev
      eyav = eyav/nev
      ezav = ezav/nev
      etav = etav/nev
      dxav = dxav/nev
      dyav = dyav/nev
      dzav = dzav/nev
      dtav = dtav/nev
      do i=1, new_npari
	emav=emav+src_em(i)
	dmav=dmav+abs(src_em(i))
      enddo
      emav=emav/new_npari
      dmav=dmav/new_npari

      if (iter.eq.1) then
         exavold = exav
         eyavold = eyav
         ezavold = ezav
         etavold = etav
	 emavold = emav
         dxavold = dxav
         dyavold = dyav
         dzavold = dzav
         dtavold = dtav
         dmavold = dmav
      endif

c     Output location statistics
      write(log,'(/,"Location summary:")')
      write(log,'(
     & " mean 2sig-error (x,y,z,t,m) [m,ms]: ",/,f7.1,f7.1,f7.1,
     & f7.1,f7.3,
     & " (",f7.1,f7.1,f7.1,f7.1,f7.3")",/,
     & " mean shift (x,y,z,t,m) [m,ms] (DX,DY,DZ,DT): ",/,
     & f7.1,f7.1,f7.1,f7.1,f7.3," (",f7.1,f7.1,f7.1,f7.1,f7.3")")')
     & exav, eyav, ezav, etav, emav,exav-exavold, eyav-eyavold, 
     & ezav-ezavold, etav-etavold, emav-emavold,
     & dxav, dyav, dzav, dtav, dmav,dxav-dxavold, dyav-dyavold, 
     & dzav-dzavold, dtav-dtavold,dmav-dmavold

c--- End of the subroutine lsfitHVFDD_lsqr_shot.f

              elseif (joint .eq.0) then ! only relocate events

                 if(iter.eq.1) then
                    resvar1= -999
                    unknowns=4*nev
                    call resstat_FDD(log,idata,ndt,unknowns,dt_res,dt_wt,dt_idx,
     &                   rms_cc,rms_ct,rms_cc0,rms_ct0,
     &                   rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
     &                   resvar1)
                 endif

                 call lsfitH_FDD_shot(log, iter, ndt, nev, nsrc,
     &                damp, eve_sta,
     &                idata, ev_cusp, src_cusp,
     &                dt_res, dt_wt,
     &                dt_ista, dt_ic1, dt_ic2,src_type,
     &                src_dx, src_dy, src_dz, src_dt, src_ex, src_ey, src_ez, src_et,
     &                exav, eyav, ezav, etav, dxav, dyav, dzav, dtav,
     &                rms_cc, rms_ct, rms_cc0, rms_ct0,
     &                rms_ccold, rms_ctold, rms_cc0old, rms_ct0old,
     &                tmp_xp, tmp_yp, tmp_zp, tmp_xs, tmp_ys, tmp_zs,dt_idx, acond) 

c              else              ! only velocity inversion
c                 write(*,*) 'Only velocity inversion....'
c                 if(iter.eq.1) then
c                    resvar1= -999
c                    unknowns=npari
c                    call resstat_FDD(log,idata,ndt,unknowns,dt_res,dt_wt,dt_idx,
c     &                   rms_cc,rms_ct,rms_cc0,rms_ct0,
c     &                   rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
c     &                   resvar1)
c                 endif
c                 call lsfitVFDD_lsqr(log,iter,ndt,nev,nsrc,damp,
c     &                idata,ev_cusp,src_cusp,
c     &                dt_res,dt_wt,
c     &                dt_ista,dt_ic1,dt_ic2, !new
c     &                src_dx,src_dy,src_dz,src_dt,src_ex,src_ey,src_ez,src_et,
c     &                exav,eyav,ezav,etav,emav,dxav,dyav,dzav,dtav,dmav,
c     &                rms_cc,rms_ct,rms_cc0,rms_ct0,
c     &                rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
c     &                tmp_xp,tmp_yp,tmp_zp,tmp_xs,tmp_ys,tmp_zs,
c     &                tmp_vp_index, tmp_vp,tmp_vs_index, tmp_vs,
c     &                dt_idx,acond, dv, 
c     &                weight1, weight2, weight3, norm_threshold)           
              endif    
           else
              write(*,*)'Singular Value decomposition is not used!'
              stop
           endif
           
c--- check for air quakes:
           mbad= 0
           k= 1
           do i= 1,nsrc
              call car2sph_ft((src_x(i)+src_dx(i))/1000, (src_y(i)+src_dy(i))/1000, (src_z(i)+src_dz(i))/1000,
     &             wlat,wlon,tlat,tlon,tdep, theta)
c              if(src_dep(i) + (src_dz(i)/1000).lt.0) then
              if(tdep.lt.air_dep) then
                 write(log,'(">>>Warning: negative depth - ",i12)')ev_cusp(i)
                 amcusp(k)= ev_cusp(i)
                 k=k+1
                 if(k.gt.1000) stop'>>> More than 1000 air quakes. Too many!'
              endif
           enddo
           mbad= k-1            ! number of neg depth events
           
c     update iteration numbers:
           if(mbad.gt.0) then
              do i= 1,niter
                 aiter(i)= aiter(i)+1
              enddo
              jiter= jiter+1    ! iteration with no update
              maxiter= maxiter+1
              
              write(log,*)'Number of air quakes (AQ) =',mbad
              if(nsrc-mbad .le. 1) then
                 write(*,*)'Warning: number of non-airquakes < 2'
                 write(*,*)'   skipping this cluster'
                 write(log,*)'Warning: number of non-airquakes < 2'
                 write(log,*)'   skipping this cluster'
                 goto 778
              endif       
              goto 500          ! skip the updating step
           endif
           
           
c---  update source parameters:
           xav= 0               ! mean centroid shift
           yav= 0
           zav= 0
           tav= 0
           alon= 0
           alat= 0
           adep= 0
           if(nsrc.eq.1) nsrc= nev
           do i= 1,nsrc
              src_cusp(i)= ev_cusp(i)
c     update absolute source parameters (cart)

              if(abs(src_dx(i)).gt.2000) then
                 src_dx(i)=2000*src_dx(i)/abs(src_dx(i))
              endif

              if(abs(src_dy(i)).gt.2000) then
                 src_dy(i)=2000*src_dy(i)/abs(src_dy(i))
              endif

              if(abs(src_dz(i)).gt.3000) then
                 src_dz(i)=3000*src_dz(i)/abs(src_dz(i))
              endif

              src_x(i)= src_x(i) + src_dx(i)
              src_y(i)= src_y(i) + src_dy(i)
              src_z(i)= src_z(i) + src_dz(i)
              src_t(i)= src_t(i) + src_dt(i)
              
c     update absolute source locations (geogr)
c     src_dep(i)= src_dep(i) + (src_dz(i)/1000)
c     call SDC2(src_x(i)/1000,src_y(i)/1000,lat,lon,1)
              call car2sph_ft(src_x(i)/1000, src_y(i)/1000, src_z(i)/1000,
     &             wlat,wlon,tlat,tlon,tdep, theta)
              src_lon(i)= tlon
              src_lat(i)= tlat
              src_dep(i)= tdep
c     alon= lon+alon	
c     alat= lat+alat
              alon= lon+src_lon(i)	
              alat= lat+src_lat(i)
              adep= adep+src_dep(i)                 
c     get mean centroid shift
              xav= xav + (src_x(i) - src_x0(i))	
              yav= yav + (src_y(i) - src_y0(i))
              zav= zav + (src_z(i) - src_z0(i))
              tav= tav + (src_t(i) - src_t0(i))
           enddo
           xav= xav/nsrc
           yav= yav/nsrc
           zav= zav/nsrc
           tav= tav/nsrc
           alon= alon/nsrc
           alat= alat/nsrc
           adep= adep/nsrc
           
c     update velocities
c     write(*,*)'In main routines...'
           if ( joint .ne. 0) then
              write(*,*)'updating the velocities'
              do n=1, npari
                 nn=mdexfx(n)         
c     calculate x and z indices of velocity grid
                 k=(nn-1)/nxy2+2
                 j=2+(nn-1+(2-k)*nxy2)/nx2
                 i=1+nn+nx2*(2-j)+nxy2*(2-k)   
                 if (k.ge.nz) k=k+2 ! if s velocity node 
                 slow=1.0/vel(i,j,k)+dv(n)/1000.0 !slowness
                 dvel=1.0/slow-vel(i,j,k)                 
                 if (k.le.nz) then !P-wave velocity
                    if(abs(dvel).lt.0.4) then
                       vel(i,j,k)=1.0/slow
                    else
                       vel(i,j,k)=vel(i,j,k)+dvel/abs(dvel)*0.4
                    endif
                    if(vel(i,j,k).lt.1.0) vel(i,j,k)=1.0
                    if(vel(i,j,k).gt.8.7) vel(i,j,k)=8.7
                 endif
                 if (k.gt.nz) then
                    if(abs(dvel).lt.0.2) then
                       vel(i,j,k)=1.0/slow
                    else
                       vel(i,j,k)=vel(i,j,k)+dvel/abs(dvel)*0.2
                    endif                     
                    if(vel(i,j,k).lt.0.8) vel(i,j,k)=0.8
                    if(vel(i,j,k).gt.4.9) vel(i,j,k)=4.9
                 endif                
              enddo
c     output velocities
              write(16,*)'P-wave velocity at iteration:',iter
              do k=1,nz
                 do j=1,ny
                    write(16, '(100f7.3)')(vel(i,j,k),i=1,nx)
                 enddo
              enddo
              write(16,*)'End of the velocity'

              if(iuses.eq.2) then
                 write(16,*)'S-wave velocity at iteration:',iter
                 do k=nz+1,2*nz
                    do j=1,ny
                       write(16, '(100f7.3)')(vel(i,j,k),i=1,nx)
                    enddo
                 enddo                 
              endif

           endif
           
           write(log,'("  cluster centroid at:",1x,f10.6,2x,f11.6,2x,f9.6)')
     &          alat,alon,adep
           write(log,'("  mean centroid (origin) shift in x,y,z,t [m,ms]: ",/
     &          f7.1,f7.1,f7.1,f7.1)'),xav,yav,zav,tav
           write(log,'("  (OS in std output gives maximum value.)")')
           
c---  get interevent distance for each observation and average signal coherency:
           cohav= 0
           picav= 0
           j= nct
           k= ncc
           ncc= 0
           nct= 0
           do i= 1,ndt
              dt_offs(i)= sqrt((src_x(dt_ic1(i))-src_x(dt_ic2(i)))**2 +
     &             (src_y(dt_ic1(i))-src_y(dt_ic2(i)))**2 +
     &             (src_z(dt_ic1(i))-src_z(dt_ic2(i)))**2)
              
              if(dt_idx(i).le.2) then
                 cohav= cohav + sqrt(dt_qual(i))
                 ncc= ncc+1
              else
                 picav= picav + dt_qual(i)
                 nct= nct+1
              endif
              
           enddo
           cohav= cohav/ncc
           picav= picav/nct
           write(log,'(/,"More:")')
           write(log,'("  mean phase coherency = ",f5.3)')cohav
           write(log,'("  mean pick quality = ",f5.3)')picav
           
c---  get number of observations and mean residual at each station
           tmpr1= 0
           tmpr2= 0
           do i= 1,nsta
              sta_np(i)= 0
              sta_ns(i)= 0
              sta_nnp(i)= 0
              sta_nns(i)= 0
              sta_rmsc(i)= 0
              sta_rmsn(i)= 0
              do j= 1,ndt
                 if(i.eq.dt_ista(j)) then
                    if(dt_idx(j).le.2) then
                       sta_rmsc(i)= sta_rmsc(i)+dt_res(j)**2
                       if(dt_idx(j).eq.1) then
                          sta_np(i)= sta_np(i)+1
                       else
                          sta_ns(i)= sta_ns(i)+1
                       endif
                    else
                       sta_rmsn(i)= sta_rmsn(i)+dt_res(j)**2
                       if(dt_idx(j).eq.3) then
                          sta_nnp(i)= sta_nnp(i)+1
                       else
                          sta_nns(i)= sta_nns(i)+1
                       endif
                    endif
                 endif
              enddo
              
              if(sta_np(i)+sta_ns(i).gt.0)
     &             sta_rmsc(i)= sqrt(sta_rmsc(i)/(sta_np(i)+sta_ns(i)))
              if(sta_nnp(i)+sta_nns(i).gt.0)
     &             sta_rmsn(i)= sqrt(sta_rmsn(i)/(sta_nnp(i)+sta_nns(i)))
              if(sta_rmsc(i).gt.tmpr1) then
                 tmpr1= sta_rmsc(i)
                 k= i
              endif
              if(sta_rmsn(i).gt.tmpr2) then
                 tmpr2= sta_rmsn(i)
                 l= i
              endif
           enddo
           tmpr1= tmpr1*1000
           tmpr2= tmpr2*1000
           if(idata.eq.1.or.idata.eq.3) then
              write(log,'("  station with largest cc rms: ",a7,"=",
     &             f7.0," ms (RMSST)")')
     &             sta_lab(k),tmpr1
           endif
           if(idata.eq.2.or.idata.eq.3) then
              write(log,'("  station with largest ct rms: ",a7,"=",
     &             f7.0," ms (RMSST)")')
     &             sta_lab(l),tmpr2
           endif
           
c---  write output scratch mdat.reloc:
           n= trimlen(fn_reloc)
           i=iter-jiter
           write(str80,'(a,".",i3.3,".",i3.3)')fn_reloc(1:n),iclust,i
           call freeunit(iunit)
           open(iunit,file=str80,status='unknown')
           write(iunit,'(i9,1x,f10.6,1x,f11.6,1x,f9.3,1x,f10.1,1x,f10.1,
     &          1x,f10.1,
     &          1x,f8.1,1x,f8.1,1x,f8.1,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3,
     &          1x,f3.1,1x,i3)')
     &          (src_cusp(i),src_lat(i),src_lon(i),src_dep(i),src_x(i),src_y(i),
     &          src_z(i),src_ex(i),src_ey(i),src_ez(i),int(ev_date(i)/10000),
     &          int(mod(ev_date(i),10000)/100),mod(ev_date(i),100),
     &          int(ev_time(i)/1000000),int(mod(ev_time(i),1000000)/10000),
     &          mod(real(ev_time(i)),10000.0)/100,ev_mag(i),iclust,i=1,nev)
           close(iunit)
c     WARNING: variable "str" is set to zero value by default
           write(log,'(/,"Relocation results for this iteration are"
     &          " stored in ",a)')str80(1:trimlen(str80))
           
 500       continue             ! case of air quakes
           
c     standard output:
           if(mbad.gt.0) then
              str3='   '
           else
              n= iter-jiter
              if(n.lt.1000) write(str3,'(i3)')n
              if(n.lt.100) write(str3,'(1x,i2)')n
              if(n.lt.10) write(str3,'(2x,i1)')n
           endif
           if(isolv.eq.1.and.idata.eq.3) then
              if(iter.eq.1) write(*,'(/,"  IT   EV  CT  CC",
     &             "    RMSCT      RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ",/,
     &             "        %   %   %",
     &             "   ms     %   ms     %    ms    m    m    m   ms    m ")')
              write(*,'(i2,a3,3(1x,i3),i5,f6.1,i5,f6.1,i6,4i5,i5,i4)')
     &             iter,str3,
     &             nint(nev*100./nevold),nint(nct*100.0/nctold),
     &             nint(ncc*100.0/nccold),
     &             nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     &             nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     &             nint(max(tmpr1,tmpr2)),
     &             nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     &             nint(max(abs(xav),abs(yav),abs(zav))),mbad
           endif
           if(isolv.eq.1.and.idata.eq.1) then
              if(iter.eq.1) write(*,'(/,"  IT   EV  CC",
     &             "    RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ",/,
     &             "        %   %",
     &             "   ms     %    ms    m    m    m   ms    m ")')
              write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4)')
     &             iter,str3,
     &             nint(nev*100./nevold),
     &             nint(ncc*100.0/nccold),
     &             nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     &             nint(max(tmpr1,tmpr2)),
     &             nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     &             nint(max(abs(xav),abs(yav),abs(zav))),mbad
           endif
           if(isolv.eq.1.and.idata.eq.2) then
              if(iter.eq.1) write(*,'(/,"  IT   EV  CT",
     &             "    RMSCT     RST   DX   DY   DZ   DT   OS  AQ",/,
     &             "        %   %",
     &             "   ms     %    ms    m    m    m   ms    m ")')
              write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4)')
     &             iter,str3,
     &             nint(nev*100./nevold),
     &             nint(nct*100.0/nctold),
     &             nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     &             nint(max(tmpr1,tmpr2)),
     &             nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     &             nint(max(abs(xav),abs(yav),abs(zav))),mbad
           endif
           
           if(isolv.eq.2.and.idata.eq.3) then
              if(iter.eq.1) write(*,'(/,"  IT   EV  CT  CC",
     &             "    RMSCT      RMSCC   RMSST   DX   DY   DZ   DT   ",
     &             "OS  AQ  CND",/,
     &                           "        %   %   %",
     &             "   ms     %   ms     %    ms    m    m    m   ms   ",
     &             " m     ")')
              write(*,'(i2,a3,3(1x,i3),i5,f6.1,i5,f6.1,i6,4i5,i5,i4,i5)')
     &             iter,str3,
     &             nint(nev*100./nevold),nint(nct*100.0/nctold),
     &             nint(ncc*100.0/nccold),
     &             nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     &             nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     &             nint(max(tmpr1,tmpr2)),
     &             nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),mbad,nint(acond)
           endif
           if(isolv.eq.2.and.idata.eq.1) then
              if(iter.eq.1) write(*,'(/,"  IT   EV  CC",
     &             "    RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ  CND",/,
     &             "        %   %",
     &             "   ms     %    ms    m    m    m   ms    m ")')
              write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4,i5)')
     &             iter,str3,
     &             nint(nev*100./nevold),
     &             nint(ncc*100.0/nccold),
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     &             nint(max(tmpr1,tmpr2)),
     &             nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     &             nint(max(abs(xav),abs(yav),abs(zav))),mbad,nint(acond)
      endif
      if(isolv.eq.2.and.idata.eq.2) then
         if(iter.eq.1) write(*,'(/,"  IT   EV  CT",
     &        "    RMSCT   RMSST   DX   DY   DZ   DT   OS  AQ  CND",/,
     &                           "        %   %",
     &        "   ms     %    ms    m    m    m   ms    m ")')
         write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4,i5)')
     &        iter,str3,
     &        nint(nev*100./nevold),
     &        nint(nct*100.0/nctold),
     &        nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     &        nint(max(tmpr1,tmpr2)),
     &        nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     &        nint(max(abs(xav),abs(yav),abs(zav))),mbad,nint(acond)
      endif
      
      call datetime(dattim)
      write(log,'("Iteration ",i2," finished ",a)') iter, dattim
      
      if(iter.eq.maxiter) goto 600 ! all iterations done.
      iter= iter+1
      goto 55                   ! next iteration
      
c--- update origin time (this is only done for final output!!)
 600  continue
      write(*,'(/,"writing out results ...")')
      do i= 1,nev
         src_t(i)= src_t(i)/1000 !from here on src_t in sec!!
c         if(src_t(i).gt.5) then
c            write(*,'("Warning: (src_t>5.0). ")')
c         endif
         iyr= int(ev_date(i)/10000)
         imo= int(mod(ev_date(i),10000)/100)
         idy= int(mod(ev_date(i),100))
         ihr= int(ev_time(i)/1000000)
         imn= int(mod(ev_time(i),1000000)/10000)
         itf= JULIAM(iyr,imo,idy,ihr,imn)

c         sc= (mod(real(ev_time(i)),10000)/100) + src_t(i)
         sc= (mod(real(ev_time(i)),10000.0)/100) - src_t(i)
         itf= itf + int(sc / 60.)
         sc=  sc  - int(sc / 60.)*60.
         if(sc.lt.0) then
            itf= itf-1
            sc= 60. + sc
         endif
         call DATUM(itf,iyr,imo,idy,ihr,imn)
         ev_date(i)= iyr*10000 + imo*100 + idy
         ev_time(i)= ihr*1000000 + imn*10000 + nint(sc*100)
      enddo
      
c---  get # of obs per event:
      do i=1,nev
         src_np(i)= 0
         src_ns(i)= 0
         src_nnp(i)= 0
         src_nns(i)= 0
         src_rmsc(i)= 0
         src_rmsn(i)= 0
      enddo
      do i=1,ndt
         if(dt_idx(i).eq.1) then
            src_np(dt_ic1(i))= src_np(dt_ic1(i))+1
            src_np(dt_ic2(i))= src_np(dt_ic2(i))+1
         endif
         if(dt_idx(i).eq.2) then
            src_ns(dt_ic1(i))= src_ns(dt_ic1(i))+1
            src_ns(dt_ic2(i))= src_ns(dt_ic2(i))+1
         endif
         if(dt_idx(i).le.2) then
            src_rmsc(dt_ic1(i))= src_rmsc(dt_ic1(i))+dt_res(i)**2
            src_rmsc(dt_ic2(i))= src_rmsc(dt_ic2(i))+dt_res(i)**2
         endif
         if(dt_idx(i).eq.3) then
            src_nnp(dt_ic1(i))= src_nnp(dt_ic1(i))+1
            src_nnp(dt_ic2(i))= src_nnp(dt_ic2(i))+1
         endif
         if(dt_idx(i).eq.4) then
            src_nns(dt_ic1(i))= src_nns(dt_ic1(i))+1
            src_nns(dt_ic2(i))= src_nns(dt_ic2(i))+1
         endif
         if(dt_idx(i).ge.3) then
            src_rmsn(dt_ic1(i))= src_rmsn(dt_ic1(i))+dt_res(i)**2
            src_rmsn(dt_ic2(i))= src_rmsn(dt_ic2(i))+dt_res(i)**2
         endif
      enddo
      do i=1,nev
         src_rmsc(i)= sqrt(src_rmsc(i)/nev)
         src_rmsn(i)= sqrt(src_rmsn(i)/nev)
      enddo
      
c---  output final residuals: mdat.res
      if(trimlen(fn_res).gt.1) then
         call freeunit(iunit)
         open(iunit,file=fn_res,status='unknown')
         write(iunit,'("STA",11x,"DT",8x,
     &        "C1",8x,"C2",4x,"IDX",5x,"QUAL",4x,"RES [ms]",3x,"WT",9x,
     &        "OFFS")')
         write(iunit,'(a7,1x,f12.7,1x,i9,1x,i9,1x,i1,1x,
     &        f9.4,1x,f12.6,1x,f11.6,1x,f8.1)')
     &        (dt_sta(j),dt_dt(j),dt_c1(j),dt_c2(j),dt_idx(j),dt_qual(j),
     &        dt_res(j)*1000,dt_wt(j),dt_offs(j),j=1,ndt)
         close(iunit)
      endif

c---  output final locations (mdat.reloc):
      write(fu1,'(i9,1x,f10.6,1x,f11.6,1x,f9.3,1x,f10.1,1x,f10.1,
     &     1x,f10.1,
     &     1x,f8.1,1x,f8.1,1x,f8.1,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3,
     &     1x,f3.1,1x,i5,1x,i5,1x,i5,1x,i5,1x,f6.3,1x,f6.3,1x,i3)')
     &     (src_cusp(i),src_lat(i),src_lon(i),src_dep(i),src_x(i),src_y(i),
     &     src_z(i),src_ex(i),src_ey(i),src_ez(i),int(ev_date(i)/10000),
     &     int(mod(ev_date(i),10000)/100),mod(ev_date(i),100),
     &     int(ev_time(i)/1000000),int(mod(ev_time(i),1000000)/10000),
     &     mod(real(ev_time(i)),10000.0)/100,ev_mag(i),
     &     src_np(i),src_ns(i),src_nnp(i),src_nns(i),
     &     src_rmsc(i),src_rmsn(i), iclust,i=1,nev)
      
c---  output stations (mdat.station):
      if(trimlen(fn_stares).gt.1) then
         write(fu3,'(a7,1x,f9.4,1x,f9.4,1x,f9.4,1x,f9.4,1x,i7,1x,
     &        i7,1x,i7,1x,i7,1x,f9.4,1x,f9.4,1x,i3)')
     &        (sta_lab(i),sta_lat(i),sta_lon(i),sta_dist(i),sta_az(i),
     &        sta_np(i),sta_ns(i),sta_nnp(i),sta_nns(i),
     &        sta_rmsc(i),sta_rmsn(i),iclust,i=1,nsta)
      endif
      
778   continue
      enddo                     ! loop over clusters (iclust)
      
      close(fu0)
      if(trimlen(fn_stares).gt.1) close(fu3)
      close(fu1)

c     output final velocity models

c      call freeunit(iunit0)
c      call freeunit(iunit1)
      iunit0=22
      iunit1=25
      open(iunit0,file=fn_vp,status='unknown')
      open(iunit1,file=fn_vs,status='unknown')
c     output the velocity model based on latitude, longitude 
c     and depth

c     !do tdep=dep0,dep1,5.0
c     do tdep=dep0,dep1,ddep
c         write(*,*)'tdep=',tdep
c        !do tlat=lat0,lat1,0.5
c        do tlat=lat0,lat1,dlat
c           i=0
c           !do tlon=lon0,lon1,0.5
c           do tlon=lon0,lon1,dlon
c              isp=0 ! P-wave velocity
c              call vel3(isp, tlon, tlat, tdep, v)
c              i=i+1
c              vp(i)=v	      
c              if(iuses.eq.2) then ! S-wave velocity
c                 isp=1
c                 call vel3(isp, tlon, tlat, tdep, v)
c                 vs(i)=v
c              endif 
c           enddo
c           write(iunit0, '(200f7.3)')(vp(j),j=1,i)
c           if(iuses.eq.2) write(iunit1, '(200f7.3)')(vs(j),j=1,i)
c       enddo
c     enddo

c     output velocities
      do k=1,nz
         do j=1,ny
            write(iunit0, '(100f7.3)')(vel(i,j,k),i=1,nx)
         enddo
      enddo

      if(iuses.eq.2) then 
         do k=nz+1,2*nz
            do j=1,ny
               write(iunit1, '(100f7.3)')(vel(i,j,k),i=1,nx)
            enddo
         enddo
      endif
      
      end                       !of main routine
