/*    punch   */ 
 
/*  from Vidale's slug3d 
    edited by j.hole 02-08-90    add ability to accept as 
       source some input values of travel times on the faces 
       of the data volume   
    edited by j.hole 06-09-90    fix bugs in Vidale's code 
       remove edges and corners, including 
       them in the sides calculation        add ability to choose 
       rate of speed of box expansion in each of 6 directions 
       **** RENAMED slugjah.c 
    edited by j.hole 17-09-90    replace source cube of constant 
       velocity with a source cube of linear velocity gradient 
    edited by j.hole 11-91       new source option: input a tt field 
       that has already been calculated, start calculations on any 
       wall of volume, change tt only if a smaller value is found 
       (allows rays to travel back into expanding cube) 
    MAJOR edits by j.hole 11/91  new fd stencils and stencil logic 
       to attempt to allow head waves to travel parallel the faces  
       of the expanding cube (and at some other angles too) 
       **** RENAMED bull3d.c 
    MAJOR edits by j.hole 01/93  extended bull3d.c logic to  
       systematically test all possible operators 
       **** RENAMED bully.c 
    MAJOR edits by j.hole 12/93  removed bug in bully.c;  only use 
       "new face" and "new edge" stencils on new faces and edges! 
       **** RENAMED punch.c 
    edited by j.hole 01/95  added ability to automatically detect  
       head waves travelling on the faces of the expanding cube  
       of already-calculated values; then restarts calculations  
       at appropriate wall of the model to allow the waves to  
       travel back into the (previous) expanding cube 
****Hole, J.A., and B.C. Zelt, 1995.  3-D finite-difference reflection  
       traveltimes.  Geophys. J. Int., 121, 427-434. 
       This reference describes the modifications for head waves  
       made 11-91, 01/93, 12/93, and 01/95.   
*/ 
 
 
/*  RECEIVED IN E-MAIL 12-02-90  */ 
/*  
Message inbox:3 -  Read 
From:    <vid@flatte.ucsc.edu> 
To:      <hole@geop.ubc.ca>, <vid@flatte.ucsc.edu> 
Cc:      <vid@rupture.ucsc.edu> 
*** NOTE THAT THESE ADDRESSES ARE OUTDATED *** 
Subject: Re:  3-D Travel Times 
 
I can send the code, but this week I am too busy to send sample input and  
output files. 
The code is included below. 
No special compiling options are necessary. 
Two input files are necessary. 
One, in ascii, has lines that look like: 
 
nx=100 
ny=100 
nz=100 
xs=10 
ys=10 
zs=10 
h=1 
timefile=run0.times 
velfile=input.vel 
 
(Notice that there are no spaces.) 
The other file is the velocity field (an nx by ny by nz array) that is 
written in C binary, and with the parameter file above would have the 
name "input.vel". 
When the program is finished, the traveltimes to all the grid points will 
appear in a file called run0.times. 
If the C code was compiled with the command "cc time3d.c -o time3d", and 
the parameter file is named "time.par", the 
program would run with the command "time3d par=time.par". 
If the program gives the correct answer for a uniform wholespace,  
it is probably working. 
Send me some mail if something doesn't work, and give me a call (408 459-4585) 
if I don't answer the e-mail (our computers are in a state of flux).   */ 
 
/* PROGRAM TO CALCULATE TRAVEL TIMES IN 3D MEDIA */ 
/* author: John E. Vidale, started Sept. 11, 1986*/ 
/*  restarted, 8-11-88, with more accurate scheme */ 
/*  last revised, 2-89, by Glenn Nelson for 5x5x5 source */ 
/* UC Santa Cruz, all rights reserved */ 
 
#include    <stdio.h> 
#include    <math.h> 
#include    <fcntl.h> 
 
#define PI 3.141592654 
#define SQR2 1.414213562 
#define SQR3 1.732050808 
#define SQR6 2.449489743 
/* #define	NCUBE	2 */	/* 1 for 3x3x3, 2 for 5x5x5	*/ 
#define t0(x,y,z)   time0[nxy*(z) + nx*(y) + (x)] 
#define s0(x,y,z)   slow0[nxy*(z) + nx*(y) + (x)] 
#define	SQR(x)	((x) * (x)) 
#define	DIST(x,y,z,x1,y1,z1)	sqrt(SQR(x-(x1))+SQR(y-(y1)) + SQR(z-(z1))) 
 
struct sorted 
	{ float time; int x1, x2;}; 
 
/* FUNCTION DECLARATIONS	*/ 
int 
	compar(); 
float 
	ex0(), fd5(), fd6(), fd7(), fd8(),         /* VIDALE'S STENCILS */ 
        fdh3d(),fdhne(),fdh2d(),fdhnf();           /* HOLE'S STENCILS */ 
/* Make sure vfield and tfield start at 0.... */ 
 
void punch_(int *nx_ft, int *ny_ft, int *nz_ft, float *h_ft, 
           int *srctype_ft, int* floatsrc_ft, float *fxs_ft, float *fys_ft,  
           float *fzs_ft, int *reverse_ft, int *ncube_ft, 
	   int *iminus_ft, int *iplus_ft, int *jminus_ft, int *jplus_ft, 
	   int *kminus_ft, int *kplus_ft, int *headpref_ft, float *maxoff_ft, 
           float *vfield, float *tfield) 
{ 
	/* NOTE THAT SEVERAL VARIABLES MUST BE SPECIFIED IN par=xxx FILE,  
	   WHILE OTHERS ARE OPTIONAL:  IF A mstpar STATEMENT READS THE  
	   VARIABLE BELOW, THEN THE VARIABLE IS REQUIRED;  IF A getpar  
	   STATEMENT IS USED BELOW, THEN THE VARIABLE IS OPTIONAL */ 
	int 
		nx,		/* x-dimension of mesh (LEFT-TO-RIGHT) */ 
		ny,		/* y-dimension of mesh (FRONT-TO-BACK) */ 
		nz,		/* z-dimension of mesh  (TOP-TO-BOTTOM) */ 
		iplus=1,	/* rate of expansion of "cube" in the */ 
		iminus=1,	/*    plus/minus x/y/z direction */ 
		jplus=1, 
		jminus=1, 
		kplus=1, 
		kminus=1, 
		igrow,		/* counter for "cube" growth */ 
		floatsrc=0,	/* IF 0, SOURCE IS ON A GRID POINT */ 
                srctype=1,      /* if 1, source is a point; 
                                      2, source is on the walls of the data volume; 
				      3, source on wall, time field known; */ 
                srcwall,        /* if 1, source on x=0 wall, if 2, on x=nx-1 wall 
                                   if 3, source on y=0 wall, if 4, on y=ny-1 wall 
                                   if 5, source on z=0 wall, if 6, on z=nz-1 wall */ 
		xs,		/* shot x position (in grid points) */ 
		ys,		/* shot y position */ 
		zs,		/* shot depth */ 
		xx, yy, zz,		/* Used to loop around xs, ys, zs coordinates	*/ 
		X1, X2, lasti, index, ii, i, j, k, radius,  
	        tfint, vfint, wfint, ofint, 
		nxy, nyz, nxz, nxyz, nwall, 
		/* counters for the position of the sides of current cube */ 
		x0, x1, y0, y1, z0, z1, 
		/* flags set to 1 until a side has been reached */ 
		dx0=1, dx1=1, dy0=1, dy1=1, dz0=1, dz1=1, rad0=1, 
                /* maximum radius to compute */ 
	        maxrad, 
		/* used in linear velocity gradient cube source */ 
		xxx,yyy,zzz, 
                /* size of source cube     1 for 3x3x3, 2 for 5x5x5...	*/ 
                NCUBE=2, 
		reverse=1,	/* will automatically do up to this number of 
			reverse propagation steps to fix waves that travel  
			back into expanding cube */ 
		headpref=6,	/* if headpref starts > 0, will determine  
			model wall closest to source and will prefer to start 
			reverse calculations on opposite wall */ 
		/* counters for detecting head waves on sides of current cube */ 
		head,headw[7]; 
	float 
		h,		/* spatial mesh interval (units consistant with vel) */ 
		fxs,	/* shot position in X (in real units)*/ 
		fys,	/* shot position in Y (in real units)*/ 
		fzs,	/* shot position in Z (in real units)*/ 
		*slow0, *time0, *wall, 
		s000, guess, try, slo, 
                /* maximum offset (real units) to compute */ 
	        maxoff = -1., 
		/* used in linear velocity gradient cube source */ 
		rx, ry, rz, dvx, dvy, dvz, dv, v0, 
		rzc, rxyc, rz1, rxy1, rho, theta1, theta2, 
		/* used to detect head waves:  if headwave operator decreases  
		   the previously-computed traveltime by at least  
		   headtest*<~time_across_cell> then the headwave counter is  
		   triggered */ 
		fhead,headtest=1.e-3, 
		/* pointers to times and slownesses */ 
		*r0, *r1, *r2, *r3, 
		*p0, *p1, *p2, *p3, *p4, *p5; 
	char 
		velfile[80],	/* file though which velocity structure is input */ 
		oldtfile[80],	/* file through which old travel times are input */ 
		timefile[80],	/* file in which travel times appear at the end */ 
                wallfile[80];   /* file containing input wall values of traveltimes */ 
 
	/* ARRAY TO ORDER SIDE FOR SOLUTION IN THE RIGHT ORDER */ 
	struct sorted *sort; 

	/* 
	fprintf(stderr,"Starting slug3d: by J. Vidale, 1988, UCSC\n"); 
	fprintf(stderr,"Starting punch: by J. Hole, 1993, UBC-Stanford\n"); 
	*/

	/* INITIALIZE PARAMETERS AND ARRAYS */ 
/*setpar(ac,av); */
/* 
	mstpar("nx",     "d", &nx); 
	mstpar("ny",     "d", &ny); 
	mstpar("nz",     "d", &nz); 
	mstpar("h",      "f", &h); 
	getpar("iminus","d",&iminus); 
	getpar("iplus","d",&iplus); 
	getpar("jminus","d",&jminus); 
	getpar("jplus","d",&jplus); 
	getpar("kminus","d",&kminus); 
	getpar("kplus","d",&kplus); 
	getpar("floatsrc","d",&floatsrc); 
	getpar("srctype","d",&srctype); 
	getpar("NCUBE","d",&NCUBE); 
	getpar("reverse","d",&reverse); 
	getpar("headpref","d",&headpref); 
	getpar("maxoff","f",&maxoff); 
*/ 
 
/* Add by H. Zhang to make punch.c is callable from Fortran */ 
	nx=*nx_ft; 
	ny=*ny_ft; 
	nz=*nz_ft; 
	h= *h_ft; 
	floatsrc=*floatsrc_ft; 
	fxs=*fxs_ft; 
	fys=*fys_ft; 
	fzs=*fzs_ft; 
	NCUBE=*ncube_ft; 
	reverse=*reverse_ft; 
	iminus=*iminus_ft; 
	iplus=*iplus_ft; 
	jminus=*jminus_ft; 
	jplus=*jplus_ft; 
	kminus=*kminus_ft; 
	kplus=*kplus_ft; 
	headpref=*headpref_ft; 
	maxoff=*maxoff_ft; 
	srctype=*srctype_ft; 
 
	 
	fxs/=h; 
	fys/=h; 
	fzs/=h;	 
	xs = (int)(fxs + 0.5); 
	ys = (int)(fys + 0.5); 
	zs = (int)(fzs + 0.5);	 
 
/* 
 
	if  (srctype==1)  { 
	   if(floatsrc==0){ 
		mstpar("xs","d",&xs); 
		mstpar("ys","d",&ys); 
		mstpar("zs","d",&zs); 
		fxs = (float)xs; 
		fys = (float)ys; 
		fzs = (float)zs; 
	   } 
	   else{ 
		mstpar("fxs","f",&fxs); 
		mstpar("fys","f",&fys); 
		mstpar("fzs","f",&fzs); 
		fxs /= h; 
		fys /= h; 
		fzs /= h; 
		xs = (int)(fxs + 0.5); 
		ys = (int)(fys + 0.5); 
		zs = (int)(fzs + 0.5); 
	   } 
	} 
	else if (srctype==2)  { 
	        mstpar("srcwall","d",&srcwall); 
		mstpar("wallfile","s",wallfile); 
	        } 
	else if (srctype==3)  { 
	        mstpar("srcwall","d",&srcwall); 
		mstpar("oldtfile","s",oldtfile); 
	        } 
	else  { 
		fprintf(stderr,"ERROR: incorrect value of srctype\n"); 
		exit(-1); 
	} 
 
	mstpar("timefile","s",timefile); 
	mstpar("velfile","s",velfile); 
	endpar(); 
*/ 

	/*
	fprintf(stderr,"Enter the time grid file name\n");
	scanf("%s", timefile);
	*/

	if(xs<2 || ys<2 || zs<2 || xs>nx-3 || ys>ny-3 || zs>nz-3){ 
		fprintf(stderr,"Source near an edge, beware of traveltime errors\n"); 
		fprintf(stderr,"for raypaths that travel parallel to edge \n"); 
		fprintf(stderr,"while wavefronts are strongly curved, (JV, 8/17/88)\n"); 
	} 
 
	/* SET MAXIMUM RADIUS TO COMPUTE */ 
	if (maxoff > 0.) { 
	  maxrad = maxoff/h + 1; 
	  fprintf(stderr,"WARNING: Computing only to max radius = %d\n",maxrad); 
	} 
	else   maxrad = 99999999; 
 
	nxy = nx * ny; 
	nyz = ny * nz; 
	nxz = nx * nz; 
	nxyz = nx * ny * nz; 
	/* FORM AND FILL TT AND SLOWNESS ARRAYS */ 
	/*
 
	if((tfint=open(timefile,O_CREAT|O_WRONLY|O_TRUNC,0664))<=1) { 
		fprintf(stderr,"cannot open %s\n",timefile); 
		exit(-1); 
	} 
	*/

	/* 
	if((vfint=open(velfile,O_RDONLY,0664))<=1) { 
		fprintf(stderr,"cannot open %s\n",velfile); 
		exit(-1); 
	} 
	if (srctype == 2)  { 
		if((wfint=open(wallfile,O_RDONLY,0664))<=1) { 
			fprintf(stderr,"cannot open %s\n",wallfile); 
			exit(-1); 
		} 
	} 
	if (srctype == 3)  { 
		if((ofint=open(oldtfile,O_RDONLY,0664))<=1) { 
			fprintf(stderr,"cannot open %s\n",oldtfile); 
			exit(-1); 
		} 
	} */ 
 
	/* ALLOCATE MAIN AND ALTERNATE GRID FOR SLOWNESSES AND TIMES */ 
	slow0 = (float *) malloc(4*nxyz); 
	time0 = (float *) malloc(4*nxyz); 
 
	/* MAKE ARRAY SORT LARGE ENOUGH FOR ANY SIDE */ 
	if(nx <= ny && nx <= nz)  { 
		sort = (struct sorted *) malloc(sizeof(struct sorted)*ny*nz); 
		nwall = nyz; 
	} 
	else if(ny <= nx && ny <= nz)  { 
		sort = (struct sorted *) malloc(sizeof(struct sorted)*nx*nz); 
		nwall = nxz; 
	} 
	else  { 
		sort = (struct sorted *) malloc(sizeof(struct sorted)*nx*ny); 
		nwall = nxy; 
	} 
	wall = (float *) malloc(4*nwall); 
	if( slow0 == NULL || time0 == NULL || sort == NULL || wall == NULL) { 
		fprintf(stderr,"cannot allocate memory\n"); 
		exit(-1); 
	} 
	/* READ IN VELOCITY FILE */ 
 
	/*read(vfint,slow0,nxyz*4);  by H. Zhang*/ 
	/* CONVERT TO SLOWNESS TIMES MESH SPACING */ 
	/* 
	for(i=0;i<nxyz;i++) slow0[i] = h/slow0[i]; */ 

	/* output velocity in punch.c--vel(5,8,i) */
	
	/*
	for(i=0;i<nx;i++){
	   printf("%f ",vfield[nx*ny*(nz-1)+nx*(ny-1)+i]);
	}
	printf("\n");
	*/
	

	for(i=0;i<nxyz;i++) slow0[i] = h/vfield[i]; 
	 
 
	/* SET TIMES TO DUMMY VALUE ***** JAH ***** BUG IN VIDALE'S CODE */ 
	for(i=0;i<nxyz;i++) time0[i] = 1.0e10; 
 
	if (srctype == 1) {			/*  VIDALE'S POINT SOURCE */ 
 
		/* FILL IN CUBE AROUND SOURCE POINT */ 
		/* VIDALE'S CONSTANT VELOCITY CUBE */ 
/*		s000 = s0(xs,ys,zs); 
		for (xx = xs-NCUBE; xx <= xs+NCUBE; xx++) { 
			if (xx < 0 || xx >= nx)	continue;  
			for (yy = ys-NCUBE; yy <= ys+NCUBE; yy++) { 
				if (yy < 0 || yy >= ny)	continue;  
				for (zz = zs-NCUBE; zz <= zs+NCUBE; zz++) { 
					if (zz < 0 || zz >= nz)	continue;  
					t0(xx, yy, zz) = s000 * DIST(fxs, fys, fzs, xx, yy, zz); 
				} 
			} 
		} 
*/ 
 
		/* HOLE'S OLD LINEAR VELOCITY GRADIENT CUBE (SEPT 1990)*/ 
/*		xx = xs; 
		while (xx>0 && xx>xs-NCUBE)	xx--; 
		xxx = xs; 
		while (xxx<nx-1 && xxx<xs+NCUBE)	xxx++; 
		dvx = (1/s0(xxx,ys,zs) - 1/s0(xx,ys,zs))/(xxx-xx); 
		fprintf(stderr, "source cube dvx = %f\n",dvx); 
		yy = ys; 
		while (yy>0 && yy>ys-NCUBE)	yy--; 
		yyy = ys; 
		while (yyy<ny-1 && yyy<ys+NCUBE)	yyy++; 
		dvy = (1/s0(xs,yyy,zs) - 1/s0(xs,yy,zs))/(yyy-yy); 
		fprintf(stderr, "source cube dvy = %f\n",dvy); 
		zz = zs; 
		while (zz>0 && zz>zs-NCUBE)	zz--; 
		zzz = zs; 
		while (zzz<nz-1 && zzz<zs+NCUBE)	zzz++; 
		dvz = (1/s0(xs,ys,zzz) - 1/s0(xs,ys,zz))/(zzz-zz); 
		fprintf(stderr, "source cube dvz = %f\n",dvz); 
		dv = sqrt(dvx*dvx + dvy*dvy + dvz*dvz); 
		v0 = h/s0(xs,ys,zs); 
		if (dv != 0.) rzc = -v0/dv; 
		fprintf(stderr, "source cube rzc = %f\n",rzc); 
		for (xx = xs-NCUBE; xx <= xs+NCUBE; xx++) { 
			if (xx < 0 || xx >= nx)	continue;  
			for (yy = ys-NCUBE; yy <= ys+NCUBE; yy++) { 
				if (yy < 0 || yy >= ny)	continue;  
				for (zz = zs-NCUBE; zz <= zs+NCUBE; zz++) { 
					if (zz < 0 || zz >= nz)	continue;  
					if (dv == 0.)  { 
					  t0(xx,yy,zz) = s0(xs,ys,zs)*DIST(fxs,fys,fzs,xx,yy,zz); 
					  continue; 
					} 
					rx = h*(xx - fxs); 
					ry = h*(yy - fys); 
					rz = h*(zz - fzs); 
					rz1 = (rx*dvx+ry*dvy+rz*dvz)/dv; 
					rxy1 = sqrt(rx*rx+ry*ry+rz*rz-rz1*rz1); 
					if (rxy1<=h/1.e6) 
					  t0(xx,yy,zz) = fabs(log((v0+dv*rz1)/v0)/dv); 
					else { 
					  rxyc = (rz1*rz1+rxy1*rxy1-2*rz1*rzc)/(2*rxy1); 
					  rho = sqrt(rzc*rzc+rxyc*rxyc); 
					  theta1 = asin(-rzc/rho); 
					  theta2 = asin((rz1-rzc)/rho); 
					  if (rxyc<0) theta1=PI-theta1; 
					  if (rxyc<rxy1) theta2=PI-theta2; 
					  t0(xx,yy,zz) = log(tan(theta2/2)/tan(theta1/2)) / dv; 
				        } 
				} 
			} 
		} 
*/ 
 
		/* HOLE'S NEW LINEAR VELOCITY GRADIENT CUBE (APRIL 1991)*/ 
		v0 = h/s0(xs,ys,zs); 
		for (xx = xs-NCUBE; xx <= xs+NCUBE; xx++) { 
			if (xx < 0 || xx >= nx)	continue;  
			for (yy = ys-NCUBE; yy <= ys+NCUBE; yy++) { 
				if (yy < 0 || yy >= ny)	continue;  
				for (zz = zs-NCUBE; zz <= zs+NCUBE; zz++) { 
					if (zz < 0 || zz >= nz)	continue;  
					if (zz == zs) 
					  dvz = 1/s0(xx,yy,zz+1)-1/s0(xs,ys,zs); 
					else 
					  dvz = (1/s0(xx,yy,zz)-1/s0(xs,ys,zs))/(zz-zs); 
					dv = fabs(dvz); 
					if (dv == 0.)  { 
					  t0(xx,yy,zz) = s0(xs,ys,zs)*DIST(fxs,fys,fzs,xx,yy,zz); 
					  continue; 
					} 
					rzc = -v0/dv; 
					rx = h*(xx - fxs); 
					ry = h*(yy - fys); 
					rz = h*(zz - fzs); 
					rz1 = rz*dvz/dv; 
					rxy1 = sqrt(rx*rx+ry*ry+rz*rz-rz1*rz1); 
					if (rxy1<=h/1.e6) 
					  t0(xx,yy,zz) = fabs(log((v0+dv*rz1)/v0)/dv); 
					else { 
					  rxyc = (rz1*rz1+rxy1*rxy1-2*rz1*rzc)/(2*rxy1); 
					  rho = sqrt(rzc*rzc+rxyc*rxyc); 
					  theta1 = asin(-rzc/rho); 
					  /* can't handle asin(1.) ! */ 
					  if (fabs(rz1-rzc)>=rho)  rho=1.0000001*fabs(rz1-rzc); 
					  theta2 = asin((rz1-rzc)/rho); 
					  if (rxyc<0) theta1=PI-theta1; 
					  if (rxyc<rxy1) theta2=PI-theta2; 
					  t0(xx,yy,zz) = log(tan(theta2/2)/tan(theta1/2)) / dv; 
				        } 
				} 
			} 
		} 
 
		/* SETS LOCATION OF THE SIDES OF THE CUBE	*/ 
		radius = NCUBE; 
		if(xs > NCUBE) x0 = xs - (NCUBE + 1); 
		else{ x0 = -1; dx0 = 0;} 
		if(xs < nx-(NCUBE + 1)) x1 = xs + (NCUBE + 1); 
		else{ x1 = nx; dx1 = 0;} 
		if(ys > NCUBE) y0 = ys - (NCUBE + 1); 
		else{ y0 = -1; dy0 = 0;} 
		if(ys < ny-(NCUBE + 1)) y1 = ys + (NCUBE + 1); 
		else{ y1 = ny; dy1 = 0;} 
		if(zs > NCUBE) z0 = zs - (NCUBE + 1); 
		else{ z0 = -1; dz0 = 0;} 
		if(zs < nz-(NCUBE + 1)) z1 = zs + (NCUBE + 1); 
		else{ z1 = nz; dz1 = 0;} 
	} 
	else if (srctype == 2) {		/*  HOLE'S EXTERNAL SOURCE */ 
 
		/* FILL IN WALLS' TIMES FROM EXTERNAL DATAFILE */ 
		read (wfint,wall,4*nwall);	/* READ X=0 WALL */ 
		if (wall[0]>-1.e-20) { 
			ii = 0; 
			for (k=0; k<nz; k++) { 
				for (j=0; j<ny; j++) { 
					t0(0,j,k) = wall[ii]; 
					ii++; 
				} 
			} 
		} 
		read (wfint,wall,4*nwall);	/* READ X=NX-1 WALL */ 
		if (wall[0]>-1.e-20) { 
			ii = 0; 
			for (k=0; k<nz; k++) { 
				for (j=0; j<ny; j++) { 
					t0(nx-1,j,k) = wall[ii]; 
					ii++; 
				} 
			} 
		} 
		read (wfint,wall,4*nwall);	/* READ Y=0 WALL */ 
		if (wall[0]>-1.e-20) { 
			ii = 0; 
			for (k=0; k<nz; k++) { 
				for (i=0; i<nx; i++) { 
					t0(i,0,k) = wall[ii]; 
					ii++; 
				} 
			} 
		} 
		read (wfint,wall,4*nwall);	/* READ Y=NY-1 WALL */ 
		if (wall[0]>-1.e-20) { 
			ii = 0; 
			for (k=0; k<nz; k++) { 
				for (i=0; i<nx; i++) { 
					t0(i,ny-1,k) = wall[ii]; 
					ii++; 
				} 
			} 
		} 
		read (wfint,wall,4*nwall);	/* READ Z=0 WALL */ 
		if (wall[0]>-1.e-20) { 
			ii = 0; 
			for (j=0; j<ny; j++) { 
				for (i=0; i<nx; i++) { 
					t0(i,j,0) = wall[ii]; 
					ii++; 
				} 
			} 
		} 
		read (wfint,wall,4*nwall);	/* READ Z=NZ-1 WALL */ 
		if (wall[0]>-1.e-20) { 
			ii = 0; 
			for (j=0; j<ny; j++) { 
				for (i=0; i<nx; i++) { 
					t0(i,j,nz-1) = wall[ii]; 
					ii++; 
				} 
			} 
		} 
 
		/* SET LOCATIONS OF SIDES OF THE CUBE SO THAT CUBE IS A FACE  */ 
		radius = 1; 
		if (srcwall == 1)	x1=1; 
		else	{  x1=nx;	dx1=0;  } 
		if (srcwall == 2)	x0=nx-2; 
		else	{  x0= -1;	dx0=0;  } 
		if (srcwall == 3)	y1=1; 
		else	{  y1=ny;	dy1=0;  } 
		if (srcwall == 4)	y0=ny-2; 
		else	{  y0= -1;	dy0=0;  } 
		if (srcwall == 5)	z1=1; 
		else	{  z1=nz;	dz1=0;  } 
		if (srcwall == 6)	z0=nz-2; 
		else	{  z0= -1;	dz0=0;  } 
	} 
	else if (srctype == 3) {                /*  HOLE'S REDO OLD TIMES */ 
	        /* READ IN OLD TIME FILE */ 
	        if (srctype == 3)  read(ofint,time0,nxyz*4); 
		/* SET LOCATIONS OF SIDES OF THE CUBE SO THAT CUBE IS A FACE */ 
		radius = 1; 
		if (srcwall == 1)	x1=1; 
		else	{  x1=nx;	dx1=0;  } 
		if (srcwall == 2)	x0=nx-2; 
		else	{  x0= -1;	dx0=0;  } 
		if (srcwall == 3)	y1=1; 
		else	{  y1=ny;	dy1=0;  } 
		if (srcwall == 4)	y0=ny-2; 
		else	{  y0= -1;	dy0=0;  } 
		if (srcwall == 5)	z1=1; 
		else	{  z1=nz;	dz1=0;  } 
		if (srcwall == 6)	z0=nz-2; 
		else	{  z0= -1;	dz0=0;  } 
	} 
	else  { 
		fprintf(stderr,"incorrect value of srctype = %d\n",srctype); 
		exit(-1); 
	} 
 
	if (headpref>0) {	/* HOLE - PREFERRED REVERSE DIRECTION */ 
		head = nx*ny*nz; 
		if (nx>5 && x1<=head)   {headpref=2;  head=x1;} 
		if (nx>5 && (nx-1-x0)<=head)   {headpref=1;  head=nx-1-x0;} 
		if (ny>5 && y1<=head)   {headpref=4;  head=y1;} 
		if (ny>5 && (ny-1-y0)<=head)   {headpref=3;  head=ny-1-y0;} 
		if (nz>5 && z1<=head)   {headpref=6;  head=z1;} 
		if (nz>5 && (nz-1-z0)<=head)   {headpref=5;  head=nz-1-z0;} 
	} 
 
	/* BIGGER LOOP - HOLE - ALLOWS AUTOMATIC REVERSE PROPAGATION IF  
		HEAD WAVES ARE ENCOUNTERED ON FACES OF EXPANDING CUBE,  
		ALLOWING WAVES TO TRAVEL BACK INTO THE CUBE */ 
	while ( reverse > -1 )  { 
 
		headw[1]=0; headw[2]=0; headw[3]=0; headw[4]=0; 
		headw[5]=0; headw[6]=0; 
 
	/* BIG LOOP */ 
	while(rad0 && (dx0 || dx1 || dy0 || dy1 || dz0 || dz1))  { 
 
		/* CALCULATE ON PRIMARY (time0) GRID */ 
 
		/* TOP SIDE */ 
      for (igrow=1;igrow<=kminus;igrow++) {   
	if(dz0){ 
		ii = 0; 
		for(j=y0+1; j<=y1-1; j++){ 
			for(i=x0+1; i<=x1-1; i++){ 
				sort[ii].time = t0(i,j,z0+1); 
				sort[ii].x1 = i; 
				sort[ii].x2 = j; 
				ii++; 
			} 
		} 
		qsort((char *)sort,ii,sizeof(struct sorted),compar); 
		for(i=0;i<ii;i++){ 
			X1 = sort[i].x1; 
			X2 = sort[i].x2; 
			index = z0*nxy + X2*nx + X1; 
			lasti = (z0+1)*nxy + X2*nx + X1; 
			fhead = 0.; 
			guess = time0[index]; 
                        if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9 
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) { 
			  try = fdh3d(              t0(X1,X2,z0+1), 
				      t0(X1+1,X2,z0+1),t0(X1+1,X2+1,z0+1),t0(X1,X2+1,z0+1), 
				      t0(X1+1,X2,z0  ),t0(X1+1,X2+1,z0  ),t0(X1,X2+1,z0  ), 
				      s0(X1,X2,z0), s0(X1,X2,z0+1), 
				      s0(X1+1,X2,z0+1),s0(X1+1,X2+1,z0+1),s0(X1,X2+1,z0+1), 
				      s0(X1+1,X2,z0  ),s0(X1+1,X2+1,z0  ),s0(X1,X2+1,z0  )); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9 
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) { 
			  try = fdh3d(              t0(X1,X2,z0+1), 
				      t0(X1-1,X2,z0+1),t0(X1-1,X2+1,z0+1),t0(X1,X2+1,z0+1), 
				      t0(X1-1,X2,z0  ),t0(X1-1,X2+1,z0  ),t0(X1,X2+1,z0  ), 
				      s0(X1,X2,z0), s0(X1,X2,z0+1), 
				      s0(X1-1,X2,z0+1),s0(X1-1,X2+1,z0+1),s0(X1,X2+1,z0+1), 
				      s0(X1-1,X2,z0  ),s0(X1-1,X2+1,z0  ),s0(X1,X2+1,z0  )); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9 
			   && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) { 
			  try = fdh3d(              t0(X1,X2,z0+1), 
				      t0(X1+1,X2,z0+1),t0(X1+1,X2-1,z0+1),t0(X1,X2-1,z0+1), 
				      t0(X1+1,X2,z0  ),t0(X1+1,X2-1,z0  ),t0(X1,X2-1,z0  ), 
				      s0(X1,X2,z0), s0(X1,X2,z0+1), 
				      s0(X1+1,X2,z0+1),s0(X1+1,X2-1,z0+1),s0(X1,X2-1,z0+1), 
				      s0(X1+1,X2,z0  ),s0(X1+1,X2-1,z0  ),s0(X1,X2-1,z0  )); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9 
			   && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) { 
			  try = fdh3d(              t0(X1,X2,z0+1), 
				      t0(X1-1,X2,z0+1),t0(X1-1,X2-1,z0+1),t0(X1,X2-1,z0+1), 
				      t0(X1-1,X2,z0  ),t0(X1-1,X2-1,z0  ),t0(X1,X2-1,z0  ), 
				      s0(X1,X2,z0), s0(X1,X2,z0+1), 
				      s0(X1-1,X2,z0+1),s0(X1-1,X2-1,z0+1),s0(X1,X2-1,z0+1), 
				      s0(X1-1,X2,z0  ),s0(X1-1,X2-1,z0  ),s0(X1,X2-1,z0  )); 
			  if (try<guess) guess = try; 
			} 
			if(guess > 1.0e9){  
			  if(time0[index+1] < 1.e9 && X1<nx-1 && X2>y0+1 && X2<y1-1 )  { 
			      try = fdhne(t0(X1,X2,z0+1),t0(X1+1,X2,z0+1),t0(X1+1,X2,z0), 
					  t0(X1+1,X2-1,z0+1),t0(X1+1,X2+1,z0+1), 
					  s0(X1,X2,z0), 
					  s0(X1,X2,z0+1),s0(X1+1,X2,z0+1),s0(X1+1,X2,z0) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-1] < 1.e9 && X1>0 && X2>y0+1 && X2<y1-1 )  { 
			      try = fdhne(t0(X1,X2,z0+1),t0(X1-1,X2,z0+1),t0(X1-1,X2,z0), 
					  t0(X1-1,X2-1,z0+1),t0(X1-1,X2+1,z0+1), 
					  s0(X1,X2,z0), 
					  s0(X1,X2,z0+1),s0(X1-1,X2,z0+1),s0(X1-1,X2,z0) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nx] < 1.e9 && X2<ny-1 && X1>x0+1 && X1<x1-1 )  { 
			      try = fdhne(t0(X1,X2,z0+1),t0(X1,X2+1,z0+1),t0(X1,X2+1,z0), 
					  t0(X1-1,X2+1,z0+1),t0(X1+1,X2+1,z0+1), 
					  s0(X1,X2,z0), 
					  s0(X1,X2,z0+1),s0(X1,X2+1,z0+1),s0(X1,X2+1,z0) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nx] < 1.e9 && X2>0 && X1>x0+1 && X1<x1-1 )  { 
			      try = fdhne(t0(X1,X2,z0+1),t0(X1,X2-1,z0+1),t0(X1,X2-1,z0), 
					  t0(X1-1,X2-1,z0+1),t0(X1+1,X2-1,z0+1), 
					  s0(X1,X2,z0), 
					  s0(X1,X2,z0+1),s0(X1,X2-1,z0+1),s0(X1,X2-1,z0) ); 
			    if (try<guess)  guess = try; 
			  } 
		        }  
			  if(time0[index+1] < 1.e9 && X1<nx-1 )  { 
			    try = fdh2d(t0(X1,X2,z0+1),t0(X1+1,X2,z0+1),t0(X1+1,X2,z0), 
					  s0(X1,X2,z0), 
					  s0(X1,X2,z0+1),s0(X1+1,X2,z0+1),s0(X1+1,X2,z0) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-1] < 1.e9 && X1>0 )  { 
			    try = fdh2d(t0(X1,X2,z0+1),t0(X1-1,X2,z0+1),t0(X1-1,X2,z0), 
					  s0(X1,X2,z0), 
					  s0(X1,X2,z0+1),s0(X1-1,X2,z0+1),s0(X1-1,X2,z0) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nx] < 1.e9 && X2<ny-1 )  { 
			    try = fdh2d(t0(X1,X2,z0+1),t0(X1,X2+1,z0+1),t0(X1,X2+1,z0), 
					  s0(X1,X2,z0), 
					  s0(X1,X2,z0+1),s0(X1,X2+1,z0+1),s0(X1,X2+1,z0) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nx] < 1.e9 && X2>0 )  { 
			    try = fdh2d(t0(X1,X2,z0+1),t0(X1,X2-1,z0+1),t0(X1,X2-1,z0), 
					  s0(X1,X2,z0), 
					  s0(X1,X2,z0+1),s0(X1,X2-1,z0+1),s0(X1,X2-1,z0) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9 
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) { 
			    try = fdh2d(t0(X1+1,X2,z0),t0(X1+1,X2+1,z0),t0(X1,X2+1,z0), 
					s0(X1,X2,z0), 
					s0(X1+1,X2,z0),s0(X1+1,X2+1,z0),s0(X1,X2+1,z0) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9 
			     && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) { 
			    try = fdh2d(t0(X1+1,X2,z0),t0(X1+1,X2-1,z0),t0(X1,X2-1,z0), 
					s0(X1,X2,z0), 
					s0(X1+1,X2,z0),s0(X1+1,X2-1,z0),s0(X1,X2-1,z0) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9 
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) { 
			    try = fdh2d(t0(X1-1,X2,z0),t0(X1-1,X2+1,z0),t0(X1,X2+1,z0), 
					s0(X1,X2,z0), 
					s0(X1-1,X2,z0),s0(X1-1,X2+1,z0),s0(X1,X2+1,z0) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9 
			     && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) { 
			    try = fdh2d(t0(X1-1,X2,z0),t0(X1-1,X2-1,z0),t0(X1,X2-1,z0), 
					s0(X1,X2,z0), 
					s0(X1-1,X2,z0),s0(X1-1,X2-1,z0),s0(X1,X2-1,z0) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if(guess > 1.0e9){  
			  if ( X1>x0+1 && X1<x1-1 && X2>y0+1 && X2<y1-1 ) { 
			    try = fdhnf(t0(X1,X2,z0+1), 
					  t0(X1+1,X2,z0+1),t0(X1,X2+1,z0+1), 
					  t0(X1-1,X2,z0+1),t0(X1,X2-1,z0+1), 
					  s0(X1,X2,z0), 
					  s0(X1,X2,z0+1) ); 
			    if (try<guess)  guess = try; 
			  } 
			}  
                          try = t0(X1,X2,z0+1) + .5*(s0(X1,X2,z0)+s0(X1,X2,z0+1)); 
			  if (try<guess)  guess = try; 
                          if ( time0[index+1]<1.e9 && X1<nx-1 )  { 
			    try = t0(X1+1,X2,z0) + .5*(s0(X1,X2,z0)+s0(X1+1,X2,z0)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-1]<1.e9 && X1>0 )  { 
			    try = t0(X1-1,X2,z0) + .5*(s0(X1,X2,z0)+s0(X1-1,X2,z0)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index+nx]<1.e9 && X2<ny-1 )  { 
			    try = t0(X1,X2+1,z0) + .5*(s0(X1,X2,z0)+s0(X1,X2+1,z0)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-nx]<1.e9 && X2>0 )  { 
			    try = t0(X1,X2-1,z0) + .5*(s0(X1,X2,z0)+s0(X1,X2-1,z0)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if (guess<time0[index])  { 
				time0[index] = guess; 
				if (fhead>headtest)  headw[5]++; 
			} 
		} 
		if(z0 == 0) dz0 = 0; 
		z0--; 
	} 
      } 
		/* BOTTOM SIDE */ 
      for (igrow=1;igrow<=kplus;igrow++) {   
	if(dz1){ 
		ii = 0; 
		for(j=y0+1; j<=y1-1; j++){ 
			for(i=x0+1; i<=x1-1; i++){ 
				sort[ii].time = t0(i,j,z1-1); 
				sort[ii].x1 = i; 
				sort[ii].x2 = j; 
				ii++; 
			} 
		} 
		qsort((char *)sort,ii,sizeof(struct sorted),compar); 
		for(i=0;i<ii;i++){ 
			X1 = sort[i].x1; 
			X2 = sort[i].x2; 
			index = z1*nxy + X2*nx + X1; 
			lasti = (z1-1)*nxy + X2*nx + X1; 
			fhead = 0.; 
			guess = time0[index]; 
                        if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9 
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) { 
			  try = fdh3d(              t0(X1,X2,z1-1), 
				      t0(X1+1,X2,z1-1),t0(X1+1,X2+1,z1-1),t0(X1,X2+1,z1-1), 
				      t0(X1+1,X2,z1  ),t0(X1+1,X2+1,z1  ),t0(X1,X2+1,z1  ), 
				      s0(X1,X2,z1), s0(X1,X2,z1-1), 
				      s0(X1+1,X2,z1-1),s0(X1+1,X2+1,z1-1),s0(X1,X2+1,z1-1), 
				      s0(X1+1,X2,z1  ),s0(X1+1,X2+1,z1  ),s0(X1,X2+1,z1  )); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9 
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) { 
			  try = fdh3d(              t0(X1,X2,z1-1), 
				      t0(X1-1,X2,z1-1),t0(X1-1,X2+1,z1-1),t0(X1,X2+1,z1-1), 
				      t0(X1-1,X2,z1  ),t0(X1-1,X2+1,z1  ),t0(X1,X2+1,z1  ), 
				      s0(X1,X2,z1), s0(X1,X2,z1-1), 
				      s0(X1-1,X2,z1-1),s0(X1-1,X2+1,z1-1),s0(X1,X2+1,z1-1), 
				      s0(X1-1,X2,z1  ),s0(X1-1,X2+1,z1  ),s0(X1,X2+1,z1  )); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9 
			   && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) { 
			  try = fdh3d(              t0(X1,X2,z1-1), 
				      t0(X1+1,X2,z1-1),t0(X1+1,X2-1,z1-1),t0(X1,X2-1,z1-1), 
				      t0(X1+1,X2,z1  ),t0(X1+1,X2-1,z1  ),t0(X1,X2-1,z1  ), 
				      s0(X1,X2,z1), s0(X1,X2,z1-1), 
				      s0(X1+1,X2,z1-1),s0(X1+1,X2-1,z1-1),s0(X1,X2-1,z1-1), 
				      s0(X1+1,X2,z1  ),s0(X1+1,X2-1,z1  ),s0(X1,X2-1,z1  )); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9 
			   && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) { 
			  try = fdh3d(              t0(X1,X2,z1-1), 
				      t0(X1-1,X2,z1-1),t0(X1-1,X2-1,z1-1),t0(X1,X2-1,z1-1), 
				      t0(X1-1,X2,z1  ),t0(X1-1,X2-1,z1  ),t0(X1,X2-1,z1  ), 
				      s0(X1,X2,z1), s0(X1,X2,z1-1), 
				      s0(X1-1,X2,z1-1),s0(X1-1,X2-1,z1-1),s0(X1,X2-1,z1-1), 
				      s0(X1-1,X2,z1  ),s0(X1-1,X2-1,z1  ),s0(X1,X2-1,z1  )); 
			  if (try<guess) guess = try; 
			} 
                        if(guess > 1.0e9){  
			  if(time0[index+1] < 1.e9 && X1<nx-1 && X2>y0+1 && X2<y1-1 )  { 
			      try = fdhne(t0(X1,X2,z1-1),t0(X1+1,X2,z1-1),t0(X1+1,X2,z1), 
					  t0(X1+1,X2-1,z1-1),t0(X1+1,X2+1,z1-1), 
					  s0(X1,X2,z1), 
					  s0(X1,X2,z1-1),s0(X1+1,X2,z1-1),s0(X1+1,X2,z1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-1] < 1.e9 && X1>0 && X2>y0+1 && X2<y1-1 )  { 
			      try = fdhne(t0(X1,X2,z1-1),t0(X1-1,X2,z1-1),t0(X1-1,X2,z1), 
					  t0(X1-1,X2-1,z1-1),t0(X1-1,X2+1,z1-1), 
					  s0(X1,X2,z1), 
					  s0(X1,X2,z1-1),s0(X1-1,X2,z1-1),s0(X1-1,X2,z1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nx] < 1.e9 && X2<ny-1 && X1>x0+1 && X1<x1-1 )  { 
			      try = fdhne(t0(X1,X2,z1-1),t0(X1,X2+1,z1-1),t0(X1,X2+1,z1), 
					  t0(X1-1,X2+1,z1-1),t0(X1+1,X2+1,z1-1), 
					  s0(X1,X2,z1), 
					  s0(X1,X2,z1-1),s0(X1,X2+1,z1-1),s0(X1,X2+1,z1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nx] < 1.e9 && X2>0 && X1>x0+1 && X1<x1-1 )  { 
			      try = fdhne(t0(X1,X2,z1-1),t0(X1,X2-1,z1-1),t0(X1,X2-1,z1), 
					  t0(X1-1,X2-1,z1-1),t0(X1+1,X2-1,z1-1), 
					  s0(X1,X2,z1), 
					  s0(X1,X2,z1-1),s0(X1,X2-1,z1-1),s0(X1,X2-1,z1) ); 
			    if (try<guess)  guess = try; 
			  } 
		        } 
			  if(time0[index+1] < 1.e9 && X1<nx-1 )  { 
			    try = fdh2d(t0(X1,X2,z1-1),t0(X1+1,X2,z1-1),t0(X1+1,X2,z1), 
					  s0(X1,X2,z1), 
					  s0(X1,X2,z1-1),s0(X1+1,X2,z1-1),s0(X1+1,X2,z1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-1] < 1.e9 && X1>0 )  { 
			    try = fdh2d(t0(X1,X2,z1-1),t0(X1-1,X2,z1-1),t0(X1-1,X2,z1), 
					  s0(X1,X2,z1), 
					  s0(X1,X2,z1-1),s0(X1-1,X2,z1-1),s0(X1-1,X2,z1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nx] < 1.e9 && X2<ny-1 )  { 
			    try = fdh2d(t0(X1,X2,z1-1),t0(X1,X2+1,z1-1),t0(X1,X2+1,z1), 
					  s0(X1,X2,z1), 
					  s0(X1,X2,z1-1),s0(X1,X2+1,z1-1),s0(X1,X2+1,z1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nx] < 1.e9 && X2>0 )  { 
			    try = fdh2d(t0(X1,X2,z1-1),t0(X1,X2-1,z1-1),t0(X1,X2-1,z1), 
					  s0(X1,X2,z1), 
					  s0(X1,X2,z1-1),s0(X1,X2-1,z1-1),s0(X1,X2-1,z1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9 
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) { 
			    try = fdh2d(t0(X1+1,X2,z1),t0(X1+1,X2+1,z1),t0(X1,X2+1,z1), 
					s0(X1,X2,z1), 
					s0(X1+1,X2,z1),s0(X1+1,X2+1,z1),s0(X1,X2+1,z1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9 
			     && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) { 
			    try = fdh2d(t0(X1+1,X2,z1),t0(X1+1,X2-1,z1),t0(X1,X2-1,z1), 
					s0(X1,X2,z1), 
					s0(X1+1,X2,z1),s0(X1+1,X2-1,z1),s0(X1,X2-1,z1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9 
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) { 
			    try = fdh2d(t0(X1-1,X2,z1),t0(X1-1,X2+1,z1),t0(X1,X2+1,z1), 
					s0(X1,X2,z1), 
					s0(X1-1,X2,z1),s0(X1-1,X2+1,z1),s0(X1,X2+1,z1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9 
			     && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) { 
			    try = fdh2d(t0(X1-1,X2,z1),t0(X1-1,X2-1,z1),t0(X1,X2-1,z1), 
					s0(X1,X2,z1), 
					s0(X1-1,X2,z1),s0(X1-1,X2-1,z1),s0(X1,X2-1,z1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if(guess > 1.0e9){  
			  if ( X1>x0+1 && X1<x1-1 && X2>y0+1 && X2<y1-1 ) { 
			    try = fdhnf(t0(X1,X2,z1-1), 
					  t0(X1+1,X2,z1-1),t0(X1,X2+1,z1-1), 
					  t0(X1-1,X2,z1-1),t0(X1,X2-1,z1-1), 
					  s0(X1,X2,z1), 
					  s0(X1,X2,z1-1) ); 
			    if (try<guess)  guess = try; 
			  } 
			}  
			  try = t0(X1,X2,z1-1) + .5*(s0(X1,X2,z1)+s0(X1,X2,z1-1)); 
			  if (try<guess)  guess = try; 
                          if ( time0[index+1]<1.e9 && X1<nx-1 )  { 
			    try = t0(X1+1,X2,z1) + .5*(s0(X1,X2,z1)+s0(X1+1,X2,z1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-1]<1.e9 && X1>0 )  { 
			    try = t0(X1-1,X2,z1) + .5*(s0(X1,X2,z1)+s0(X1-1,X2,z1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index+nx]<1.e9 && X2<ny-1 )  { 
			    try = t0(X1,X2+1,z1) + .5*(s0(X1,X2,z1)+s0(X1,X2+1,z1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-nx]<1.e9 && X2>0 )  { 
			    try = t0(X1,X2-1,z1) + .5*(s0(X1,X2,z1)+s0(X1,X2-1,z1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if (guess<time0[index]) { 
				time0[index] = guess; 
				if (fhead>headtest)  headw[6]++; 
			} 
		} 
		if(z1 == nz-1) dz1 = 0; 
		z1++; 
	} 
      } 
		/* FRONT SIDE */ 
      for (igrow=1;igrow<=jminus;igrow++) {   
	if(dy0){ 
		ii = 0; 
		for(k=z0+1; k<=z1-1; k++){ 
			for(i=x0+1; i<=x1-1; i++){ 
				sort[ii].time = t0(i,y0+1,k); 
				sort[ii].x1 = i; 
				sort[ii].x2 = k; 
				ii++; 
			} 
		} 
		qsort((char *)sort,ii,sizeof(struct sorted),compar); 
		for(i=0;i<ii;i++){ 
			X1 = sort[i].x1; 
			X2 = sort[i].x2; 
			index = X2*nxy + y0*nx + X1; 
			lasti = X2*nxy + (y0+1)*nx + X1; 
			fhead = 0.; 
			guess = time0[index]; 
			if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9 
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) { 
			  try = fdh3d(              t0(X1,y0+1,X2), 
				      t0(X1+1,y0+1,X2),t0(X1+1,y0+1,X2+1),t0(X1,y0+1,X2+1), 
				      t0(X1+1,y0  ,X2),t0(X1+1,y0  ,X2+1),t0(X1,y0  ,X2+1), 
				      s0(X1,y0,X2), s0(X1,y0+1,X2), 
				      s0(X1+1,y0+1,X2),s0(X1+1,y0+1,X2+1),s0(X1,y0+1,X2+1), 
				      s0(X1+1,y0  ,X2),s0(X1+1,y0  ,X2+1),s0(X1,y0  ,X2+1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9 
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { 
			  try = fdh3d(              t0(X1,y0+1,X2), 
				      t0(X1-1,y0+1,X2),t0(X1-1,y0+1,X2+1),t0(X1,y0+1,X2+1), 
				      t0(X1-1,y0  ,X2),t0(X1-1,y0  ,X2+1),t0(X1,y0  ,X2+1), 
				      s0(X1,y0,X2), s0(X1,y0+1,X2), 
				      s0(X1-1,y0+1,X2),s0(X1-1,y0+1,X2+1),s0(X1,y0+1,X2+1), 
				      s0(X1-1,y0  ,X2),s0(X1-1,y0  ,X2+1),s0(X1,y0  ,X2+1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9 
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) { 
			  try = fdh3d(              t0(X1,y0+1,X2), 
				      t0(X1+1,y0+1,X2),t0(X1+1,y0+1,X2-1),t0(X1,y0+1,X2-1), 
				      t0(X1+1,y0  ,X2),t0(X1+1,y0  ,X2-1),t0(X1,y0  ,X2-1), 
				      s0(X1,y0,X2), s0(X1,y0+1,X2), 
				      s0(X1+1,y0+1,X2),s0(X1+1,y0+1,X2-1),s0(X1,y0+1,X2-1), 
				      s0(X1+1,y0  ,X2),s0(X1+1,y0  ,X2-1),s0(X1,y0  ,X2-1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9 
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) { 
			  try = fdh3d(              t0(X1,y0+1,X2), 
				      t0(X1-1,y0+1,X2),t0(X1-1,y0+1,X2-1),t0(X1,y0+1,X2-1), 
				      t0(X1-1,y0  ,X2),t0(X1-1,y0  ,X2-1),t0(X1,y0  ,X2-1), 
				      s0(X1,y0,X2), s0(X1,y0+1,X2), 
				      s0(X1-1,y0+1,X2),s0(X1-1,y0+1,X2-1),s0(X1,y0+1,X2-1), 
				      s0(X1-1,y0  ,X2),s0(X1-1,y0  ,X2-1),s0(X1,y0  ,X2-1)); 
			  if (try<guess) guess = try; 
			} 
			if(guess > 1.0e9){  
			  if(time0[index+1] < 1.e9 && X1<nx-1 && X2>z0+1 && X2<z1-1 )  { 
			      try = fdhne(t0(X1,y0+1,X2),t0(X1+1,y0+1,X2),t0(X1+1,y0,X2), 
					  t0(X1+1,y0+1,X2-1),t0(X1+1,y0+1,X2+1), 
					  s0(X1,y0,X2), 
					  s0(X1,y0+1,X2),s0(X1+1,y0+1,X2),s0(X1+1,y0,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-1] < 1.e9 && X1>0 && X2>z0+1 && X2<z1-1 )  { 
			      try = fdhne(t0(X1,y0+1,X2),t0(X1-1,y0+1,X2),t0(X1-1,y0,X2), 
					  t0(X1-1,y0+1,X2-1),t0(X1-1,y0+1,X2+1), 
					  s0(X1,y0,X2), 
					  s0(X1,y0+1,X2),s0(X1-1,y0+1,X2),s0(X1-1,y0,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>x0+1 && X1<x1-1 )  { 
			      try = fdhne(t0(X1,y0+1,X2),t0(X1,y0+1,X2+1),t0(X1,y0,X2+1), 
					  t0(X1-1,y0+1,X2+1),t0(X1+1,y0+1,X2+1), 
					  s0(X1,y0,X2), 
					  s0(X1,y0+1,X2),s0(X1,y0+1,X2+1),s0(X1,y0,X2+1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nxy] < 1.e9 && X2>0 && X1>x0+1 && X1<x1-1 )  { 
			      try = fdhne(t0(X1,y0+1,X2),t0(X1,y0+1,X2-1),t0(X1,y0,X2-1), 
					  t0(X1-1,y0+1,X2-1),t0(X1+1,y0+1,X2-1), 
					  s0(X1,y0,X2), 
					  s0(X1,y0+1,X2),s0(X1,y0+1,X2-1),s0(X1,y0,X2-1) ); 
			    if (try<guess)  guess = try; 
			  } 
		        }  
			  if(time0[index+1] < 1.e9 && X1<nx-1 )  { 
			    try = fdh2d(t0(X1,y0+1,X2),t0(X1+1,y0+1,X2),t0(X1+1,y0,X2), 
					  s0(X1,y0,X2), 
					  s0(X1,y0+1,X2),s0(X1+1,y0+1,X2),s0(X1+1,y0,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-1] < 1.e9 && X1>0 )  { 
			    try = fdh2d(t0(X1,y0+1,X2),t0(X1-1,y0+1,X2),t0(X1-1,y0,X2), 
					  s0(X1,y0,X2), 
					  s0(X1,y0+1,X2),s0(X1-1,y0+1,X2),s0(X1-1,y0,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 )  { 
			    try = fdh2d(t0(X1,y0+1,X2),t0(X1,y0+1,X2+1),t0(X1,y0,X2+1), 
					  s0(X1,y0,X2), 
					  s0(X1,y0+1,X2),s0(X1,y0+1,X2+1),s0(X1,y0,X2+1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nxy] < 1.e9 && X2>0 )  { 
			    try = fdh2d(t0(X1,y0+1,X2),t0(X1,y0+1,X2-1),t0(X1,y0,X2-1), 
					  s0(X1,y0,X2), 
					  s0(X1,y0+1,X2),s0(X1,y0+1,X2-1),s0(X1,y0,X2-1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9 
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) { 
			    try = fdh2d(t0(X1+1,y0,X2),t0(X1+1,y0,X2+1),t0(X1,y0,X2+1), 
					s0(X1,y0,X2), 
					s0(X1+1,y0,X2),s0(X1+1,y0,X2+1),s0(X1,y0,X2+1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9 
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) { 
			    try = fdh2d(t0(X1+1,y0,X2),t0(X1+1,y0,X2-1),t0(X1,y0,X2-1), 
					s0(X1,y0,X2), 
					s0(X1+1,y0,X2),s0(X1+1,y0,X2-1),s0(X1,y0,X2-1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9 
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { 
			    try = fdh2d(t0(X1-1,y0,X2),t0(X1-1,y0,X2+1),t0(X1,y0,X2+1), 
					s0(X1,y0,X2), 
					s0(X1-1,y0,X2),s0(X1-1,y0,X2+1),s0(X1,y0,X2+1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9 
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) { 
			    try = fdh2d(t0(X1-1,y0,X2),t0(X1-1,y0,X2-1),t0(X1,y0,X2-1), 
					s0(X1,y0,X2), 
					s0(X1-1,y0,X2),s0(X1-1,y0,X2-1),s0(X1,y0,X2-1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if(guess > 1.0e9){  
			  if ( X1>x0+1 && X1<x1-1 && X2>z0+1 && X2<z1-1 ) { 
			    try = fdhnf(t0(X1,y0+1,X2), 
					  t0(X1+1,y0+1,X2),t0(X1,y0+1,X2+1), 
					  t0(X1-1,y0+1,X2),t0(X1,y0+1,X2-1), 
					  s0(X1,y0,X2), 
					  s0(X1,y0+1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			}  
			  try = t0(X1,y0+1,X2) + .5*(s0(X1,y0,X2)+s0(X1,y0+1,X2)); 
			  if (try<guess)  guess = try; 
                          if ( time0[index+1]<1.e9 && X1<nx-1 )  { 
			    try = t0(X1+1,y0,X2) + .5*(s0(X1,y0,X2)+s0(X1+1,y0,X2)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-1]<1.e9 && X1>0 )  { 
			    try = t0(X1-1,y0,X2) + .5*(s0(X1,y0,X2)+s0(X1-1,y0,X2)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index+nxy]<1.e9 && X2<nz-1 )  { 
			    try = t0(X1,y0,X2+1) + .5*(s0(X1,y0,X2)+s0(X1,y0,X2+1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-nxy]<1.e9 && X2>0 )  { 
			    try = t0(X1,y0,X2-1) + .5*(s0(X1,y0,X2)+s0(X1,y0,X2-1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if (guess<time0[index]) { 
				time0[index] = guess; 
				if (fhead>headtest)  headw[3]++; 
			} 
		} 
		if(y0 == 0) dy0 = 0; 
		y0--; 
	} 
      } 
		/* BACK SIDE */ 
      for (igrow=1;igrow<=jplus;igrow++) {   
	if(dy1){ 
		ii = 0; 
		for(k=z0+1; k<=z1-1; k++){ 
			for(i=x0+1; i<=x1-1; i++){ 
				sort[ii].time = t0(i,y1-1,k); 
				sort[ii].x1 = i; 
				sort[ii].x2 = k; 
				ii++; 
			} 
		} 
		qsort((char *)sort,ii,sizeof(struct sorted),compar); 
		for(i=0;i<ii;i++){ 
			X1 = sort[i].x1; 
			X2 = sort[i].x2; 
			index = X2*nxy + y1*nx + X1; 
			lasti = X2*nxy + (y1-1)*nx + X1; 
			fhead = 0.; 
			guess = time0[index]; 
			if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9 
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) { 
			  try = fdh3d(              t0(X1,y1-1,X2), 
				      t0(X1+1,y1-1,X2),t0(X1+1,y1-1,X2+1),t0(X1,y1-1,X2+1), 
				      t0(X1+1,y1  ,X2),t0(X1+1,y1  ,X2+1),t0(X1,y1  ,X2+1), 
				      s0(X1,y1,X2), s0(X1,y1-1,X2), 
				      s0(X1+1,y1-1,X2),s0(X1+1,y1-1,X2+1),s0(X1,y1-1,X2+1), 
				      s0(X1+1,y1  ,X2),s0(X1+1,y1  ,X2+1),s0(X1,y1  ,X2+1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9 
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { 
			  try = fdh3d(              t0(X1,y1-1,X2), 
				      t0(X1-1,y1-1,X2),t0(X1-1,y1-1,X2+1),t0(X1,y1-1,X2+1), 
				      t0(X1-1,y1  ,X2),t0(X1-1,y1  ,X2+1),t0(X1,y1  ,X2+1), 
				      s0(X1,y1,X2), s0(X1,y1-1,X2), 
				      s0(X1-1,y1-1,X2),s0(X1-1,y1-1,X2+1),s0(X1,y1-1,X2+1), 
				      s0(X1-1,y1  ,X2),s0(X1-1,y1  ,X2+1),s0(X1,y1  ,X2+1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9 
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) { 
			  try = fdh3d(              t0(X1,y1-1,X2), 
				      t0(X1+1,y1-1,X2),t0(X1+1,y1-1,X2-1),t0(X1,y1-1,X2-1), 
				      t0(X1+1,y1  ,X2),t0(X1+1,y1  ,X2-1),t0(X1,y1  ,X2-1), 
				      s0(X1,y1,X2), s0(X1,y1-1,X2), 
				      s0(X1+1,y1-1,X2),s0(X1+1,y1-1,X2-1),s0(X1,y1-1,X2-1), 
				      s0(X1+1,y1  ,X2),s0(X1+1,y1  ,X2-1),s0(X1,y1  ,X2-1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9 
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) { 
			  try = fdh3d(              t0(X1,y1-1,X2), 
				      t0(X1-1,y1-1,X2),t0(X1-1,y1-1,X2-1),t0(X1,y1-1,X2-1), 
				      t0(X1-1,y1  ,X2),t0(X1-1,y1  ,X2-1),t0(X1,y1  ,X2-1), 
				      s0(X1,y1,X2), s0(X1,y1-1,X2), 
				      s0(X1-1,y1-1,X2),s0(X1-1,y1-1,X2-1),s0(X1,y1-1,X2-1), 
				      s0(X1-1,y1  ,X2),s0(X1-1,y1  ,X2-1),s0(X1,y1  ,X2-1)); 
			  if (try<guess) guess = try; 
			} 
			if(guess > 1.0e9){  
			  if(time0[index+1] < 1.e9 && X1<nx-1 && X2>z0+1 && X2<z1-1 )  { 
			      try = fdhne(t0(X1,y1-1,X2),t0(X1+1,y1-1,X2),t0(X1+1,y1,X2), 
					  t0(X1+1,y1-1,X2-1),t0(X1+1,y1-1,X2+1), 
					  s0(X1,y1,X2), 
					  s0(X1,y1-1,X2),s0(X1+1,y1-1,X2),s0(X1+1,y1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-1] < 1.e9 && X1>0 && X2>z0+1 && X2<z1-1 )  { 
			      try = fdhne(t0(X1,y1-1,X2),t0(X1-1,y1-1,X2),t0(X1-1,y1,X2), 
					  t0(X1-1,y1-1,X2-1),t0(X1-1,y1-1,X2+1), 
					  s0(X1,y1,X2), 
					  s0(X1,y1-1,X2),s0(X1-1,y1-1,X2),s0(X1-1,y1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>x0+1 && X1<x1-1 )  { 
			      try = fdhne(t0(X1,y1-1,X2),t0(X1,y1-1,X2+1),t0(X1,y1,X2+1), 
					  t0(X1-1,y1-1,X2+1),t0(X1+1,y1-1,X2+1), 
					  s0(X1,y1,X2), 
					  s0(X1,y1-1,X2),s0(X1,y1-1,X2+1),s0(X1,y1,X2+1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nxy] < 1.e9 && X2>0 && X1>x0+1 && X1<x1-1 )  { 
			      try = fdhne(t0(X1,y1-1,X2),t0(X1,y1-1,X2-1),t0(X1,y1,X2-1), 
					  t0(X1-1,y1-1,X2-1),t0(X1+1,y1-1,X2-1), 
					  s0(X1,y1,X2), 
					  s0(X1,y1-1,X2),s0(X1,y1-1,X2-1),s0(X1,y1,X2-1) ); 
			    if (try<guess)  guess = try; 
			  } 
		        }  
			  if(time0[index+1] < 1.e9 && X1<nx-1 )  { 
			    try = fdh2d(t0(X1,y1-1,X2),t0(X1+1,y1-1,X2),t0(X1+1,y1,X2), 
					  s0(X1,y1,X2), 
					  s0(X1,y1-1,X2),s0(X1+1,y1-1,X2),s0(X1+1,y1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-1] < 1.e9 && X1>0 )  { 
			    try = fdh2d(t0(X1,y1-1,X2),t0(X1-1,y1-1,X2),t0(X1-1,y1,X2), 
					  s0(X1,y1,X2), 
					  s0(X1,y1-1,X2),s0(X1-1,y1-1,X2),s0(X1-1,y1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 )  { 
			    try = fdh2d(t0(X1,y1-1,X2),t0(X1,y1-1,X2+1),t0(X1,y1,X2+1), 
					  s0(X1,y1,X2), 
					  s0(X1,y1-1,X2),s0(X1,y1-1,X2+1),s0(X1,y1,X2+1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nxy] < 1.e9 && X2>0 )  { 
			    try = fdh2d(t0(X1,y1-1,X2),t0(X1,y1-1,X2-1),t0(X1,y1,X2-1), 
					  s0(X1,y1,X2), 
					  s0(X1,y1-1,X2),s0(X1,y1-1,X2-1),s0(X1,y1,X2-1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9 
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) { 
			    try = fdh2d(t0(X1+1,y1,X2),t0(X1+1,y1,X2+1),t0(X1,y1,X2+1), 
					s0(X1,y1,X2), 
					s0(X1+1,y1,X2),s0(X1+1,y1,X2+1),s0(X1,y1,X2+1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9 
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) { 
			    try = fdh2d(t0(X1+1,y1,X2),t0(X1+1,y1,X2-1),t0(X1,y1,X2-1), 
					s0(X1,y1,X2), 
					s0(X1+1,y1,X2),s0(X1+1,y1,X2-1),s0(X1,y1,X2-1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9 
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { 
			    try = fdh2d(t0(X1-1,y1,X2),t0(X1-1,y1,X2+1),t0(X1,y1,X2+1), 
					s0(X1,y1,X2), 
					s0(X1-1,y1,X2),s0(X1-1,y1,X2+1),s0(X1,y1,X2+1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9 
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) { 
			    try = fdh2d(t0(X1-1,y1,X2),t0(X1-1,y1,X2-1),t0(X1,y1,X2-1), 
					s0(X1,y1,X2), 
					s0(X1-1,y1,X2),s0(X1-1,y1,X2-1),s0(X1,y1,X2-1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if(guess > 1.0e9){  
			  if ( X1>x0+1 && X1<x1-1 && X2>z0+1 && X2<z1-1 ) { 
			    try = fdhnf(t0(X1,y1-1,X2), 
					  t0(X1+1,y1-1,X2),t0(X1,y1-1,X2+1), 
					  t0(X1-1,y1-1,X2),t0(X1,y1-1,X2-1), 
					  s0(X1,y1,X2), 
					  s0(X1,y1-1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			}  
			  try = t0(X1,y1-1,X2) + .5*(s0(X1,y1,X2)+s0(X1,y1-1,X2)); 
			  if (try<guess)  guess = try; 
                          if ( time0[index+1]<1.e9 && X1<nx-1 )  { 
			    try = t0(X1+1,y1,X2) + .5*(s0(X1,y1,X2)+s0(X1+1,y1,X2)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-1]<1.e9 && X1>0 )  { 
			    try = t0(X1-1,y1,X2) + .5*(s0(X1,y1,X2)+s0(X1-1,y1,X2)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index+nxy]<1.e9 && X2<nz-1 )  { 
			    try = t0(X1,y1,X2+1) + .5*(s0(X1,y1,X2)+s0(X1,y1,X2+1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-nxy]<1.e9 && X2>0 )  { 
			    try = t0(X1,y1,X2-1) + .5*(s0(X1,y1,X2)+s0(X1,y1,X2-1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if (guess<time0[index]) { 
				time0[index] = guess; 
				if (fhead>headtest)  headw[4]++; 
			} 
		} 
		if(y1 == ny-1) dy1 = 0; 
		y1++; 
	} 
      } 
		/* LEFT SIDE */ 
      for (igrow=1;igrow<=iminus;igrow++) {   
	if(dx0){ 
		ii = 0; 
		for(k=z0+1; k<=z1-1; k++){ 
			for(j=y0+1; j<=y1-1; j++){ 
				sort[ii].time = t0(x0+1,j,k); 
				sort[ii].x1 = j; 
				sort[ii].x2 = k; 
				ii++; 
			} 
		} 
		qsort((char *)sort,ii,sizeof(struct sorted),compar); 
		for(i=0;i<ii;i++){ 
			X1 = sort[i].x1; 
			X2 = sort[i].x2; 
			index = X2*nxy + X1*nx + x0; 
			lasti = X2*nxy + X1*nx + (x0+1); 
			fhead = 0.; 
			guess = time0[index]; 
			if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9 
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) { 
			  try = fdh3d(              t0(x0+1,X1,X2), 
				      t0(x0+1,X1+1,X2),t0(x0+1,X1+1,X2+1),t0(x0+1,X1,X2+1), 
				      t0(x0  ,X1+1,X2),t0(x0  ,X1+1,X2+1),t0(x0  ,X1,X2+1), 
				      s0(x0,X1,X2), s0(x0+1,X1,X2), 
				      s0(x0+1,X1+1,X2),s0(x0+1,X1+1,X2+1),s0(x0+1,X1,X2+1), 
				      s0(x0  ,X1+1,X2),s0(x0  ,X1+1,X2+1),s0(x0  ,X1,X2+1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9 
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { 
			  try = fdh3d(              t0(x0+1,X1,X2), 
				      t0(x0+1,X1-1,X2),t0(x0+1,X1-1,X2+1),t0(x0+1,X1,X2+1), 
				      t0(x0  ,X1-1,X2),t0(x0  ,X1-1,X2+1),t0(x0  ,X1,X2+1), 
				      s0(x0,X1,X2), s0(x0+1,X1,X2), 
				      s0(x0+1,X1-1,X2),s0(x0+1,X1-1,X2+1),s0(x0+1,X1,X2+1), 
				      s0(x0  ,X1-1,X2),s0(x0  ,X1-1,X2+1),s0(x0  ,X1,X2+1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9 
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) { 
			  try = fdh3d(              t0(x0+1,X1,X2), 
				      t0(x0+1,X1+1,X2),t0(x0+1,X1+1,X2-1),t0(x0+1,X1,X2-1), 
				      t0(x0  ,X1+1,X2),t0(x0  ,X1+1,X2-1),t0(x0  ,X1,X2-1), 
				      s0(x0,X1,X2), s0(x0+1,X1,X2), 
				      s0(x0+1,X1+1,X2),s0(x0+1,X1+1,X2-1),s0(x0+1,X1,X2-1), 
				      s0(x0  ,X1+1,X2),s0(x0  ,X1+1,X2-1),s0(x0  ,X1,X2-1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9 
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) { 
			  try = fdh3d(              t0(x0+1,X1,X2), 
				      t0(x0+1,X1-1,X2),t0(x0+1,X1-1,X2-1),t0(x0+1,X1,X2-1), 
				      t0(x0  ,X1-1,X2),t0(x0  ,X1-1,X2-1),t0(x0  ,X1,X2-1), 
				      s0(x0,X1,X2), s0(x0+1,X1,X2), 
				      s0(x0+1,X1-1,X2),s0(x0+1,X1-1,X2-1),s0(x0+1,X1,X2-1), 
				      s0(x0  ,X1-1,X2),s0(x0  ,X1-1,X2-1),s0(x0  ,X1,X2-1)); 
			  if (try<guess) guess = try; 
			} 
			if(guess > 1.0e9){  
			  if(time0[index+nx] < 1.e9 && X1<ny-1 && X2>z0+1 && X2<z1-1 )  { 
			      try = fdhne(t0(x0+1,X1,X2),t0(x0+1,X1+1,X2),t0(x0,X1+1,X2), 
					  t0(x0+1,X1+1,X2-1),t0(x0+1,X1+1,X2+1), 
					  s0(x0,X1,X2), 
					  s0(x0+1,X1,X2),s0(x0+1,X1+1,X2),s0(x0,X1+1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nx] < 1.e9 && X1>0 && X2>z0+1 && X2<z1-1 )  { 
			      try = fdhne(t0(x0+1,X1,X2),t0(x0+1,X1-1,X2),t0(x0,X1-1,X2), 
					  t0(x0+1,X1-1,X2-1),t0(x0+1,X1-1,X2+1), 
					  s0(x0,X1,X2), 
					  s0(x0+1,X1,X2),s0(x0+1,X1-1,X2),s0(x0,X1-1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>y0+1 && X1<y1-1 )  { 
			      try = fdhne(t0(x0+1,X1,X2),t0(x0+1,X1,X2+1),t0(x0,X1,X2+1), 
					  t0(x0+1,X1-1,X2+1),t0(x0+1,X1+1,X2+1), 
					  s0(x0,X1,X2), 
					  s0(x0+1,X1,X2),s0(x0+1,X1,X2+1),s0(x0,X1,X2+1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nxy] < 1.e9 && X2>0 && X1>y0+1 && X1<y1-1 )  { 
			      try = fdhne(t0(x0+1,X1,X2),t0(x0+1,X1,X2-1),t0(x0,X1,X2-1), 
					  t0(x0+1,X1-1,X2-1),t0(x0+1,X1+1,X2-1), 
					  s0(x0,X1,X2), 
					  s0(x0+1,X1,X2),s0(x0+1,X1,X2-1),s0(x0,X1,X2-1) ); 
			    if (try<guess)  guess = try; 
			  } 
		        }  
			  if(time0[index+nx] < 1.e9 && X1<ny-1 )  { 
			    try = fdh2d(t0(x0+1,X1,X2),t0(x0+1,X1+1,X2),t0(x0,X1+1,X2), 
					  s0(x0,X1,X2), 
					  s0(x0+1,X1,X2),s0(x0+1,X1+1,X2),s0(x0,X1+1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nx] < 1.e9 && X1>0 )  { 
			    try = fdh2d(t0(x0+1,X1,X2),t0(x0+1,X1-1,X2),t0(x0,X1-1,X2), 
					  s0(x0,X1,X2), 
					  s0(x0+1,X1,X2),s0(x0+1,X1-1,X2),s0(x0,X1-1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 )  { 
			    try = fdh2d(t0(x0+1,X1,X2),t0(x0+1,X1,X2+1),t0(x0,X1,X2+1), 
					  s0(x0,X1,X2), 
					  s0(x0+1,X1,X2),s0(x0+1,X1,X2+1),s0(x0,X1,X2+1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nxy] < 1.e9 && X2>0 )  { 
			    try = fdh2d(t0(x0+1,X1,X2),t0(x0+1,X1,X2-1),t0(x0,X1,X2-1), 
					  s0(x0,X1,X2), 
					  s0(x0+1,X1,X2),s0(x0+1,X1,X2-1),s0(x0,X1,X2-1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9 
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) { 
			    try = fdh2d(t0(x0,X1+1,X2),t0(x0,X1+1,X2+1),t0(x0,X1,X2+1), 
					s0(x0,X1,X2), 
					s0(x0,X1+1,X2),s0(x0,X1+1,X2+1),s0(x0,X1,X2+1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9 
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) { 
			    try = fdh2d(t0(x0,X1+1,X2),t0(x0,X1+1,X2-1),t0(x0,X1,X2-1), 
					s0(x0,X1,X2), 
					s0(x0,X1+1,X2),s0(x0,X1+1,X2-1),s0(x0,X1,X2-1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9 
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { 
			    try = fdh2d(t0(x0,X1-1,X2),t0(x0,X1-1,X2+1),t0(x0,X1,X2+1), 
					s0(x0,X1,X2), 
					s0(x0,X1-1,X2),s0(x0,X1-1,X2+1),s0(x0,X1,X2+1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9 
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) { 
			    try = fdh2d(t0(x0,X1-1,X2),t0(x0,X1-1,X2-1),t0(x0,X1,X2-1), 
					s0(x0,X1,X2), 
					s0(x0,X1-1,X2),s0(x0,X1-1,X2-1),s0(x0,X1,X2-1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if(guess > 1.0e9){  
			  if ( X1>y0+1 && X1<y1-1 && X2>z0+1 && X2<z1-1 ) { 
			    try = fdhnf(t0(x0+1,X1,X2), 
					  t0(x0+1,X1+1,X2),t0(x0+1,X1,X2+1), 
					  t0(x0+1,X1-1,X2),t0(x0+1,X1,X2-1), 
					  s0(x0,X1,X2), 
					  s0(x0+1,X1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			}  
			  try = t0(x0+1,X1,X2) + .5*(s0(x0,X1,X2)+s0(x0+1,X1,X2)); 
			  if (try<guess)  guess = try; 
                          if ( time0[index+nx]<1.e9 && X1<ny-1 )  { 
			    try = t0(x0,X1+1,X2) + .5*(s0(x0,X1,X2)+s0(x0,X1+1,X2)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-nx]<1.e9 && X1>0 )  { 
			    try = t0(x0,X1-1,X2) + .5*(s0(x0,X1,X2)+s0(x0,X1-1,X2)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index+nxy]<1.e9 && X2<nz-1 )  { 
			    try = t0(x0,X1,X2+1) + .5*(s0(x0,X1,X2)+s0(x0,X1,X2+1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-nxy]<1.e9 && X2>0 )  { 
			    try = t0(x0,X1,X2-1) + .5*(s0(x0,X1,X2)+s0(x0,X1,X2-1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if (guess<time0[index]) { 
				time0[index] = guess; 
				if (fhead>headtest)  headw[1]++; 
			} 
		} 
		if(x0 == 0) dx0 = 0; 
		x0--; 
	} 
      } 
		/* RIGHT SIDE */ 
      for (igrow=1;igrow<=iplus;igrow++) {   
	if(dx1){ 
		ii = 0; 
		for(k=z0+1; k<=z1-1; k++){ 
			for(j=y0+1; j<=y1-1; j++){ 
				sort[ii].time = t0(x1-1,j,k); 
				sort[ii].x1 = j; 
				sort[ii].x2 = k; 
				ii++; 
			} 
		} 
		qsort((char *)sort,ii,sizeof(struct sorted),compar); 
		for(i=0;i<ii;i++){ 
			X1 = sort[i].x1; 
			X2 = sort[i].x2; 
			index = X2*nxy + X1*nx + x1; 
			lasti = X2*nxy + X1*nx + (x1-1); 
			fhead = 0.; 
			guess = time0[index]; 
			if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9 
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) { 
			  try = fdh3d(              t0(x1-1,X1,X2), 
				      t0(x1-1,X1+1,X2),t0(x1-1,X1+1,X2+1),t0(x1-1,X1,X2+1), 
				      t0(x1  ,X1+1,X2),t0(x1  ,X1+1,X2+1),t0(x1  ,X1,X2+1), 
				      s0(x1,X1,X2), s0(x1-1,X1,X2), 
				      s0(x1-1,X1+1,X2),s0(x1-1,X1+1,X2+1),s0(x1-1,X1,X2+1), 
				      s0(x1  ,X1+1,X2),s0(x1  ,X1+1,X2+1),s0(x1  ,X1,X2+1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9 
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { 
			  try = fdh3d(              t0(x1-1,X1,X2), 
				      t0(x1-1,X1-1,X2),t0(x1-1,X1-1,X2+1),t0(x1-1,X1,X2+1), 
				      t0(x1  ,X1-1,X2),t0(x1  ,X1-1,X2+1),t0(x1  ,X1,X2+1), 
				      s0(x1,X1,X2), s0(x1-1,X1,X2), 
				      s0(x1-1,X1-1,X2),s0(x1-1,X1-1,X2+1),s0(x1-1,X1,X2+1), 
				      s0(x1  ,X1-1,X2),s0(x1  ,X1-1,X2+1),s0(x1  ,X1,X2+1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9 
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) { 
			  try = fdh3d(              t0(x1-1,X1,X2), 
				      t0(x1-1,X1+1,X2),t0(x1-1,X1+1,X2-1),t0(x1-1,X1,X2-1), 
				      t0(x1  ,X1+1,X2),t0(x1  ,X1+1,X2-1),t0(x1  ,X1,X2-1), 
				      s0(x1,X1,X2), s0(x1-1,X1,X2), 
				      s0(x1-1,X1+1,X2),s0(x1-1,X1+1,X2-1),s0(x1-1,X1,X2-1), 
				      s0(x1  ,X1+1,X2),s0(x1  ,X1+1,X2-1),s0(x1  ,X1,X2-1)); 
			  if (try<guess) guess = try; 
			} 
			if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9 
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) { 
			  try = fdh3d(              t0(x1-1,X1,X2), 
				      t0(x1-1,X1-1,X2),t0(x1-1,X1-1,X2-1),t0(x1-1,X1,X2-1), 
				      t0(x1  ,X1-1,X2),t0(x1  ,X1-1,X2-1),t0(x1  ,X1,X2-1), 
				      s0(x1,X1,X2), s0(x1-1,X1,X2), 
				      s0(x1-1,X1-1,X2),s0(x1-1,X1-1,X2-1),s0(x1-1,X1,X2-1), 
				      s0(x1  ,X1-1,X2),s0(x1  ,X1-1,X2-1),s0(x1  ,X1,X2-1)); 
			  if (try<guess) guess = try; 
			} 
			if(guess > 1.0e9){  
			  if(time0[index+nx] < 1.e9 && X1<ny-1 && X2>z0+1 && X2<z1-1 )  { 
			      try = fdhne(t0(x1-1,X1,X2),t0(x1-1,X1+1,X2),t0(x1,X1+1,X2), 
					  t0(x1-1,X1+1,X2-1),t0(x1-1,X1+1,X2+1), 
					  s0(x1,X1,X2), 
					  s0(x1-1,X1,X2),s0(x1-1,X1+1,X2),s0(x1,X1+1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nx] < 1.e9 && X1>0 && X2>z0+1 && X2<z1-1 )  { 
			      try = fdhne(t0(x1-1,X1,X2),t0(x1-1,X1-1,X2),t0(x1,X1-1,X2), 
					  t0(x1-1,X1-1,X2-1),t0(x1-1,X1-1,X2+1), 
					  s0(x1,X1,X2), 
					  s0(x1-1,X1,X2),s0(x1-1,X1-1,X2),s0(x1,X1-1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>y0+1 && X1<y1-1 )  { 
			      try = fdhne(t0(x1-1,X1,X2),t0(x1-1,X1,X2+1),t0(x1,X1,X2+1), 
					  t0(x1-1,X1-1,X2+1),t0(x1-1,X1+1,X2+1), 
					  s0(x1,X1,X2), 
					  s0(x1-1,X1,X2),s0(x1-1,X1,X2+1),s0(x1,X1,X2+1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nxy] < 1.e9 && X2>0 && X1>y0+1 && X1<y1-1 )  { 
			      try = fdhne(t0(x1-1,X1,X2),t0(x1-1,X1,X2-1),t0(x1,X1,X2-1), 
					  t0(x1-1,X1-1,X2-1),t0(x1-1,X1+1,X2-1), 
					  s0(x1,X1,X2), 
					  s0(x1-1,X1,X2),s0(x1-1,X1,X2-1),s0(x1,X1,X2-1) ); 
			    if (try<guess)  guess = try; 
			  } 
		        }  
			  if(time0[index+nx] < 1.e9 && X1<ny-1 )  { 
			    try = fdh2d(t0(x1-1,X1,X2),t0(x1-1,X1+1,X2),t0(x1,X1+1,X2), 
					  s0(x1,X1,X2), 
					  s0(x1-1,X1,X2),s0(x1-1,X1+1,X2),s0(x1,X1+1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nx] < 1.e9 && X1>0 )  { 
			    try = fdh2d(t0(x1-1,X1,X2),t0(x1-1,X1-1,X2),t0(x1,X1-1,X2), 
					  s0(x1,X1,X2), 
					  s0(x1-1,X1,X2),s0(x1-1,X1-1,X2),s0(x1,X1-1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 )  { 
			    try = fdh2d(t0(x1-1,X1,X2),t0(x1-1,X1,X2+1),t0(x1,X1,X2+1), 
					  s0(x1,X1,X2), 
					  s0(x1-1,X1,X2),s0(x1-1,X1,X2+1),s0(x1,X1,X2+1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index-nxy] < 1.e9 && X2>0 )  { 
			    try = fdh2d(t0(x1-1,X1,X2),t0(x1-1,X1,X2-1),t0(x1,X1,X2-1), 
					  s0(x1,X1,X2), 
					  s0(x1-1,X1,X2),s0(x1-1,X1,X2-1),s0(x1,X1,X2-1) ); 
			    if (try<guess)  guess = try; 
			  } 
			  if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9 
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) { 
			    try = fdh2d(t0(x1,X1+1,X2),t0(x1,X1+1,X2+1),t0(x1,X1,X2+1), 
					s0(x1,X1,X2), 
					s0(x1,X1+1,X2),s0(x1,X1+1,X2+1),s0(x1,X1,X2+1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9 
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) { 
			    try = fdh2d(t0(x1,X1+1,X2),t0(x1,X1+1,X2-1),t0(x1,X1,X2-1), 
					s0(x1,X1,X2), 
					s0(x1,X1+1,X2),s0(x1,X1+1,X2-1),s0(x1,X1,X2-1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9 
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) { 
			    try = fdh2d(t0(x1,X1-1,X2),t0(x1,X1-1,X2+1),t0(x1,X1,X2+1), 
					s0(x1,X1,X2), 
					s0(x1,X1-1,X2),s0(x1,X1-1,X2+1),s0(x1,X1,X2+1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9 
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) { 
			    try = fdh2d(t0(x1,X1-1,X2),t0(x1,X1-1,X2-1),t0(x1,X1,X2-1), 
					s0(x1,X1,X2), 
					s0(x1,X1-1,X2),s0(x1,X1-1,X2-1),s0(x1,X1,X2-1) ); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if(guess > 1.0e9){  
			  if ( X1>y0+1 && X1<y1-1 && X2>z0+1 && X2<z1-1 ) { 
			    try = fdhnf(t0(x1-1,X1,X2), 
					  t0(x1-1,X1+1,X2),t0(x1-1,X1,X2+1), 
					  t0(x1-1,X1-1,X2),t0(x1-1,X1,X2-1), 
					  s0(x1,X1,X2), 
					  s0(x1-1,X1,X2) ); 
			    if (try<guess)  guess = try; 
			  } 
			}  
			  try = t0(x1-1,X1,X2) + .5*(s0(x1,X1,X2)+s0(x1-1,X1,X2)); 
			  if (try<guess)  guess = try; 
                          if ( time0[index+nx]<1.e9 && X1<ny-1 )  { 
			    try = t0(x1,X1+1,X2) + .5*(s0(x1,X1,X2)+s0(x1,X1+1,X2)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-nx]<1.e9 && X1>0 )  { 
			    try = t0(x1,X1-1,X2) + .5*(s0(x1,X1,X2)+s0(x1,X1-1,X2)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index+nxy]<1.e9 && X2<nz-1 )  { 
			    try = t0(x1,X1,X2+1) + .5*(s0(x1,X1,X2)+s0(x1,X1,X2+1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			  if ( time0[index-nxy]<1.e9 && X2>0 )  { 
			    try = t0(x1,X1,X2-1) + .5*(s0(x1,X1,X2)+s0(x1,X1,X2-1)); 
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;} 
			  } 
			if (guess<time0[index]) { 
				time0[index] = guess; 
				if (fhead>headtest)  headw[2]++; 
			} 
		} 
		if(x1 == nx-1) dx1 = 0; 
		x1++; 
	} 
      } 
 
		/* UPDATE RADIUS */ 
		radius++; 
		/* Remove by HZ
		if(radius%5 == 0) fprintf(stderr,"Completed radius = %d\n",radius); 
		*/
                if(radius == maxrad) rad0 = 0; 
 
	}	/* END BIG LOOP */ 
 
 
	/* TEST IF REVERSE PROPAGATION IS NEEDED */ 
 
	if (headw[1]==0 && headw[2]==0 && headw[3]==0 && headw[4]==0  
		     && headw[5]==0 && headw[6]==0) 
		reverse=0; 
	else { 
		head=0; 
		if (headw[1]>0) { 
		  /* HZ
			fprintf(stderr,"Head waves found on left: %d\n",headw[1]); 
			*/
			if (headw[1]>head)  { 
				head = headw[1]; 
				srcwall = 1; 
			} 
		} 
		if (headw[2]>0) { 
		  /* HZ
			fprintf(stderr,"Head waves found on right: %d\n",headw[2]);
			*/ 
			if (headw[2]>head)  { 
				head = headw[2]; 
				srcwall = 2; 
			} 
		} 
		if (headw[3]>0) { 
		  /* HZ
			fprintf(stderr,"Head waves found on front: %d\n",headw[3]); 
			*/
			if (headw[3]>head)  { 
				head = headw[3]; 
				srcwall = 3; 
			} 
		} 
		if (headw[4]>0) { 
		  /* HZ
			fprintf(stderr,"Head waves found on back: %d\n",headw[4]); 
			*/
			if (headw[4]>head)  { 
				head = headw[4]; 
				srcwall = 4; 
			} 
		} 
		if (headw[5]>0) { 
		  /* HZ
			fprintf(stderr,"Head waves found on top: %d\n",headw[5]); 
			*/
			if (headw[5]>head)  { 
				head = headw[5]; 
				srcwall = 5; 
			} 
		} 
		if (headw[6]>0) { 
		  /* HZ
			fprintf(stderr,"Head waves found on bottom: %d\n",headw[6]); 
			*/
			if (headw[6]>head)  { 
				head = headw[6]; 
				srcwall = 6; 
			} 
		} 
		if (headpref>0 && headw[headpref]>0) { 
//			fprintf(stderr,"Preference to restart on wall opposite source\n"); 
			srcwall = headpref; 
		} 
		/* SET LOCATIONS OF SIDES OF THE CUBE SO THAT CUBE IS A FACE */ 
		dx0=1; dx1=1; dy0=1; dy1=1; dz0=1; dz1=1; rad0=1; 
		radius = 1; 
		if (srcwall == 1)	{  x1=1; 
		
		//	fprintf(stderr,"RESTART at left side of model\n");
		} 
			
		else	{  x1=nx;	dx1=0;  } 
		if (srcwall == 2)	{ x0=nx-2; 
		
		//	fprintf(stderr,"RESTART at right side of model\n"); 
		} 
		else	{  x0= -1;	dx0=0;  } 
		if (srcwall == 3)	{ y1=1; 
		//	fprintf(stderr,"RESTART at front side of model\n"); 
		} 
		else	{  y1=ny;	dy1=0;  } 
		if (srcwall == 4)	{ y0=ny-2; 
		//	fprintf(stderr,"RESTART at back side of model\n"); 
		} 
		else	{  y0= -1;	dy0=0;  } 
		if (srcwall == 5)	{ z1=1; 
		//	fprintf(stderr,"RESTART at top side of model\n"); 
		} 
		else	{  z1=nz;	dz1=0;  } 
		if (srcwall == 6)	{ z0=nz-2; 
		//	fprintf(stderr,"RESTART at bottom side of model\n");
		} 
		else	{  z0= -1;	dz0=0;  } 
		/*if (reverse == 0)   
			fprintf(stderr,"WARNING:  RESTART CANCELLED by choice of 
                        input parameter 'reverse'\n"); */
	} 
	reverse--; 
 
	}	/* END BIGGER LOOP - HOLE */ 
 
 
	/* OUTPUT COMPLETED WAVEFRONT */ 

	/*
	write(tfint,time0,nxyz*4);  
	*/

	for(i=0;i<nxyz;i++) 
	  tfield[i]=time0[i]; 
	/*
	fprintf(stderr,"wavefront done \n"); 
	*/
	
	/*
	for(j=0;j<ny;j++){
	  for(i=0;i<nx;i++)
	    printf("%f ",tfield[nx*ny*0+j*nx+i]);
	  printf("\n");
	}
	*/
	/* free memory */
	   free(time0);
	   free(slow0);
	   free(sort);
	   free(wall);	  
} 
 
/* -------------------------------------------------------------------------- */ 
 
compar(a,b) 
struct sorted *a, *b; 
{ 
	if(a->time > b->time) return(1); 
	if(b->time > a->time) return(-1); 
	else return(0); 
} 
 
/* RICHARD'S EXTENDED STENCIL */ 
float fd5(t1,t2,t4,t3,t5,t6,t7,slo) 
float t1, t2, t3, t4, t5, t6, t7, slo; 
{ 
	float x, inc1; 
	double sqrt(); 
	x = 6.0*slo*slo - (t1-t2)*(t1-t2) - (t2-t3)*(t2-t3) - (t3-t1)*(t3-t1); 
	x -= (t4-t5)*(t4-t5) + (t5-t6)*(t5-t6) + (t6-t4)*(t6-t4); 
	if( x < 0 ) { 
	/*	fprintf(stderr,"Warning: x<0 in fd5: Richard: x= %f\n",x); 
		fprintf(stderr,"      slo= %f\n",slo);*/ 
		x = 0.0; 
	} 
	inc1 = sqrt(x)/1.41428; 
 	x = t7 + inc1; 
/* FOR STABILITY, ENSURE THAT NEW POINT IS LATER THAN OLD POINTS */ 
	if( x < t1 ) x = t1; 
	if( x < t2 ) x = t2; 
	if( x < t3 ) x = t3; 
	if( x < t4 ) x = t4; 
	if( x < t5 ) x = t5; 
	if( x < t6 ) x = t6; 
	if( x < t7 ) x = t7; 
	return(x); 
} 
 
/* DIFFERENT ORDER USED IN SIDES, STUPID MISTAKE */ 
float fd6(t1,t2,t3,t4,t5,t6,t7,slo) 
float t1, t2, t3, t4, t5, t6, t7, slo; 
{ 
	float x, inc1; 
	double sqrt(); 
	x = 6.0*slo*slo - (t1-t2)*(t1-t2) - (t2-t3)*(t2-t3) - (t3-t1)*(t3-t1); 
	x -= (t4-t5)*(t4-t5) + (t5-t6)*(t5-t6) + (t6-t4)*(t6-t4); 
	if( x < 0 ) { 
/*		fprintf(stderr,"Warning: x<0 in fd6: different: x= %f\n",x); 
		fprintf(stderr,"      slo= %f\n",slo); */ 
		x = 0.0; 
	} 
	inc1 = sqrt(x)/1.41428; 
 	x = t7 + inc1; 
/* FOR STABILITY, ENSURE THAT NEW POINT IS LATER THAN OLD POINTS */ 
	if( x < t1 ) x = t1; 
	if( x < t2 ) x = t2; 
	if( x < t3 ) x = t3; 
	if( x < t4 ) x = t4; 
	if( x < t5 ) x = t5; 
	if( x < t6 ) x = t6; 
	if( x < t7 ) x = t7; 
	return(x); 
}  
 
/* FIND EXACT SOLUTION ON MAIN GRID */ 
float ex0(nx,ny,nz,xs,ys,zs,index) 
int nx, ny, nz, xs, ys, zs, index; 
{ 
	int nxr, nyr, nzr; 
	float try; 
	double sqrt(); 
	nxr = ((index%(nx*ny))%nx); 
	nyr = (index%(nx*ny))/nx; 
	nzr = index/(nx*ny); 
	try = sqrt((float)((xs-nxr)*(xs-nxr) + (ys-nyr)*(ys-nyr) + (zs-nzr)*(zs-nzr))); 
	return(try); 
} 
 
/* NEW FACE STENCIL, JEV, 11-15-88 */ 
float fd7(t1,t2,t3,t4,t5,slo) 
float t1, t2, t3, t4, t5, slo; 
{ 
	float x; 
	double sqrt(); 
	x = slo*slo - 0.25*((t1-t3)*(t1-t3) + (t2-t4)*(t2-t4)); 
	if( x < 0 ) { 
	/*	fprintf(stderr,"Warning: x<0 in fd7: new face \n");*/ 
		x = 0.0; 
	} 
	x = t5 + sqrt(x); 
/* FOR STABILITY, ENSURE THAT NEW POINT IS LATER THAN OLD POINTS */ 
	if( x < t1 ) x = t1; 
	if( x < t2 ) x = t2; 
	if( x < t3 ) x = t3; 
	if( x < t4 ) x = t4; 
	return(x); 
} 
 
/* NEW EDGE STENCIL, JEV, 11-15-88 */ 
float fd8(t3,t4,t1,t2,t5,slo) 
float t1, t2, t3, t4, t5, slo; 
{ 
	float x; 
	double sqrt(); 
	x = slo*slo*2.0 - (t1-t2)*(t1-t2)*0.5 - (t3-t4)*(t3-t4); 
	if( x < 0 ) { 
	/*	fprintf(stderr,"Warning: x<0 in fd8: new edge \n");*/ 
		x = 0.0; 
	} 
	x = t5 + sqrt(x); 
/* FOR STABILITY, ENSURE THAT NEW POINT IS LATER THAN OLD POINTS */ 
	if( x < t1 ) x = t1; 
	if( x < t2 ) x = t2; 
	if( x < t3 ) x = t3; 
	if( x < t4 ) x = t4; 
	return(x); 
} 
 
/* 3D TRANSMISSION STENCIL 
   STENCIL FROM VIDALE; CONDITIONS AND OTHER OPTIONS FROM HOLE 
   JAH 11/91 */ 
float fdh3d(t1,t2,t3,t4,t5,t6,t7,ss0,s1,s2,s3,s4,s5,s6,s7) 
     float  t1,t2,t3,t4,t5,t6,t7,ss0,s1,s2,s3,s4,s5,s6,s7; 
     /* ss0 at newpoint; s1,t1 adjacent on oldface; 
	s2,t2 and s4,t4 on oldface adjacent to s1; 
	s3,t3 on oldface diametrically opposite newpoint; 
	s5,t5 on newface adjacent to newpoint AND to s2; 
	s6,t6 on newface diagonal to newpoint (adjacent to s3); 
	s7,t7 on newface adjacent to newpoint AND to s4 
	*/ 
{ 
  float x,slo; 
  double sqrt(); 
  slo = .125*(ss0+s1+s2+s3+s4+s5+s6+s7); 
  x = 6.*slo*slo - (t4-t2)*(t4-t2) - (t2-t6)*(t2-t6) - (t6-t4)*(t6-t4) 
                 - (t7-t5)*(t7-t5) - (t5-t1)*(t5-t1) - (t1-t7)*(t1-t7); 
  if (x>=0.)  { 
    x = t3 + sqrt(.5*x); 
    if ( (x<t1) || (x<t2) || (x<t4) || (x<t5) || (x<t6) || (x<t7) )   
      x = 1.e11;   /* ACAUSAL; ABORT */ 
  } 
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */ 
  return(x); 
} 
 
/* 3D STENCIL FOR NEW EDGE 
   STENCIL FROM VIDALE; CONDITIONS AND OTHER OPTIONS FROM HOLE 
   JAH 11/91 */ 
float fdhne(t1,t2,t3,t4,t5,ss0,s1,s2,s3) 
     float  t1,t2,t3,t4,t5,ss0,s1,s2,s3; 
     /* ss0 at newpoint; s1,t1 adjacent on oldface; 
	s2,t2 diagonal on oldface; s3,t3 adjacent on newface; 
	t4,t5 beside t2 on old face opposite each other */ 
{ 
  float x,slo; 
  double sqrt(); 
  slo = .25*(ss0+s1+s2+s3); 
  x = 2.*slo*slo - (t3-t1)*(t3-t1) - .5*(t5-t4)*(t5-t4); 
  if (x>=0.)  { 
    x = t2 + sqrt(x); 
    if ( (x<t1) || (x<t3) || (x<t4) || (x<t5) )     /* ACAUSAL; ABORT */ 
      x = 1.e11; 
  } 
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */ 
  return(x); 
} 
 
/* 2D TRANSMISSION STENCIL (FOR HEAD WAVES ON FACES OF GRID CELLS) 
   STENCIL FROM VIDALE (1988 2D PAPER); CONDITIONS AND OTHER OPTIONS FROM HOLE 
   JAH 11/91 */ 
float fdh2d(t1,t2,t3,ss0,s1,s2,s3) 
     float  t1,t2,t3,ss0,s1,s2,s3; 
     /* ss0 at newpoint; s1,t1 & s3,t3 adjacent; s2,t2 diagonal 
      */ 
{ 
  float x,slo; 
  double sqrt(); 
  slo = .25*(ss0+s1+s2+s3); 
  x = 2.*slo*slo - (t3-t1)*(t3-t1); 
  if (x>=0.)  { 
    x = t2 + sqrt(x); 
    if ( (x<t1) || (x<t3) )  x = 1.e11;   /* ACAUSAL; ABORT */ 
  } 
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */ 
  return(x); 
} 
 
/* 3D STENCIL FOR NEW FACE 
   STENCIL FROM VIDALE; CONDITIONS AND OTHER OPTIONS FROM HOLE 
   JAH 11/91 */ 
float fdhnf(t1,t2,t3,t4,t5,ss0,s1) 
     float  t1,t2,t3,t4,t5,ss0,s1; 
     /* ss0 at newpoint; s1,t1 adjacent on old face; 
	t2,t4 beside t1 on old face and opposite each other; 
	t3,t5 beside t1 on old face and opposite each other 
	*/ 
{ 
  float x,slo; 
  double sqrt(); 
  slo = .5*(ss0+s1); 
  x = slo*slo - .25*( (t4-t2)*(t4-t2) + (t5-t3)*(t5-t3) ); 
  if (x>=0.)  { 
    x = t1 + sqrt(x); 
    if ( (x<t2) || (x<t3) || (x<t4) || (x<t5) )     /* ACAUSAL; ABORT */ 
      x = 1.e11; 
  } 
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */ 
  return(x); 
} 
 
 
/* ==================================================================== 
   ==================================================================== */ 
 
/* copyright (c) Robert W. Clayton 
 *		 Seismological Laboratory 
 *		 Caltech 
 *		 Pasadena, CA 91125 
 * 
 * Getpar routines: 
 * 
 * Externally visable routines: 
 * 
 *		setpar(argc,argv) 
 *		getpar(name,type,valptr) 
 *		mstpar(name,type,valptr) 
 *		endpar() 
 * 
 * To get C-version: 
 *		cc -c getpar.c 
 * 
 * To get F77-version: 
 *		cp getpar.c fgetpar.c 
 *		cc -c -DFORTRAN fgetpar.c 
 *		rm fgetpar.c 
 * 
 * To get the environment processing stuff add the flag 
 *-DENVIRONMENT to each of the cc's above. 
 */ 
#include	<stdio.h> 
 
#define MAXLINE		1024	/* max length of line in par file */ 
#define MAXNAME		64	/* max length of name */ 
#define MAXVALUE	1024	/* max length of value */ 
#define MAXFILENAME	64	/* max length of par file name */ 
#define MAXVECTOR	10	/* max # of elements for unspecified vectors */ 
#define GETPAR_ERROR	100	/* exit status for getpar error */ 
#define GETPAR_STOP	101	/* exit status for STOP or mstpar */ 
#define MAXPARLEVEL	4	/* max recurrsion level for par files */ 
 
#ifdef FORTRAN 
#define GETPAR	getpar_ 
#define MSTPAR	mstpar_ 
#define ENDPAR	endpar_ 
#else 
#define GETPAR	getpar 
#define MSTPAR	mstpar 
#define ENDPAR	endpar 
#endif 
 
#define INIT	 1	/* bits for FLAGS (ext_par.argflags) */ 
#define STOP	 2 
#define LIST	 4 
#define END_PAR	 8 
#define VERBOSE	16 
 
#define LISTINC		32	/* increment size for arglist */ 
#define BUFINC		1024	/* increment size for argbuf */ 
 
struct arglist		/* structure of list set up by setpar */ 
   { 
	char *argname; 
	char *argval; 
	int hash; 
   }; 
struct ext_par		/* global variables for getpar */ 
   { 
	char *progname; 
	int argflags; 
	struct arglist *arglist; 
	struct arglist *arghead; 
	char *argbuf; 
	int nlist; 
	int nbuf; 
	int listmax; 
	int bufmax; 
	FILE *listout; 
   }	ext_par; 
 
/* abbreviations: */ 
#define AL 		struct arglist 
#define PROGNAME	ext_par.progname 
#define FLAGS		ext_par.argflags 
#define ARGLIST		ext_par.arglist 
#define ARGHEAD		ext_par.arghead 
#define ARGBUF		ext_par.argbuf 
#define NLIST		ext_par.nlist 
#define NBUF		ext_par.nbuf 
#define LISTMAX		ext_par.listmax 
#define BUFMAX		ext_par.bufmax 
#define LISTFILE	ext_par.listout 
 
#ifdef FORTRAN 
setpar_() 
#else 
setpar(ac,av)		/* set up arglist & process INPUT command */ 
int ac; char **av; 
#endif 
   { 
	register char *pl, *pn, *pv; 
	char  t, name[MAXNAME], value[MAXVALUE]; 
	FILE *file, *gp_create_dump(); 
	int i, addflags, nevlist, testav, testae; 
	struct arglist *alptr; 
#ifdef FORTRAN 
	int ac; char **av; 
	extern int xargc; extern char **xargv; 
	ac= xargc; av= xargv; 
#endif 
 
	PROGNAME= *av; 
	FLAGS= INIT; 
	LISTFILE= stderr; 
 
	ARGLIST= NULL; 
	ARGBUF = NULL; 
	NLIST= NBUF= LISTMAX= BUFMAX= 0; 
#ifdef ENVIRONMENT 
	gp_do_environment(ac,av); 
#endif 
	nevlist= NLIST; 
	while(--ac>0) 
	   { 
		av++; 
		pl= *av; 
		while(*pl == ' ' || *pl == '\t') pl++; 
		/* get name */ 
		pn= name; 
		while(*pl != '=' && *pl != '\0') *pn++ = *pl++; 
		*pn++ = '\0'; 
		/* get value */ 
		if(*pl == '=') pl++; 
		pv= value; 
		if(*pl == '"' || *pl == '\'') 
		   { 
			t= *pl++; 
			while(*pl != '\0') 
			   { 
				if(*pl == t) 
				   { 
					if(pl[-1] != '\\') break; 
					pv[-1]= t; 
					pl++; 
				   } 
				 else	*pv++ = *pl++; 
			   } 
		   } 
		 else	while(*pl) *pv++ = *pl++; 
		*pv= '\0'; 
		if(name[0] == '-') gp_add_entry("SWITCH",&name[1]); 
		 else		gp_add_entry(name,value); 
		if(strcmp("par",name)==0) /* par file */ 
			gp_do_par_file(value,1); 
	   } 
	/* do not internally call getpar before this point because 
	   ARGHEAD is not set. The search will have no stopping point */ 
	ARGHEAD= ARGLIST; 
#ifdef ENVIRONMENT 
	*value= '\0'; 
	if(GETPAR("NOENV","s",value)) ARGHEAD= ARGLIST+ nevlist; 
#endif 
	addflags= 0; 
	*value= '\0'; 
	if(GETPAR("STOP","s",value)) addflags |= STOP; 
	*value= '\0'; 
	if(GETPAR("VERBOSE","s",value)) addflags |= VERBOSE; 
	*value= '\0'; 
	if(GETPAR("LIST","s",value)) 
	   { 
		addflags |= LIST; 
		LISTFILE =gp_create_dump(value,"list"); 
	   } 
	*value= '\0'; 
	if(GETPAR("INPUT","s",value)) 
	   { 
		file =gp_create_dump(value,"list input"); 
		fprintf(file,"%s: getpar input listing\n",PROGNAME); 
		for(i=0, alptr=ARGLIST; i<NLIST; i++, alptr++) 
		   { 
			fprintf(file,"%3d: %16s = %s\n", 
				i,alptr->argname,alptr->argval); 
		   } 
		gp_close_dump(file); 
	   } 
	FLAGS |= addflags; 
   } 
 
gp_add_entry(name,value)	/* add an entry to arglist, expanding memory */ 
register char *name, *value;	/* if necessary */ 
   { 
	struct arglist *alptr; 
	int len; 
	register char *ptr; 
 
	/* check arglist memory */ 
	if(NLIST >= LISTMAX) 
	   { 
		LISTMAX += LISTINC; 
		if(ARGLIST == NULL) 
			ARGLIST= (AL *)malloc(LISTMAX * sizeof(AL)); 
		 else	ARGLIST= (AL *)realloc(ARGLIST,LISTMAX * sizeof(AL)); 
	   } 
	/* check argbuf memory */ 
	len= strlen(name) + strlen(value) + 2; /* +2 for terminating nulls */ 
	if(NBUF+len >= BUFMAX) 
	   { 
		BUFMAX += BUFINC; 
		if(ARGBUF == NULL) 
			ARGBUF= (char *)malloc(BUFMAX); 
		 else	ARGBUF= (char *)realloc(ARGBUF,BUFMAX); 
	   } 
	if(ARGBUF == NULL || ARGLIST == NULL) 
		gp_getpar_err("setpar","cannot allocate memory"); 
 
	/* add name */ 
	alptr= ARGLIST + NLIST; 
	alptr->hash= gp_compute_hash(name); 
	ptr= alptr->argname= ARGBUF + NBUF; 
	do *ptr++ = *name; while(*name++); 
 
	/* add value */ 
	NBUF += len; 
	alptr->argval= ptr; 
	do *ptr++ = *value; while(*value++); 
	NLIST++; 
   } 
#ifdef ENVIRONMENT 
gp_do_environment(ac,av) 
int ac; char **av; 
   { 
	char **ae; 
	register char *pl, *pn, *pv; 
	char name[MAXNAME], value[MAXVALUE], t; 
 
	/* The environ pointer ae, is assumed to have a specific relation 
	   to the arg pointer av. This may not be portable. */ 
	ae= av +(ac+1); 
	if(ae == NULL) return; 
 
	while(*ae != NULL) 
	   { 
		pl= *ae++; 
		while(*pl == ' ' || *pl == '\t') pl++; 
		/* get name */ 
		pn= name; 
		while(*pl != '=' && *pl != '\0') *pn++ = *pl++; 
		*pn = '\0'; 
		if(strcmp("NOENV",pn) == 0) return; 
 
		/* get value */ 
		if(*pl == '=') pl++; 
		pv= value; 
		if(*pl == '"' || *pl == '\'') 
		   { 
			t= *pl++; 
			while(*pl != '\0') 
			   { 
				if(*pl == t) 
				   { 
					if(pl[-1] != '\\') break; 
					pv[-1]= t; 
					pl++; 
				   } 
				 else	*pv++ = *pl++; 
			   } 
		   } 
		 else	while(*pl) *pv++ = *pl++; 
		*pv= '\0'; 
		gp_add_entry(name,value); 
	   } 
   } 
#endif 
 
ENDPAR()	/* free arglist & argbuf memory, & process STOP command */ 
   { 
	if(ARGLIST != NULL) free(ARGLIST); 
	if(ARGBUF  != NULL) free(ARGBUF); 
	ARGBUF=  NULL; 
	ARGLIST= NULL; 
	if(FLAGS & STOP) 
	   { 
		fprintf(stderr,"%s[endpar]: stop due to STOP in input\n", 
			PROGNAME); 
		exit(GETPAR_STOP); 
	   } 
	FLAGS= END_PAR;	/* this stops further getpar calls */ 
   } 
 
#ifdef FORTRAN 
mstpar_(name,type,val,dum1,dum2) 
int dum1, dum2;	/* dum1 & dum2 are extra args that fortran puts in */ 
#else 
mstpar(name,type,val) 
#endif 
char *name, *type; 
int *val; 
   { 
	int cnt; 
	char *typemess; 
 
	if( (cnt= GETPAR(name,type,val)) > 0) return(cnt); 
 
	/* The following line corrects a common input error */ 
	if(type[1]=='v') { type[1]= type[0]; type[0]='v'; } 
 
	switch(*type) 
	   { 
		case 'd': typemess= "an integer";	break; 
		case 'f': typemess= "a float";		break; 
		case 'F': typemess= "a double";		break; 
		case 's': typemess= "a string";		break; 
		case 'b': typemess= "a boolean";	break; 
		case 'v': switch(type[1]) 
			   { 
				case 'd': typemess= "an integer vector"; break; 
				case 'f': typemess= "a float vector"; 	 break; 
				case 'F': typemess= "a double vector";	 break; 
				default : typemess= "unknow vectorn (error)"; 
					break; 
			   } 
			  break; 
		default : typemess= "unknown (error)";	break; 
	   } 
	gp_getpar_err("mstpar","must specify value for '%s', expecting %s", 
		name,typemess); 
   } 
 
#ifdef FORTRAN 
getpar_(name,type,val,dum1,dum2) 
int dum1, dum2;	/* dum1 & dum2 are extra args that fortran puts in */ 
#else 
getpar(name,type,val) 
#endif 
char *name, *type; 
int *val; 
   { 
	register char *sptr; 
	register struct arglist *alptr; 
	register int i; 
	double atof(), *dbl; 
	float *flt; 
	int h, hno, hyes, found; 
	char line[MAXLINE], *str, *noname; 
 
	if(FLAGS & END_PAR) 
		gp_getpar_err("getpar","called after endpar"); 
	if( (FLAGS & INIT) == 0) 
		gp_getpar_err("getpar","not initialized with setpar"); 
	if(FLAGS & VERBOSE) 
		fprintf(stderr,"getpar: looking for %s\n",name); 
 
	/* The following line corrects a common input error */ 
	if(type[1]=='v') { type[1]= type[0]; type[0]='v'; } 
 
 
	if(*type == 'b') goto boolean; 
 
	h= gp_compute_hash(name); 
	found=0; 
	/* search list backwards, stopping at first find */ 
	for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--) 
	   { 
		if(alptr->hash != h) continue; 
		if(strcmp(alptr->argname,name) != 0) continue; 
		str= alptr->argval; 
		switch(*type) 
		   { 
			case 'd': 
				*val= atoi(str); 
				found=1; 
				break; 
			case 'f': 
				flt= (float *) val; 
				*flt= atof(str); 
				found=1; 
				break; 
			case 'F': 
				dbl= (double *) val; 
				*dbl= atof(str); 
				found=1; 
				break; 
			case 's': 
				sptr= (char *) val; 
				while(*str) *sptr++ = *str++; 
				*sptr= '\0'; 
				found=1; 
				break; 
			case 'v': 
				found= gp_getvector(str,type,val); 
				break; 
			default: 
				gp_getpar_err("getpar", 
					"unknown conversion type %s",type); 
				break; 
		   } 
		break; 
	   } 
	goto list; 
boolean: 
	noname= line; 
	sprintf(noname,"no%s",name); 
	hno = gp_compute_hash(noname); 
	hyes= gp_compute_hash(  name); 
	found=0; 
	/* search list backwards, stopping at first find */ 
	for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--) 
	   { 
		if(alptr->hash != hno && alptr->hash != hyes) continue; 
		if(strcmp(alptr->argname,  name)== 0) 
		   { 
			if(alptr->argval[0] == '\0') *val= 1; 
			 else *val= atol(alptr->argval); 
			found++; 
			break; 
		   } 
		if(strcmp(alptr->argname,noname)== 0) 
		   {	*val= 0; found++; break; } 
	   } 
   list: 
	if(FLAGS & LIST) 
	   { 
		switch(*type) 
		   { 
			case 'd': sprintf(line,"(int) = %d",*val); break; 
			case 'f': flt= (float *)val; 
				  sprintf(line,"(flt) = %14.6e",*flt); break; 
			case 'F': dbl= (double *)val; 
				  sprintf(line,"(dbl) = %14.6e",*dbl); break; 
			case 's': sprintf(line,"(str) = %s",val); break; 
			case 'b': sprintf(line,"(boo) = %d",*val); break; 
			case 'v': switch(type[1]) 
				   { 
					/* should list these out */ 
					case 'd': sprintf(line,"(int vec)"); 
						break; 
					case 'f': sprintf(line,"(flt vec)"); 
						break; 
					case 'F': sprintf(line,"(dbl vec)"); 
						break; 
					default : sprintf(line," vec type error"); 
						break; 
				   } 
				  break; 
			default : sprintf(line," type error"); break; 
		   } 
		fprintf(LISTFILE,"%16s (%s) %s \n",name, 
			(found ? "set":"def"),line); 
	   } 
	return(found); 
   } 
FILE *gp_create_dump(fname,filetype) 
char *fname; 
char *filetype; 
   { 
	FILE *temp; 
 
	if(*fname == '\0') return(stderr); 
	if(strcmp(fname,"stderr") == 0) return(stderr); 
	if(strcmp(fname,"stdout") == 0) return(stdout); 
	if( (temp= fopen(fname,"w")) != NULL) return(temp); 
	fprintf(stderr,"%s[setpar]: cannot create %s file %s\n", 
		PROGNAME,filetype,fname); 
	return(stderr); 
   } 
 
gp_close_dump(file) 
FILE *file; 
   { 
	if(file == stderr || file == stdout) return; 
	fclose(file); 
   } 
 
gp_compute_hash(s) 
register char *s; 
   { 
	register int h; 
	h= s[0]; 
	if(s[1]) h |= (s[1])<<8;	else return(h); 
	if(s[2]) h |= (s[2])<<16;	else return(h); 
	if(s[3]) h |= (s[3])<<24; 
	return(h); 
   } 
 
gp_do_par_file(fname,level) 
char *fname; 
int level; 
   { 
	register char *pl, *pn, *pv; 
	char t1, t2, line[MAXLINE], name[MAXNAME], value[MAXVALUE]; 
	FILE *file, *fopen(); 
 
	if(level > MAXPARLEVEL) 
		gp_getpar_err("setpar","%d (too many) recursive par file",level); 
		 
	if( (file=fopen(fname,"r"))==NULL) 
		gp_getpar_err("setpar","cannot open par file %s",fname); 
 
	while( fgets(line,MAXLINE,file) != NULL ) 
	   { 
		pl= line; 
		/* loop over entries on each line */ 
	loop:	while(*pl==' ' || *pl=='\t') pl++; 
		if(*pl=='\0'|| *pl=='\n') continue; 
		if(*pl=='#') continue; /* comments on rest of line */ 
 
		/* get name */ 
		pn= name; 
		while(*pl != '=' && *pl != '\0' && *pl != ' ' 
			&& *pl != '\t') *pn++ = *pl++; 
		*pn = '\0'; 
		if(*pl == '=') pl++; 
 
		/* get value */ 
		*value= '\0'; 
		pv= value; 
		if(*pl=='"' || *pl=='\'')	{ t1= t2= *pl++; } 
		 else				{ t1= ' '; t2= '\t'; } 
		while(*pl!=t1 && *pl!=t2 && 
			*pl!='\0' && *pl!='\n') *pv++= *pl++; 
		*pv= '\0'; 
		if(*pl=='"' || *pl=='\'') pl++; 
		gp_add_entry(name,value); 
		if(strcmp("par",name) == 0) 
			gp_do_par_file(value,level+1); 
		goto loop; 
	   } 
	fclose(file); 
   } 
 
gp_getpar_err(subname,mess,a1,a2,a3,a4) 
char *subname, *mess; 
int a1, a2, a3, a4; 
   { 
	fprintf(stderr,"\n***** ERROR in %s[%s] *****\n\t", 
		(PROGNAME == NULL ? "(unknown)" : PROGNAME),subname); 
	fprintf(stderr,mess,a1,a2,a3,a4); 
	fprintf(stderr,"\n"); 
	exit(GETPAR_ERROR); 
   } 
gp_getvector(list,type,val) 
char *list, *type; 
int *val; 
   { 
	register char *p; 
	register int index, cnt; 
	char *valptr; 
	int limit; 
	int ival, *iptr; 
	float fval, *fptr; 
	double dval, *dptr, atof(); 
 
	limit= MAXVECTOR; 
	if(type[2] == '(' || type[2] == '[') limit= atol(&type[3]); 
	if(limit <= 0) 
		gp_getpar_err("getpar","bad limit=%d specified",limit); 
	index= 0; 
	p= list; 
	while(*p != '\0'  && index < limit) 
	   { 
		cnt=1; 
	 backup: /* return to here if we find a repetition factor */ 
		while(*p == ' ' || *p == '\t') p++; 
		if(*p == '\0') return(index); 
		valptr= p; 
		while( *p != ',' && *p != '*' && *p != 'x' && *p != 'X' && 
			*p != '\0') p++; 
		if(*p == '*' || *p == 'x' || *p == 'X') 
		   { 
			cnt= atol(valptr); 
			if(cnt <= 0) 
				gp_getpar_err("getpar", 
					"bad repetition factor=%d specified", 
					 cnt); 
			if(index+cnt > limit) cnt= limit - index; 
			p++; 
			goto backup; 
		   } 
		switch(type[1]) 
		   { 
			case 'd': 
				iptr= (int *) val; 
				ival= atol(valptr); 
				while(cnt--) iptr[index++] = ival; 
				break; 
			case 'f': 
				fptr= (float *) val; 
				fval= atof(valptr); 
				while(cnt--) fptr[index++] = fval; 
				break; 
			case 'F': 
				dptr= (double *) val; 
				dval= atof(valptr); 
				while(cnt--) dptr[index++] = dval; 
				break; 
			default: 
				gp_getpar_err("getpar", 
					"bad vector type=%c specified",type[1]); 
				break; 
		   } 
		if(*p != '\0') p++; 
	   } 
	return(index); 
   } 
