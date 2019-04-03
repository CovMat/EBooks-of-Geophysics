#include        <stdio.h> 
#include        <math.h> 
#include	"geometry.h" 

int sph2car(float del,float az,float dep,float *x,float *y,float *z); 
int car2sph(float x,float y,float z,float *del,float *az,float *dep);
void fpos(float del,float az,float wlat, float wlon,  
          float *lat,float *lon); 
void locate(double oldlat,double oldlon,double angd,double az,
	    double *newlat,double *newlon);
float get_az(float,float);
geome mygrt(double elt,double eln,double slt,double sln);
void rot2(double oldlat,double oldlon,double angle,
	  double *newlat,double *newlon);
void rot3(double oldlat,double oldlon,double angle,
	  double *newlat,double *newlon);
void carsph(double x1,double x2,double x3,double *lat,double *lon);
void sphcar(double lat,double lon,double *x1,double *x2,double *x3);
void rotate(float x, float y, float theta, float *xr, float *yr);


/* For linux GNU gcc */
/*
void sph2car_ft__(float *lat_ft, float *lon_ft, float *dep_ft,  
                float *wlat_ft,float *wlon_ft, 
                float *x_ft, float *y_ft, float *z_ft, float *theta_ft)
*/

/* For Intel compiler */ 
void sph2car_ft_(float *lat_ft, float *lon_ft, float *dep_ft,  
                float *wlat_ft,float *wlon_ft, 
                float *x_ft, float *y_ft, float *z_ft, float *theta_ft)

{
	geome dat; 
	double lat,lon,wlat, wlon;
        float x, y, z, dep; 
	float del, az;
	float xr,yr;
	float theta;
	 
	wlat=(double)(*wlat_ft); 
	wlon=(double)(*wlon_ft); 
	dep=(double)(*dep_ft); 
	lat=(double)(*lat_ft); 
	lon=(double)(*lon_ft);
	theta=*theta_ft;
	dat=mygrt(wlat, wlon, lat, lon);
	del=dat.del;
	az=dat.az;
	sph2car(del, az, dep, &x, &y, &z); 
	
	/* rotate */
	rotate(x,y,theta,&xr,&yr);
	*x_ft=xr; 
	*y_ft=yr; 
	*z_ft=z; 
	 
} 
 
/* For Linux GNU GCC */
/*
void car2sph_ft__(float *x_ft, float *y_ft, float *z_ft, 
                float *wlat_ft,float *wlon_ft, 
                float *lat_ft, float *lon_ft, float *dep_ft, 
                float *theta_ft) 
*/

/*For Sun compiler */          
/*
void car2sph_ft_(float *x_ft, float *y_ft, float *z_ft, 
                float *wlat_ft,float *wlon_ft, 
                float *lat_ft, float *lon_ft, float *dep_ft, 
                float *theta_ft)
*/

/*For Absoft compiler */          
void car2sph_ft_(float *x_ft, float *y_ft, float *z_ft, 
                float *wlat_ft,float *wlon_ft, 
                float *lat_ft, float *lon_ft, float *dep_ft, 
                float *theta_ft)

{ 
	float lat, lon, wlat, wlon, x, y, z; 
	float del, az, dep; 
	float xr, yr;
	float theta;

	wlat=*wlat_ft; 
	wlon=*wlon_ft; 
	x=*x_ft; 
	y=*y_ft; 
	z=*z_ft; 
	theta=*theta_ft;
	/* rotate first */
	theta=-theta;
	rotate(x,y,theta,&xr, &yr);

	car2sph(xr,yr,z,&del,&az,&dep); 
	fpos(del,az,wlat, wlon, &lat,&lon); 
	 
	*lat_ft=lat; 
	*lon_ft=lon; 
	*dep_ft=dep; 
	 
}

void rotate(float x, float y, float theta, float *xr, float *yr)
{
  *xr=x*cos(theta)-y*sin(theta);
  *yr=x*sin(theta)+y*cos(theta);
}

/*
* Returns PREM P-velocity at a particular depth
*/
float prem_p_(float *dep_ft)
{
        float dep;
        float r,x;

	dep=*dep_ft;

        r=6371.0-dep;
        x=r/6371.0;

        if(r<1221.15){
                return(11.2622-6.3640*x*x);
        } else if((r>=1221.5)&&(r<3480.0)){
                return(11.0487-4.0362*x+4.8023*x*x-13.5732*x*x*x);
        } else if((r>=3480.0)&&(r<3630.1)){
                return(15.3891-5.3181*x+5.5242*x*x-2.5514*x*x*x);
        } else if((r>=3630.1)&&(r<5600.1)){
                return(24.9520-40.4673*x+51.4832*x*x-26.6419*x*x*x);
        } else if((r>=5600.1)&&(r<5701.1)){
                return(29.2766-23.6027*x+5.5242*x*x-2.5514*x*x*x);
        } else if((r>=5701.1)&&(r<5771.1)){
                return(19.0957-9.8672*x);
        } else if((r>=5771.1)&&(r<5971.1)){
                return(39.7027-32.6166*x);
        } else if((r>=5971.1)&&(r<6151.1)){
                return(20.3926-12.2569*x);
        } else if((r>=6151.1)&&(r<6346.7)){
                return(.8317+7.218*x);
        } else if((r>=6346.7)&&(r<6356.1)){
                return(6.8);
        } else if((r>=6356.1)&&(r<=6371.0)){
                return(5.8);
        } else {
                return(-1.0);
        }
}

 
/* 
* Converts spherical to Cartesian coords, center point 
* is  del=0.0, dep=0.0 (az anything) -> x=0.0, y=0.0, z=0.0 
*/ 
int sph2car(float del,float az,float dep,float *x,float *y,float *z) 
{ 
        double d1,d2,d3,sfac,dist,ddel,daz,ddep; 
 
 
	ddel=(double)del; 
	daz=(double)az; 
	ddep=(double)dep; 
 
        sfac=180.0/3.1415926; 
 
        d1=6371.0; 
        d3=d1-ddep; 
        d2=sqrt(d1*d1+d3*d3-2.0*d1*d3*cos(ddel/sfac)); 
 
        *z=(float)(0.5*(d1*d1+d2*d2-d3*d3)/d1); 
 
        dist=sqrt(d2*d2-(double)((*z)*(*z))); 
 
        *x=(float)dist*sin(daz/sfac); 
        *y=(float)dist*cos(daz/sfac); 
 
        return(0); 
} 
 
/* 
* Converts cartesian to spherical coords, center point 
* is at x=0.0, y=0.0, z=0.0 --> del=0.0, dep=0.0 (az anything) 
*/ 
 
int car2sph(float x,float y,float z,float *del,float *az,float *dep) 
{ 
        double d1,d2,d3,sfac,dx,dy,dz; 
 
        sfac=180.0/3.1415926; 
	dx=(double)x; 
	dy=(double)y; 
	dz=(double)z; 
 
        *az=get_az(x,y); 
 
        d1=6371.0; 
        d2=sqrt(dx*dx+dy*dy+dz*dz); 
        d3=sqrt(dx*dx+dy*dy+(dz-6371.0)*(dz-6371.0)); 
        *del=(float)(sfac*acos((d2*d2-d1*d1-d3*d3)/(-2.0*d1*d3))); 
 
        *dep=(float)d1-d3; 
 
        return(0); 
} 
 
void fpos(float del,float az,float wlat, float wlon,  
          float *lat,float *lon) 
{ 
	float f1; 
	double lt,ln; 
 
	f1=3.1415926/180.0; 
 
        locate(f1*wlat,f1*wlon,f1*del,f1*az,&lt,&ln); 
 
	*lat=(float)lt/f1; 
	*lon=(float)ln/f1; 
}

void locate(double oldlat,double oldlon,double angd,double az,
	    double *newlat,double *newlon)
{
	double tlat1,tlon1;

	rot2( PI / 2.-angd,PI-az, PI / 2.-oldlat,&tlat1,&tlon1);
	rot3(tlat1,tlon1,oldlon,newlat,newlon);
}

/*
* Returns azimuth with respect to North assuming that
* positive y is North and positive x is East.
*/
float get_az(float x,float y)
{
        double sfac,az,dx,dy;

	dx=(double)x;
	dy=(double)y;
        sfac=180.0/3.1415926;

        if(fabs(dx)<1.0e-10){
            if(dy>0.0){
                return(0.0);
            } else {
                return(180.0);
            }
        }

        if(dx>=0.0){
            return(90.0-(float)(sfac*atan(dy/dx)));
        } else {
            return(270.0-(float)(sfac*atan(dy/dx)));
        }
}

/*This function modeled after fortran routine in delaz.f.
The delaz.f that is in fprog, not the mild one. Enter
latitude [-90,90] and longitude [-180,180]. */
geome mygrt(double elt,double eln,double slt,double sln)
{
        double slat,slon,elat,elon;
        double a,b,c,a1,b1,c1;
        double tmp1,tmp2,tmp2a,tmp3,z,cd,bz;
        geome info;

/*Go to radians*/
        slat=slt*PI/180.0; slon=sln*PI/180.0;
        elat=elt*PI/180.0; elon=eln*PI/180.0;

/*Correct for ellipticity*/
        slat=atan(.996647*tan(slat)); elat=atan(.996647*tan(elat));

/*Got to colatitudes*/
        slat=PI/2.0-slat; elat=PI/2.0-elat;

/*Make all longitudes postive*/
        if(slon<0.0){ slon+=2.0*PI; }
        if(elon<0.0){ elon+=2.0*PI; }

/*compute direction cosines*/
        a=sin(elat)*cos(elon); b=sin(elat)*sin(elon); c=cos(elat);
        a1=sin(slat)*cos(slon); b1=sin(slat)*sin(slon); c1=cos(slat);

        cd=a*a1+b*b1+c*c1;
/*Make sure acos won't barf*/
        if(cd>1.0) { cd=1.0; }
        if(cd<-1.0) { cd=-1.0; }
        info.del=acos(cd)*180.0/PI;
        info.dist=info.del*PI*6371.0/180.0;

        tmp1=cos(elon)*cos(slon)+sin(elon)*sin(slon);
        tmp2a=1.0-cd*cd;
        if(tmp2a<=0.0){
	    tmp2=0.0;
	    tmp3=1.0;
	} else {
            tmp2=sqrt(tmp2a);
            tmp3=(sin(elat)*cos(slat)-cos(elat)*sin(slat)*tmp1)/tmp2;
	}
/*Make sure acos won't barf*/
        if(tmp3>1.0) { tmp3=1.0; }
        if(tmp3<-1.0) { tmp3=-1.0; }
        z=acos(tmp3);

/*This test gets correct orientation for az. */
        if((sin(slon)*cos(elon)-cos(slon)*sin(elon))<0.0){
            z=2.0*PI-z;
        }
        info.az=180.0*z/PI;

        tmp1=cos(slon)*cos(elon)+sin(slon)*sin(elon);
        tmp2a=1.0-cd*cd;
        if(tmp2a<=0.0){
	    tmp2=0.0;
	    tmp3=1.0;
	} else {
            tmp2=sqrt(tmp2a);
            tmp3=(sin(slat)*cos(elat)-cos(slat)*sin(elat)*tmp1)/tmp2;
/*            tmp3=(sin(elat)*cos(slat)-cos(elat)*sin(slat)*tmp1)/tmp2;*/
	}
/*Make sure acos won't barf*/
        if(tmp3>1.0) { tmp3=1.0; }
        if(tmp3<-1.0) { tmp3=-1.0; }
        bz=acos(tmp3);
/*This test gets correct orientation for baz. */
        if((sin(elon)*cos(slon)-cos(elon)*sin(slon))<0.0){
           bz=2.0*PI-bz;
        }

        info.baz=180.0*bz/PI;

        return(info);
}

void rot2(double oldlat,double oldlon,double angle,
	  double *newlat,double *newlon)
{
	double x1,x2,x3,temp;

	sphcar(oldlat,oldlon,&x1,&x2,&x3);
	temp = x1;
	x1 = x1 * cos(angle) + x3 * sin(angle);
	x3 = x3 * cos(angle) - temp * sin(angle);
	carsph(x1,x2,x3,newlat,newlon);
}

void rot3(double oldlat,double oldlon,double angle,
	  double *newlat,double *newlon)
{

	*newlat=oldlat;
	*newlon=oldlon+angle;
	if (*newlon>2.0*PI){
             *newlon-=2.0*PI;
	}else if(*newlon<-2.0*PI){
             *newlon+=2.0*PI;
	}
}

void sphcar(double lat,double lon,double *x1,double *x2,double *x3)
{
	*x1=cos(lat)*cos(lon);
	*x2=cos(lat)*sin(lon);
	*x3=sin(lat);
}

void carsph(double x1,double x2,double x3,double *lat,double *lon)
{
	if(x3>1.0){ x3=1.0; }
	if(x3<-1.0){ x3=-1.0; }
	*lat=asin(x3);
	if(x1==0.0){
	  if(x2>0.0){ *lon=PI/2.0 ; }
	  if(x2<0.0){ *lon=3.0*PI/2.0; }
	} else {
	  *lon=atan(x2/x1);
	  if(x1<0.0){ *lon+=PI; }
	  if(x1>0.0&&x2<0.0){ *lon+=2.0*PI; }
	}
}

