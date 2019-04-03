c hypoDD.inc: Stores parameters that define array dimensions in hypoDD.
c             Modify to fit size of problem and available computer memory.
c Parameter Description:
c MAXEVE:   Max number of events (must be at least the size of the number 
c           of events listed in the event file)
c MAXDATA:  Max number of observations (must be at least the size of the 
c           number of observations).  
c MAXEVE0:  Max number of events used for SVD. If only LSQR is used, 
c           MAXEVE0 can be set to 2 to free up memory. 
c MAXDATA0: Max number of observations used for SVD. If only LSQR is used, 
c           MAXDATA0 can be set to 1 to free up memory. 
c MAXLAY:   Max number of model layers.
c MAXSTA:   Max number of stations.
c MAXCL:    Max number of clusters allowed. 
	integer*4 MAXEVE, MAXLAY, MAXDATA, MAXSTA, MAXEVE0, MAXDATA0, MAXCL
	integer*4 MAXGRID, MAXSEG

cc parameters for small size problems (e.g. SUN ULTRA-5, 256 MB RAM)
c      parameter(MAXEVE=   1200,    
c     &          MAXDATA=  1620000,
c     &          MAXEVE0=  30,	
c     &          MAXDATA0= 5000,	
c     &          MAXLAY=   12,		
c     &          MAXSTA=   520,		
c     &          MAXCL=    20)		

c parameters for medium size problems (e.g. : SUN ULTRA-2, 768 MB RAM)
       parameter(MAXEVE=  200,    
     &          MAXDATA=  9000,
     &          MAXEVE0=  180,	
     &          MAXDATA0= 8000,	
     &          MAXLAY=   12,		
     &          MAXSTA=   150,	
     &          MAXGRID= 10000000,	
     &          MAXSEG=  1000000,
     &          MAXCL=    50)		

cc parameters for large problems (e.g. SUN ULTRA SPARC, 1.4 GB RAM):
c      parameter(MAXEVE=   30000,    
c     &          MAXDATA=  4000000,
c     &          MAXEVE0=  2,	
c     &          MAXDATA0= 1,	
c     &          MAXLAY=   12,		
c     &          MAXSTA=   2000,		
c     &          MAXCL=    1)		
