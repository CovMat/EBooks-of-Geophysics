#define PI 3.1415926
#define ROT 20.0

typedef struct{
        double del;
        double dist;
        double az;
        double baz;
} geome;

/*My stuff*/
geome mygrt(double,double,double,double);
void cart2sph(double,double,double,double,double,double,
              double,double*,double*,double*);
void sph2cart(double,double,double,double,double,double,
              double,double*,double*,double*);

/*Michael function translated from fortran*/
void carsph(double,double,double,double*,double*);
void locate(double,double,double,double,double*,double*);
void rot2(double,double,double,double*,double*);
void rot3(double,double,double,double*,double*);
void sphcar(double,double,double*,double*,double*);
