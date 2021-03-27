#ifndef __ObsData_h__
#define __ObsData_h__

#include "Vec_DP.h"
#include<vector>

struct ObsData {// observation data
    int year,month;
    double day;//fractional day
    double t;// time (dynamical)
    double d;// distance from earth
    double RA;// right ascension
    double dec;// declination
    Vec_3D pos;// heliocentric position of observer
    Vec_3D dir;// geocentric unit vector to planet
    void set_time(int, int, double, const char * =0);
    void get_date(int&, int&, double&);
    void set_RA_DEC(int, int, double,
                    int, int, double);
};

int read(std::vector<ObsData>& o, const char *f);
void RAdec(double&, double&, const Vec_3D&);

int iauCal2jd(int, int, int, double*, double*);
int iauDat(int, int, int, double, double*);
int iauEpv00(double, double, double[2][3], double[2][3]);
void iauPvtob(double, double, double, double,
		      double, double, double, double[2][3]);
double iauEra00(double, double);
int iauJd2cal(double dj1, double dj2,
		      int *iy, int *im, int *id, double *fd);
int iauDat(int, int, int, double, double*);

#endif // __ObsData_h__