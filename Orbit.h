#ifndef __Orbit_h__
#define __Orbit_h__

#include "ObsData.h"

struct Orbit {
    double p;// semi-latus rectum
    double e;// eccentricity
    double t0;// time of perihelion passage
    Vec_3D Z;// unit vec normal to orbital plane
    Vec_3D X;// unit vec from sun to perihelion
    Vec_3D Y;// z cross x
    ObsData *o;// pointer to observation data
    double sigma[6];// error in t0,q,e,Omega,i,omega
    void from_rvt(const Vec_3D&, const Vec_3D&, double);
    void from_3points(const Vec_3D[3]);
    void from_2points(const Vec_3D[2], double);
    double time(const Vec_3D&);
    void position(Vec_3D&, double);
    void direction(Vec_3D&, ObsData&);
    void gauss_eq3(const Vec_DP&, Vec_DP&);
    void gauss_eq4(const Vec_DP&, Vec_DP&);
    void residual(const Vec_DP&, Vec_DP&);
    double gauss_eq4a(double);
    double inclination() const;
    double ascend_node() const;
    double arg_perihel() const;
    double perihel_dist() const;
    void perihel_date(int&, int&, double&) const;
    void set_param(const double[6]);
    void get_param(double[6]) const;
    void get_param(double[6], double[6]) const;
    Orbit(ObsData *o1) : o(o1) {};
    Orbit() {};
};

int determine_by_3obs(std::vector<Orbit>&, ObsData[3], double=0);
int determine_by_4obs(std::vector<Orbit>&, ObsData[4]);
int determine(std::vector<Orbit>&, std::vector<ObsData>&);
int determine_lsq(std::vector<Orbit>&, std::vector<ObsData>&);

void CSZ(double&, double&, double);
void CSZ(double&, double&, double&, double&, double);

void DateFromTime(int&, int&, double&, double);

void gaussj(Mat_DP&, Vec_DP&);
bool newt(Vec_DP&, void(const Vec_DP&, Vec_DP&));
double zbrent(double(const double), const double,
		      const double, const double);
void mrqmin(const Vec_DP&, Vec_DP&, const Vec<bool>&,
		    Mat_DP&, double&, void(const Vec_DP&, Vec_DP&));

#endif // __Orbit_h__