#ifndef __Planet_h__
#define __Planet_h__

#include "Orbit.h"

struct Planet {
    Vec_3D r;// position / AU
    Vec_3D v;// velocity / AU/day
    Mat_DP P;// covariance of r,v
    double t;// time (dynamical) / sec
    Planet(const Orbit&, double, double=0, double=0);
    void update(ObsData&, double, double);
    void get_param(double[6]) const;
};

int determine_kf(std::vector<Orbit>&,
                 std::vector<ObsData>&,
                 double, double, double=0, double=0);

void CSZ(double&, double&, double);
void kalman(void f(const Vec_DP&, Vec_DP&),
            void h(const Vec_DP&, Vec_DP&),
            Mat_DP&, const Mat_DP&, const Mat_DP&,
            Vec_DP&, const Vec_DP&);

#endif // __Planet_h__