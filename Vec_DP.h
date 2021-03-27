#ifndef __Vec_DP_h__
#define __Vec_DP_h__

#include "Vec.h"
#include<iostream>

struct Vec_3D : Vec_DP {
    Vec_3D() : Vec_DP(3) {};
    Vec_3D(const double *a) : Vec_DP(a,3) {};
    Vec_3D(const double& a) : Vec_DP(a,3) {};
};

struct Mat_3D : Mat_DP {
    Mat_3D() : Mat_DP(3,3) {};
    Mat_3D(const double *a) : Mat_DP(a,3,3) {};
    Mat_3D(const double& a) : Mat_DP(a,3,3) {};
};

double dot(const Vec_DP&, const Vec_DP&);
double abs(const Vec_DP&);
double angle(const Vec_DP&, const Vec_DP&);

std::ostream& operator<<(std::ostream&, const Vec_DP&);

void normalize(Vec_DP&, const Vec_DP&);
void cross(Vec_3D&, const Vec_3D&, const Vec_3D&);
void mul(Mat_DP&, const Mat_DP&, const Mat_DP&);
void mul(Vec_DP&, const Mat_DP&, const Vec_DP&);
void mul(Vec_DP&, const Vec_DP&, const Mat_DP&);
void add(Vec_DP&, const Vec_DP&, const Vec_DP&);
void sub(Vec_DP&, const Vec_DP&, const Vec_DP&);
void add(Mat_DP&, const Mat_DP&, const Mat_DP&);
void sub(Mat_DP&, const Mat_DP&, const Mat_DP&);
void negate(Vec_DP&, const Vec_DP&);
void negate(Mat_DP&, const Mat_DP&);
void transpose(Mat_DP&, const Mat_DP&);
void ident(Mat_DP&);
void inv(Mat_DP&, const Mat_DP&);
void RotMat(Mat_3D&, const Vec_3D&, double);
void RotMat(Mat_3D&, double, double, double);

void gaussj(Mat_DP&, Mat_DP&);

#endif // __Vec_DP_h__
