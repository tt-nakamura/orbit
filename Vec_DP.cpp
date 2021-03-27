#include "Vec_DP.h"
#include<cmath>
#include<iostream>

std::ostream& operator<<(std::ostream& s, const Vec_DP& a) {
    s << '[' << a[0];
    for(int i=1; i<a.size(); i++) s << ' ' << a[i];
    s << ']';
    return s;
}

double dot(const Vec_DP& a, const Vec_DP& b)
// dot product and a and b
{
    double c(0);
    int i;
    for(i=0; i<a.size(); i++) c += a[i]*b[i];
    return c;
}

double abs(const Vec_DP& a)
// absolute value of a
{
    double c(0);
    int i;
    for(i=0; i<a.size(); i++) c += a[i]*a[i];
    return sqrt(c);
}

double angle(const Vec_DP& a, const Vec_DP& b)
// angle between a and b
{
    return acos(dot(a,b)/(abs(a)*abs(b)));
}

void normalize(Vec_DP& b, const Vec_DP& a)
// b=a/|a|; &b==&a is allowed
{
    double c(abs(a));
    int i;
    for(i=0; i<a.size(); i++) b[i] = a[i]/c;
}

void cross(Vec_3D& c, const Vec_3D& a, const Vec_3D& b)
// c = a x b; &c,&a,&b need not be distinct
{
    double u,v;
    u = a[1]*b[2] - a[2]*b[1];
    v = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    c[0] = u;
    c[1] = v;
}

void mul(Mat_DP& C, const Mat_DP& A, const Mat_DP& B)
// C = AB; &C,&A,&B need not be distinct
{
    if(&C==&A) {
        Mat_DP D(A);
        mul(C,D,B);
    }
    else if(&C==&B) {
        Mat_DP D(B);
        mul(C,A,D);
    }
    else {
        int i,j,k;
        for(i=0; i<C.nrows(); i++) {
            for(j=0; j<C.ncols(); j++) {
                C[i][j] = 0;
                for(k=0; k<A.ncols(); k++)
                    C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}

void mul(Vec_DP& v, const Mat_DP& A, const Vec_DP& u)
// v = Au; &v==&u is allowed
{
    if(&v==&u) {
        Vec_DP a(u);
        mul(v,A,a);
    }
    else {
        int i,j;
        for(i=0; i<v.size(); i++) {
            v[i] = 0;
            for(j=0; j<A.ncols(); j++)
                v[i] += A[i][j]*u[j];
        }
    }
}

void mul(Vec_DP& v, const Vec_DP& u, const Mat_DP& A)
// v = uA; &v==&u is allowed
{
    if(&v==&u) {
        Vec_DP a(u);
        mul(v,a,A);
    }
    else {
        int i,j;
        for(j=0; j<v.size(); j++) {
            v[j] = 0;
            for(i=0; i<A.nrows(); i++)
                v[j] += u[i]*A[i][j];
        }
    }
}

void transpose(Mat_DP& B, const Mat_DP& A)
// B = A^T; &A==&B is allowed
{
    if(&B==&A) {
        Mat_DP C(A);
        transpose(B,C);
    }
    else {
        int i,j;
        for(i=0; i<B.nrows(); i++)
            for(j=0; j<B.ncols(); j++)
                B[i][j] = A[j][i];
    }
}

void add(Vec_DP& c, const Vec_DP& a, const Vec_DP& b)
// c=a+b; &a,&b,&c need not be distinct
{
    int i;
    for(i=0; i<c.size(); i++) c[i] = a[i] + b[i];
}

void sub(Vec_DP& c, const Vec_DP& a, const Vec_DP& b)
// c=a-b; &a,&b,&c need not be distinct
{
    int i;
    for(i=0; i<c.size(); i++) c[i] = a[i] - b[i];
}

void add(Mat_DP& C, const Mat_DP& A, const Mat_DP& B)
// C=A+B; &C,&A,&B need not be distinct
{
    int i,j;
    for(i=0; i<C.nrows(); i++)
        for(j=0; j<C.ncols(); j++)
            C[i][j] = A[i][j] + B[i][j];
}

void sub(Mat_DP& C, const Mat_DP& A, const Mat_DP& B)
// C=A+B; &C,&A,&B need not be distinct
{
    int i,j;
    for(i=0; i<C.nrows(); i++)
        for(j=0; j<C.ncols(); j++)
            C[i][j] = A[i][j] - B[i][j];
}

void negate(Vec_DP& b, const Vec_DP& a)
// b=-a; &b==&a is allowed
{
    int i;
    for(i=0; i<b.size(); i++) b[i] = -a[i];
}

void negate(Mat_DP& B, const Mat_DP& A)
// B=-A; &B==&A is allowed
{
    int i,j;
    for(i=0; i<B.nrows(); i++)
        for(j=0; j<B.ncols(); j++)
            B[i][j] = -A[i][j];
}

void ident(Mat_DP& A) {
	int i,j;
    for(i=0; i<A.nrows(); i++)
        for(j=0; j<A.ncols(); j++) A[i][j] = (i==j);
}

void inv(Mat_DP& B, const Mat_DP& A)
// B = A^{-1}; &B==&A is allowed
{
    Mat_DP C;
    if(&B!=&A) B=A;
    gaussj(B,C);
}

void RotMat(Mat_3D& R, const Vec_3D& axis, double angle)
// R = matrix of rotation about axis by angle (radian)
// rotated = rotatee * R (multiply from right)
{
    int i,j;
    double c(cos(angle)), s(sin(angle));
    Vec_3D n;
    normalize(n, axis);
    for(i=0; i<3; i++) {
        R[i][i] = c + n[i]*n[i]*(1-c);
        for(j=0; j<i; j++) R[i][j] = R[j][i] = n[i]*n[j]*(1-c);
    }
    for(i=0; i<3; i++) n[i] *= s;
    R[0][1] += n[2]; R[1][0] -= n[2];
    R[1][2] += n[0]; R[2][1] -= n[0];
    R[2][0] += n[1]; R[0][2] -= n[1];
}

void RotMat(Mat_3D& R, double phi, double theta, double psi)
// phi, theta, psi = Euler angles (radian)
// first rotate about z axis by phi and let new y axis be y'
// then rotate about y' axis by theta and let new z axis be z'
// finally rotate about z' axis by psi
{
    double cph(cos(phi)), sph(sin(phi));
    double cth(cos(theta)), sth(sin(theta));
    double cps(cos(psi)), sps(sin(psi));
    double cc(cph*cth), sc(sph*cth);
    R[0][0] = cc*cps - sph*sps;
    R[1][0] = -cc*sps - sph*cps;
    R[2][0] = cph*sth;
    R[0][1] = sc*cps + cph*sps;
    R[1][1] = -sc*sps + cph*cps;
    R[2][1] = sph*sth;
    R[0][2] = -sth*cps;
    R[1][2] = sth*sps;
    R[2][2] = cth;
}