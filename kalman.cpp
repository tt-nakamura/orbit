#include "Vec_DP.h"

void fdjac(const Vec_DP&, const Vec_DP&, Mat_DP&,
           void f(const Vec_DP&, Vec_DP&));

void kalman(void f(const Vec_DP&, Vec_DP&),
            void h(const Vec_DP&, Vec_DP&),
            Mat_DP& P,
            const Mat_DP& Q,
            const Mat_DP& R,
            Vec_DP& x,
            const Vec_DP& y)
// one step of kalman filtering
// input:
//   f(x,y) = map old state x to new state y
//   h(x,y) = map state x to observation y
//   P = covariance matrix of x
//   Q = covariance matrix of noise in f(x,y)
//   R = covariance matrix of noise in h(x,y)
//   y = observation data
// output:
//   P = filtered covariance matrix of x
//   x = filtered state
{
    int n(x.size()), m(y.size());
    Mat_DP F(n,n),H(m,n),A(n,n),B(n,n);
    Mat_DP K(n,m),L(n,m),M(m,m);
    Vec_DP v(n),z(m),w(n);

    f(x,v);
    fdjac(x,v,F,f);// numerical derivative
    transpose(A,F);
    mul(B,P,A);
    mul(P,F,B);
    add(P,P,Q);
    h(v,z);
    fdjac(v,z,H,h);// numerical derivative
    transpose(K,H);
    mul(L,P,K);
    mul(M,H,L);
    add(M,M,R);
    inv(M,M);
    mul(K,L,M);
    sub(z,y,z);
    mul(w,K,z);
    add(x,v,w);
    mul(B,K,H);
    ident(A);
    sub(A,A,B);
    mul(P,A,P);
}
