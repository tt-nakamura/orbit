#include "util.h"
#include "Orbit.h"
#include "const.h"
#include<cmath>
#include<vector>

// q = distance from sun to perihelion / AU
double Orbit::perihel_dist() const { return p/(1+e); }

// i = angle between ecliptic and orbital plane
double Orbit::inclination() const { return acos(Z[2]); }

double Orbit::ascend_node() const
// Omega = angle between ascending node and vernal equinox
{
    double a(atan2(Z[0], -Z[1]));
    return (a<0 ? a+PI2 : a);
}

double Orbit::arg_perihel() const
// omega = angle between ascending node and perihelion
{
    double a;
    Vec_3D x,y;
    x[0]=-Z[1]; x[1]=Z[0]; x[2]=0;
    cross(y,Z,x);
    a = atan2(dot(X,y),dot(X,x));
    return (a<0 ? a+PI2 : a);
}

void Orbit::perihel_date(int &y, int &m, double& d) const
// y,m,d = time of perihelion passege (UTC)
{   DateFromTime(y,m,d,t0); }

void Orbit::from_3points(const Vec_3D r[3])
// determine conic curve from three points
// r = heliocentric ecliptic coordinates
{
    int i;
    double r1, r2, r3, u, v;
    double theta2, theta3, theta23;
    theta2 = angle(r[0], r[1]);
    theta3 = angle(r[0], r[2]);
    theta23 = theta3 - theta2;
    r1 = abs(r[0]);
    r2 = abs(r[1]);
    r3 = abs(r[2]);
    u = (sin(theta2) - sin(theta3))/r1
       + sin(theta3)/r2 - sin(theta2)/r3;
    v = (cos(theta3) - cos(theta2))/r1
       + (1 - cos(theta3))/r2 - (1 - cos(theta2))/r3;
    e = sqrt(u*u + v*v);
    p = sin(theta23)/r1 - sin(theta3)/r2 + sin(theta2)/r3;
    e /= p;
    p = (sin(theta23) - sin(theta3) + sin(theta2))/p;
    cross(Z, r[0], r[2]);
    normalize(Z,Z);
    cross(Y, Z, r[0]);
    for(i=0; i<3; i++) X[i] = u*r[0][i] + v*Y[i];
    normalize(X,X);
    cross(Y,Z,X);
}

void Orbit::from_2points(const Vec_3D r[2], double t)
// determine conic curve from two points and time
// r = heliocentric ecliptic coordinates
// t = time to go from r[0] to r[1]
{
    static double eps(1.e-12);
    static int IMAX(20);
    int i;
    double r1,r2,f,rrf,v,w;
    double u,x,c,s,du,dc,ds,z,dz(1);
    r1 = abs(r[0]);
    r2 = abs(r[1]);
    f = angle(r[0], r[1]);
    rrf = 2*sqrt(r1*r2)*cos(0.5*f);
    v = 0.5*((r1+r2)/rrf - 1);
    w = GM*t*t/pow(rrf, 3);
    z = 1+4*v/3;
    z = 16*(w/z - v*z)/(1 + v*(4-3.2*v));// init guess
    for(i=0; i<IMAX; i++) {// solve gauss eq
        CSZ(c,s,z/4);
        u = z*c/8 + v;
        du = (1 - z*s/4)/16;
        CSZ(c,s,dc,ds,z);
        if(fabs(dz) <= eps*fabs(z)) break;
        x = u*s/pow(0.5*c, 1.5);
        dz = 1+x;
        dz = (w/dz - u*dz)
            /(dz*du + 2*x*(du + u*(ds/s - 1.5*dc/c)));
        z += dz;
    }
    p = w/u*pow(r1*r2*sin(f)/t, 2)/GM;
    c = p/r1 - 1;
    s = (p/r2 - 1 - c*cos(f))/sin(f);
    e = sqrt(c*c + s*s);
    cross(Z, r[0], r[1]);
    normalize(Z,Z);
    cross(Y, Z, r[0]);
    for(i=0; i<3; i++) X[i] = c*r[0][i] + s*Y[i];
    normalize(X,X);
    cross(Y,Z,X);
}

void Orbit::from_rvt(const Vec_3D& r, const Vec_3D& v, double t)
// construct orbit from position r and velocity v at time t
{
    int i;
    double h,rp,rvh,c,s;
    cross(Z,r,v);
    h = abs(Z);
    normalize(Z,Z);
    p = h*h/GM;
    rp = abs(r)/p;
    rvh = dot(r,v)/h;
    c = 1/rp-1;
    s = -rvh/rp;
    e = sqrt(c*c + s*s);
    cross(Y,Z,r);
    for(i=0; i<3; i++) X[i] = c*r[i] + s*Y[i];
    normalize(X,X);
    cross(Y,Z,X);
    t0 = t - time(r);
}

void Orbit::position(Vec_3D& r, double t)
// r = planet's position at dynamical time t
{
    static double eps(1.e-12);
    static int IMAX(20);
    int i;
    double nt((t-t0)*GRAV/pow(p, 1.5));
    double e1(1/(1+e)),e2(1-e*e);
    double x(nt*fabs(e2)),dx(x),u,z,s,c;
    for(i=0; i<IMAX; i++) {
        u = x*x;
        z = e2*u;
        CSZ(c,s,z);
        if(fabs(dx) <= eps*fabs(x)) break;
        u *= e;
        dx = (nt - x*(e1 + u*s))/(e1 + u*c);
        x += dx;
    }
    c = p*(e1 - u*c);
    s = p*x*(1 - z*s);
    for(i=0; i<3; i++) r[i] = c*X[i] + s*Y[i];
}

double Orbit::time(const Vec_3D& r)
// time to go from perihelion to r
{
    int i;
    static int IMAX(20);
    static double eps(1.e-12);
    double e1(1/(1+e)),e2(1-e*e);
    double x,y,z,u,du,c,s;
    x = dot(r,X)/p;
    y = dot(r,Y)/p;
    u = y*y;
    for(i=0; i<IMAX; i++) {
        z = e2*u;
        CSZ(c,s,z);
        du = (e1-u*c-x)/(1-z*s);
        u += du;
        if(fabs(du) <= eps*fabs(u)) break;
    }
    x = sqrt(u);
    if(y<0) x=-x;
    return x*(e1 + e*u*s)*pow(p, 1.5)/GRAV;
}

void Orbit::direction(Vec_3D& r, ObsData& o)
// r = unit vector from earth to planet
//     at observation data o
//     in ecliptic coordinates
// update o.d = distance from earth to planet
// finite velocity of light is taken into account
{
    static int IMAX(20);
    double d(0);
    for(int i=0; i<IMAX; i++) {
        position(r, o.t - d/C_LIGHT);
        for(int j=0; j<3; j++) r[j] -= o.pos[j];
        o.d = d;
        if((d = abs(r)) == o.d) break;
    }
    normalize(r,r);
}

Orbit *orb_p;
void gauss_eq3(const Vec_DP& x, Vec_DP& f) { orb_p->gauss_eq3(x,f); }
void gauss_eq4(const Vec_DP& x, Vec_DP& f) { orb_p->gauss_eq4(x,f); }
void residual(const Vec_DP& x, Vec_DP& f) { orb_p->residual(x,f); }
double gauss_eq4a(double d2) { orb_p->gauss_eq4a(d2); }

void Orbit::gauss_eq3(const Vec_DP& x, Vec_DP& f)
// x = distance from earth to planet at time t1 and t3
// f = gauss equations to be zeroed (when inclination!=0)
// finite velocity of light is taken into account
{
    int i;
    double d1(x[0]), d3(x[1]), t[3];
    Vec_3D r[3];
    Mat_DP A(3,3);
    Vec_DP d(3);
    for(i=0; i<3; i++) r[0][i] = (o[0].pos)[i] + d1*(o[0].dir)[i];
    for(i=0; i<3; i++) r[2][i] = (o[2].pos)[i] + d3*(o[2].dir)[i];
    for(i=0; i<3; i++) A[i][0] = r[0][i];
    for(i=0; i<3; i++) A[i][2] = r[2][i];
    for(i=0; i<3; i++) A[i][1] =-(o[1].dir)[i];
    for(i=0; i<3; i++) d[i] = (o[1].pos)[i];
    gaussj(A,d);
    o[1].d = d[1];
    for(i=0; i<3; i++) r[1][i] = (o[1].pos)[i] + d[1]*(o[1].dir)[i];
    from_3points(r);
    for(i=0; i<3; i++) t[i] = time(r[i]);
    t0 = o[0].t - d1/C_LIGHT - t[0];
    f[0] = o[1].t - o[0].t - (t[1] - t[0]) - (d[1]-d1)/C_LIGHT;
    f[1] = o[2].t - o[1].t - (t[2] - t[1]) - (d3-d[1])/C_LIGHT;
}

void Orbit::gauss_eq4(const Vec_DP& x, Vec_DP& f)
// x = distance from earth to planet at time t1 and t4
// f = gauss equations to be zeroed (when inclination==0)
// finite velocity of light is taken into account
{
    double d1(x[0]), d4(x[1]);
    int i;
    Vec_3D r[2];
    for(i=0; i<3; i++) r[0][i] = (o[0].pos)[i] + d1*(o[0].dir)[i];
    for(i=0; i<3; i++) r[1][i] = (o[3].pos)[i] + d4*(o[3].dir)[i];
    from_2points(r, o[3].t - o[0].t - (d4-d1)/C_LIGHT);
    t0 = o[0].t - d1/C_LIGHT - time(r[0]);
    direction(r[0], o[1]);
    direction(r[1], o[2]);
    f[0] = r[0][0]*(o[1].dir)[1] - r[0][1]*(o[1].dir)[0];
    f[1] = r[1][0]*(o[2].dir)[1] - r[1][1]*(o[2].dir)[0];
}

int determine_by_3obs(std::vector<Orbit>& b,
                      ObsData o[3], double D_MIN)
// b = orbital elements; multiple solutions are possible
// o = observation data at three times
// D_MIN = minimum determinant of three direction vectors
// o[i].d = distances are updated by the last solution
// return number of solutions
{
    static double eps(1.e-8);
    static int IMAX(20), N(100), n(3);
    int i,j,k,l,m;
    double a1,a3,b1,b3,c1,c3;
    double R,D,E,A,B,J[3],r,y,y1,x[n],dx;
    Mat_DP U(3,3);
    Vec_DP d(2);
    Vec_3D v;
    a1 = o[2].t - o[1].t;
    a3 = o[1].t - o[0].t;
    A = a1 + a3;
    a1 /= A;
    a3 /= A;
    A *= A*GM/6;
    b1 = A*a1*(1+a1)*(1-a1);
    b3 = A*a3*(1+a3)*(1-a3);
    cross(v, o[0].dir, o[2].dir);
    D = dot(o[1].dir, v);
    E = dot(o[1].pos, o[1].dir);
    R = dot(o[1].pos, o[1].pos);
    for(i=0; i<3; i++) J[i] = dot(o[i].pos, v);
    if(fabs(D) <= D_MIN) return 0;
    A = (a1*J[0] + a3*J[2] - J[1])/D;
    B = -(b1*J[0] + b3*J[2])/D;
    if(A <= 0) return 0;
    if(B <= 0) return 0;
    dx = A/N;
    for(i=m=0; i<N && m<n; i++) {
        x[m] = (N-1-i)*dx;
        y1 = y;
        y = x[m]*(x[m] + 2*E) + R - pow(B/(A-x[m]), 2./3);
        if(i && y*y1<=0) m++;// bracketing solutions
    }
    b.clear();
    for(k=l=0; k<m; k++) {
        // estimate initial value for gauss equation
        for(i=0, y=1; i<IMAX; i++) {
            r = x[k]*(x[k] + 2*E) + R;
            if(fabs(y) <= eps*fabs(x[k])) break;
            y = pow(B/(A-x[k]), 2./3);
            y = 0.5*(y-r)/(x[k] + E - y/(3*(A-x[k])));
            x[k] += y;
        }
        r = pow(r, 1.5);
        c1 = a1 + b1/r;
        c3 = a3 + b3/r;
        for(i=0; i<3; i++) 
            v[i] = (o[1].pos)[i]
                - c1*(o[0].pos)[i]
                - c3*(o[2].pos)[i];
        for(i=0; i<3; i++)
            for(j=0; j<3; j++) U[i][j] = (o[j].dir)[i];
        gaussj(U,v);
        if((d[0] = v[0]/c1) <= 0) continue;
        if((d[1] = v[2]/c3) <= 0) continue;
        b.resize(l+1);
        b[l].o = o;
        orb_p = &b[l];
        if(newt(d, gauss_eq3) ||// solve gauss equation
           std::isnan(b[l].e)) continue;
        l++;
    }
    o[0].d = d[0];
    o[2].d = d[1];
    return l;
}

double Orbit::gauss_eq4a(double d2)
// equation to be zeroed to estimate the initial value
// for gauss equation (when inclination==0)
// d2 = distance from earth to planet at time t2
{
    double c1,c2,r,s;
    double& d3(o[2].d);
    int i;
    Mat_DP A(2,2);
    Vec_DP b(2);
    s = dot(o[1].pos, o[1].dir);
    r = dot(o[1].pos, o[1].pos);
    r += d2*(d2 + 2*s);
    r = pow(r, 1.5);
    c1 = o[2].t - o[1].t;
    c2 = o[1].t - o[0].t;
    s = c1 + c2;
    c1 /= s;
    c2 /= s;
    s *= s*GM/6;
    c1 += s*c1*(1+c1)*(1-c1)/r;
    c2 += s*c2*(1+c2)*(1-c2)/r;
    for(i=0; i<2; i++) A[i][0] = (o[0].dir)[i];
    for(i=0; i<2; i++) A[i][1] = (o[2].dir)[i];
    for(i=0; i<2; i++)
        b[i] = (o[1].pos)[i] + d2*(o[1].dir)[i]
            - c1*(o[0].pos)[i] - c2*(o[2].pos)[i];
    gaussj(A,b);
    o[0].d = b[0]/c1;
    o[2].d = b[1]/c2;
    s = dot(o[2].pos, o[2].dir);
    r = dot(o[2].pos, o[2].pos);
    r += d3*(d3 + 2*s);
    r = pow(r, 1.5);
    c1 = o[3].t - o[2].t;
    c2 = o[2].t - o[1].t;
    s = c1 + c2;
    c1 /= s;
    c2 /= s;
    s *= s*GM/6;
    c1 += s*c1*(1+c1)*(1-c1)/r;
    c2 += s*c2*(1+c2)*(1-c2)/r;
    for(i=0; i<2; i++) A[i][0] = (o[1].dir)[i];
    for(i=0; i<2; i++) A[i][1] = (o[3].dir)[i];
    for(i=0; i<2; i++)
        b[i] = (o[2].pos)[i] + d3*(o[2].dir)[i]
            - c1*(o[1].pos)[i] - c2*(o[3].pos)[i];
    gaussj(A,b);
    o[1].d = b[0]/c1;
    o[3].d = b[1]/c2;
    return o[1].d - d2;
}

int determine_by_4obs(std::vector<Orbit>& b, ObsData o[4])
// b = orbital elements; multiple solutions are possible
// o = observation data at four times
// o[i].d = distance are updated by the last solution
// return number of solutions
{
    static double eps(1.e-8), D_MIN(1.e-8);
    static int N(100), n(3);
    int i,j,m;
    double a1,a2,a3,a4,dx,y,y1,x[n];
    Vec_DP d(2),u(4);
    Mat_DP A(4,4);
    a1 = o[2].t - o[1].t;
    a3 = o[1].t - o[0].t;
    a2 = o[3].t - o[2].t;
    a4 = a1; y = a1 + a3;
    a1 /= y;
    a3 /= y; y = a2 + a4;
    a2 /= y;
    a4 /= y;
    for(i=0; i<2; i++) A[i][0] = (o[0].dir)[i];
    for(i=0; i<2; i++) A[i][1] =-(o[1].dir)[i];
    for(i=0; i<2; i++) A[i][2] = (o[2].dir)[i]*a3;
    for(i=0; i<2; i++) A[i][3] = A[i+2][0] = 0;
    for(i=0; i<2; i++) A[i+2][1] = (o[1].dir)[i]*a2;
    for(i=0; i<2; i++) A[i+2][2] =-(o[2].dir)[i];
    for(i=0; i<2; i++) A[i+2][3] = (o[3].dir)[i];
    u[0] = (o[1].pos)[0] - a1*(o[0].pos)[0] - a3*(o[2].pos)[0];
    u[1] = (o[1].pos)[1] - a1*(o[0].pos)[1] - a3*(o[2].pos)[1];
    u[2] = (o[2].pos)[0] - a2*(o[1].pos)[0] - a4*(o[3].pos)[0];
    u[3] = (o[2].pos)[1] - a2*(o[1].pos)[1] - a4*(o[3].pos)[1];
    gaussj(A,u);
    if(u[1] <= D_MIN) return 0;
    dx = 2*u[1]/N;
    Orbit orb(o);
    for(i=m=0; i<=N && m<n; i++) {
        x[m] = (N-i)*dx;
        y1 = y;
        y = orb.gauss_eq4a(x[m]);
        if(i && y*y1<=0) m++;// bracketing solutions
    }
    b.clear();
    for(i=j=0; i<m; i++) {
        b.resize(j+1);
        b[j].o = o;
        orb_p = &b[j];
        // estimate initial value for gauss equation
        zbrent(gauss_eq4a, x[i], x[i]+dx, eps);
        if((d[0] = o[0].d) <= 0) continue;
        if((d[1] = o[3].d) <= 0) continue;
        if(newt(d, gauss_eq4) ||// solve gauss equation
           std::isnan(b[j].e)) continue;
        j++;
    }
    o[0].d = d[0];
    o[3].d = d[1];
    return j;
}

int determine(std::vector<Orbit>& b, std::vector<ObsData>& o)
// input:
//   o = observation data at more than 3 times
// output:
//   b = orbital elements; multiple solutions are possible
//       determined by picking up 3 or 4 data from o
//   o.d = distances in picked-up data updated by the last solution
// return:
//   number of solutions
{
    // max and min angles between first and last data
    static double ANGL_MAX(60*DEGREE);
    static double ANGL_MIN(1*DEGREE), D_MIN(1.e-4);
    int i,j,k,l,n(o.size());
    double u[n],u1,u2;
    ObsData o1[4];
    for(k=1, u[0]=0; k<n;) {
        u[k] = u[k-1] + angle(o[k-1].dir, o[k].dir);
        if(u[k++] >= ANGL_MAX) break;
    }
    if(k<3 || u[k-1] < ANGL_MIN)
        error("need more data to determine orbit");
    u1 = 0.5*u[k-1];
    for(i=1; i<k;) if(u[i++] >= u1) break;
    if(i==k) i--;
    o1[0] = o[0];
    o1[1] = o[i-1];
    o1[2] = o[k-1];
    l = determine_by_3obs(b, o1, D_MIN);
    if(l>0 || k<=3) return l;
    u1 = u[k-1]/3;
    u2 = 2*u1;
    for(i=1; i<k;) if(u[i++] >= u1) break;
    for(j=i; j<k;) if(u[j++] >= u2) break;
    if(j==k) j--;
    while(i>=j) i--;
    o1[1] = o[i-1];
    o1[2] = o[j-1];
    o1[3] = o[k-1];
    l = determine_by_4obs(b, o1);
    return l;
}

void Orbit::residual(const Vec_DP& x, Vec_DP& y)
// x = orbital elements t0,q,e,Omega,i,omega
// y = residuals to be minimized (to be squared)
{
    int i,j,k,n(y.size()/3);
    Vec_3D r;
    set_param(&x[0]);
    for(i=k=0; i<n; i++) {
        direction(r, o[i]);
        for(j=0; j<3; j++, k++)
            y[k] = r[j] - (o[i].dir)[j];
    }
}

int determine_lsq(std::vector<Orbit>& b, std::vector<ObsData>& o)
// input:
//   o = observation data at more than 3 times
// output:
//   b = orbital elements; multiple solutions are possible
//       determined by LEAST SQUARE fit to all data in o
//   errors in orbital elements are stored in sigma[i]
//   o.d = distances updated by the last solution
// return:
//   number of solutions
{
    int i,j,m,n(o.size());
    double chisq;
    Vec_DP p(6), sigma(1., n*3);
    Mat_DP covar(6,6);// covariances of fit
    Vec<bool> ip(true,6);
    m = determine(b,o);
    for(i=0; i<m; i++) {
        b[i].get_param(&p[0]);
        b[i].o = o.data();
        orb_p = &b[i];
        // Levenberg-Marquardt minimization of residual
        mrqmin(sigma, p, ip, covar, chisq, residual);
        b[i].set_param(&p[0]);
        for(j=0; j<6; j++)
            b[i].sigma[j] = sqrt(covar[j][j]*chisq/n);
    }
    return m;
}

void Orbit::set_param(const double x[6])
// x[0] = time of perihelion passage
// x[1] = distance from sun to perihelion
// x[2] = eccentricity
// x[3] = longitude of acsending node
// x[4] = inclination of orbital plane
// x[5] = argument of perihelion
{
    int i;
    Mat_3D R;
    t0 = x[0];
    p = x[1]*(1 + x[2]);
    e = x[2];
    RotMat(R, x[3]-PI/2, x[4], x[5]+PI/2);
    for(i=0; i<3; i++) X[i] = R[0][i];
    for(i=0; i<3; i++) Y[i] = R[1][i];
    for(i=0; i<3; i++) Z[i] = R[2][i];
}

void Orbit::get_param(double p[6]) const {
    p[0] = t0;
    p[1] = perihel_dist();
    p[2] = e;
    p[3] = ascend_node();
    p[4] = inclination();
    p[5] = arg_perihel();
}

void Orbit::get_param(double p[6], double err[6]) const {
    int i;
    get_param(p);
    for(i=0; i<6; i++) err[i] = sigma[i];
}