#include<cmath>
#include "Planet.h"
#include "const.h"

double move_dt;

void move(const Vec_DP& X, Vec_DP& Y)
// input:
//   X[0:3],X[3:6] = current position and velocity
//   dt = time increment is input from global variable
// output:
//   Y[0:3],Y[3:6] = next position and velocity
{
    static double eps(1e-12);
    static int IMAX(20);
    const double& dt(move_dt);
    int i;
    double h,p,rp,rvh,c,s,e,e2,ht;
    double x,dx,z,f,g,df,dg;
    Vec_3D r,v,Z;
    for(i=0; i<3; i++) r[i] = X[i];
    for(i=0; i<3; i++) v[i] = X[i+3];
	cross(Z,r,v);
    h = abs(Z);
    p = h*h/GM;
    rp = abs(r)/p;
    rvh = dot(r,v)/h;
    c = 1/rp-1;
    s = rvh/rp;
    e = sqrt(c*c + s*s);
    e2 = 1-e*e;
    ht = dt*h/p/p;
    x = ht/rp;
    for(i=0; i<IMAX; i++) {
        z = e2*x*x;
        CSZ(c,s,z);
        dx = -(x*(x*(x*s + rvh*c) + rp*(1-z*s)) - ht)
             /(x*(x*c + rvh*(1-z*s)) + rp*(1-z*c));
        x += dx;
        if(fabs(dx) <= eps*fabs(x)) break;
    }
    f = 1 - x*x*c/rp;
    g = dt - p*p/h*pow(x,3)*s;
    for(i=0; i<3; i++) Z[i] = f*r[i] + g*v[i];
    rp = abs(Z)/p;
    dg = 1 - x*x*c/rp;
    df = (f*dg - 1)/g;
    for(i=0; i<3; i++) Y[i] = Z[i];
    for(i=0; i<3; i++) Y[i+3] = df*r[i] + dg*v[i];
}

ObsData *obs_p;

void observe(const Vec_DP& X, Vec_DP& Y)
// input:
//   X[0:3],X[3:6] = position and velocity
//   o.pos = observer's position is input from global variable
// output:
//   Y[0:2] = right ascension and declination
//   finite velocity of light is taken into account
//   o.d = distance from earth is updated
{
    static double eps(1e-12);
    static int IMAX(20),JMAX(20);
    int i,j;
    double h,p,rp,rvh,c,s,e,e2,d,dt,ht;
    double x,dx,z,f,g;
    Vec_3D r,v,Z;
    ObsData &o(*obs_p);
    for(i=0; i<3; i++) r[i] = X[i];
    for(i=0; i<3; i++) v[i] = X[i+3];
    cross(Z,r,v);
    h = abs(Z);
    p = h*h/GM;
    rp = abs(r)/p;
    rvh = dot(r,v)/h;
    c = 1/rp-1;
    s = rvh/rp;
    e = sqrt(c*c + s*s);
    e2 = 1-e*e;
    sub(Z,r,o.pos);
    d = abs(Z);
    for(i=0; i<IMAX; i++) {
		dt = -d/C_LIGHT;
        ht = dt*h/p/p;
        x = ht/rp;
        for(j=0; j<JMAX; j++) {
            z = e2*x*x;
            CSZ(c,s,z);
            dx = -(x*(x*(x*s + rvh*c) + rp*(1-z*s)) - ht)
                 /(x*(x*c + rvh*(1-z*s)) + rp*(1-z*c));
            x += dx;
            if(fabs(dx) <= eps*fabs(x)) break;
        }
        f = 1 - x*x*c/rp;
        g = dt - p*p/h*pow(x,3)*s;
        for(j=0; j<3; j++) Z[j] = f*r[j] + g*v[j];
        for(j=0; j<3; j++) Z[j] -= o.pos[j];
        o.d = d;
        if((d = abs(Z)) == o.d) break;
    }
    RAdec(Y[0],Y[1],Z);
}

Planet::Planet(const Orbit& b, double t1,
               double v_pos, double v_vel)
: t(t1), P(0.,6,6)
// construct Planet from orbit b at time t1
// v_pos = variance in position r / AU^2
// v_vel = variance in velocity v / (AU/day)^2
{
    static double eps(1.e-12);
    static int IMAX(20);
    int i;
    const Vec_3D &X(b.X), &Y(b.Y);
    const double &e(b.e), &p(b.p), &t0(b.t0);
    double nt((t-t0)*GRAV/pow(p, 1.5));
    double e1(1/(1+e)),e2(1-e*e);
    double x(nt*fabs(e2)), dx(x), u,z,s,c,y;
    for(i=0; i<IMAX; i++) {
        u = x*x;
        z = e2*u;
        CSZ(c,s,z);
        if(fabs(dx) <= eps*fabs(x)) break;
        u *= e;
        dx = (nt - x*(e1 + u*s))/(e1 + u*c);
        x += dx;
    }
    y = p*x*(1 - z*s);
    x = p*(e1 - u*c);
    for(i=0; i<3; i++) r[i] = x*X[i] + y*Y[i];
    u = GRAV*sqrt(p/(x*x + y*y));
    x = -y*u/p;
    y = (1 - z*c)*u;
    for(i=0; i<3; i++) v[i] = x*X[i] + y*Y[i];
    for(i=0; i<3; i++) P[i][i] = v_pos;
    for(i=3; i<6; i++) P[i][i] = v_vel;
}

void Planet::update(ObsData& o, double v_accel, double v_angle)
// update state (r,v) of planet by observation data o
// v_accel = variance in acceleration // / (AU/day^2)^2
// v_angle = variance in observation angles / radian^2
{
    if(o.t <= t) return;
    int i;
    double dt(o.t-t),dt2(dt*dt),dt4(dt2*dt2);
    double vr(v_accel*dt4/4),vv(v_accel*dt2);
    double vRA(v_angle/pow(cos(o.dec),2));
    Vec_DP x(6),y(2);
    Mat_DP Q(0.,6,6), R(0.,2,2);
    for(i=0; i<3; i++) x[i] = r[i];
    for(i=0; i<3; i++) x[i+3] = v[i];
    y[0] = o.RA;
    y[1] = o.dec;
    for(i=0; i<3; i++) Q[i][i] = vr;
    for(i=3; i<6; i++) Q[i][i] = vv;
    R[0][0] = vRA;
    R[1][1] = v_angle;
    move_dt = dt;
    obs_p = &o;
    kalman(move,observe,P,Q,R,x,y);// kalman filter
    for(i=0; i<3; i++) r[i] = x[i];
    for(i=0; i<3; i++) v[i] = x[i+3];
    t = o.t;
}

void Planet::get_param(double p[6]) const
// p[0] = time of perihelion passage
// p[1] = distance from sun to perihelion
// p[2] = eccentricity
// p[3] = longitude of acsending node
// p[4] = inclination of orbital plane
// p[5] = argument of perihelion
{
    Orbit b;
    b.from_rvt(r,v,t);
    b.get_param(p);
}

int determine_kf(std::vector<Orbit>& b,
                 std::vector<ObsData>& o,
                 double v_accel, double v_angle,
                 double v_pos, double v_vel)
// input:
//   o = observation data at more than 3 times
//   v_accel = variance in acceleration / (AU/day^2)^2
//   v_angle = variance in observation angles / radian^2
//   v_pos = variance in position r / AU^2
//   v_vel = variance in velocity v / (AU/day)^2
// output:
//   b = orbital elements; multiple solutions are possible
//       determined by KALMAN FILTERING of all data in o
//   o.d = distances from earth are updated
//   errors in elements are not estimated
// return:
//   number of solutions
{
    int i,j,m,n(o.size());
    Vec_DP y(2);
    m = determine(b,o);
    for(i=0; i<m; i++) {
        Planet p(b[i], o[0].t, v_pos, v_vel);
        for(j=1; j<n; j++)
            p.update(o[j], v_accel, v_angle);
        b[i].from_rvt(p.r, p.v, o[n-1].t);
    }
    return m;
}