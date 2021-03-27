#include<cmath>
#include<fstream>
#include<string>
#include<vector>
#include "ObsData.h"
#include "Observatory.h"
#include "util.h"
#include "const.h"

void ObsData::set_time(int y, int m, double d, const char *code)
// y,m,d = time of observation (UTC)
// code = observation site
// set observer's position in heliocentric ecliptic coordinates
// earth rotation is taken into account
{
    double pv[2][3], pv1[2][3];
    int dd;
    double u, dt;
    year = y;
    month = m;
    day = d;
    d = modf(d, &u);
    dd = int(u);
    iauCal2jd(y,m,dd,&t,&u);
    iauDat(y,m,dd,d,&dt);// dt = leap seconds
    t = (u+=d) + (dt + TT_TAI)/DAYSEC;// dynamical time
    iauEpv00(DJM0,t,pv,pv1);
    pos[0] = pv[0][0];
    pos[1] = pv[0][1]*cos_o + pv[0][2]*sin_o;
    pos[2] = pv[0][2]*cos_o - pv[0][1]*sin_o;
    if(code==0) return;
    Observatory *o = Observatory::lookup(code);    
    if(o==0) return;
    u = iauEra00(DJM0,u);// earth rotation angle
    iauPvtob(o->longitude, o->latitude, 0,0,0,0,u,pv);
    pos[0] += pv[0][0]/DAU;
    pos[1] += (pv[0][1]*cos_o + pv[0][2]*sin_o)/DAU;
    pos[2] += (pv[0][2]*cos_o - pv[0][1]*sin_o)/DAU;
}

void ObsData::set_RA_DEC(int h, int min, double s,
                         int deg, int am, double as)
// h,min,s = right ascension in hour,minite,second
// deg,am,as = declination in degree,arcmin,arcsec
// set planet's direction in geocentric ecliptic coordinates
{
    double u,v;
    RA = hms2h(h,min,s)*15*DEGREE;
    dec = hms2h(deg,am,as)*DEGREE;
    dir[0] = cos(dec)*cos(RA);
    u = cos(dec)*sin(RA);
    v = sin(dec);
    dir[1] = u*cos_o + v*sin_o;
    dir[2] = v*cos_o - u*sin_o;
}

int read(std::vector<ObsData>& o, const char *f)
// o = data read from file f (filename)
// return number of data
{
    int y,m,h,mm,d,am,i;
    double dd,s,as;
    char c;
    std::string str;
    std::ifstream fs(f);
    if(!fs.is_open()) error("file not exist");
    for(i=0; !fs.eof(); c=0) {
        getline(fs,str);
        // skip lines unless beginning with digit
        if(!isdigit(str[0])) continue;
        y = atoi(str.substr(15,5).c_str());
        m = atoi(str.substr(20,3).c_str());
        dd = atof(str.substr(23,9).c_str());
        h = atoi(str.substr(32,3).c_str());
        mm = atoi(str.substr(35,3).c_str());
        s = atof(str.substr(38,6).c_str());
        d = atoi(str.substr(44,4).c_str());
        am = atoi(str.substr(48,3).c_str());
        as = atof(str.substr(51,5).c_str());
        str = str.substr(77,3);
        o.resize(i+1);
        o[i].set_time(y,m,dd, str.c_str());
        o[i].set_RA_DEC(h,mm,s,d,am,as);
        i++;
    }
    return i;
}

void RAdec(double& RA, double& dec, const Vec_3D& r)
// right ascension and declination of direction to r
// r = geocentric ecliptic coordinates
{
    double u,v;
    u = r[1]*cos_o - r[2]*sin_o;
    v = r[1]*sin_o + r[2]*cos_o;
    dec = asin(v/abs(r));
    RA = atan2(u, r[0]);
    if(RA<0) RA += PI2;
}

void DateFromTime(int &y, int &m, double& d, double t)
// input: t = time (dynamical)
// output: y,m,d = calendar date (UTC)
{
    int dd;
    double s;
    t -= TT_TAI/DAYSEC;
    iauJd2cal(DJM0, t, &y, &m, &dd, &d);
    iauDat(y,m,dd,d,&s);// s = leap seconds
    s /= DAYSEC;
    if((d += dd-s) >= 1) return;
    iauJd2cal(DJM0, t-s, &y, &m, &dd, &d);
    d += dd;
}

void ObsData::get_date(int &y, int &m, double& d)
// y,m,d = calendar date (UTC)
{   y = year;
    m = month;
    d = day;
}