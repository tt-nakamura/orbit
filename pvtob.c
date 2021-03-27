#include "math.h"

#define D2PI (6.283185307179586476925287)/* 2Pi */
#define DAYSEC (86400.0)/* Seconds per day. */
#define WGS84 1/* Reference ellipsoids */
#define GRS80 2
#define WGS72 3

int iauEform ( int n, double *a, double *f )
{
	
	/* Look up a and f for the specified reference ellipsoid. */
	switch ( n ) {
			
		case WGS84:
			*a = 6378137.0;
			*f = 1.0 / 298.257223563;
			break;
			
		case GRS80:
			*a = 6378137.0;
			*f = 1.0 / 298.257222101;
			break;
			
		case WGS72:
			*a = 6378135.0;
			*f = 1.0 / 298.26;
			break;
			
		default:
			
			/* Invalid identifier. */
			*a = 0.0;
			*f = 0.0;
			return -1;
			
	}
	
	/* OK status. */
	return 0;
}	

int iauGd2gce ( double a, double f, double elong, double phi,
	       double height, double xyz[3] )
{
	double sp, cp, w, d, ac, as, r;
	
	
	/* Functions of geodetic latitude. */
	sp = sin(phi);
	cp = cos(phi);
	w = 1.0 - f;
	w = w * w;
	d = cp*cp + w*sp*sp;
	if ( d <= 0.0 ) return -1;
	ac = a / sqrt(d);
	as = w * ac;
	
	/* Geocentric vector. */
	r = (ac + height) * cp;
	xyz[0] = r * cos(elong);
	xyz[1] = r * sin(elong);
	xyz[2] = (as + height) * sp;
	
	/* Success. */
	return 0;
}	

void iauZp(double p[3])
{
	p[0] = 0.0;
	p[1] = 0.0;
	p[2] = 0.0;
	
	return;
}

int iauGd2gc ( int n, double elong, double phi, double height,
	      double xyz[3] )
{
	int j;
	double a, f;
	
	
	/* Obtain reference ellipsoid parameters. */
	j = iauEform ( n, &a, &f );
	
	/* If OK, transform longitude, geodetic latitude, height to x,y,z. */
	if ( j == 0 ) {
		j = iauGd2gce ( a, f, elong, phi, height, xyz );
		if ( j != 0 ) j = -2;
	}
	
	/* Deal with any errors. */
	if ( j != 0 ) iauZp ( xyz );
	
	/* Return the status. */
	return j;
}

void iauIr(double r[3][3])
{
	r[0][0] = 1.0;
	r[0][1] = 0.0;
	r[0][2] = 0.0;
	r[1][0] = 0.0;
	r[1][1] = 1.0;
	r[1][2] = 0.0;
	r[2][0] = 0.0;
	r[2][1] = 0.0;
	r[2][2] = 1.0;
	
	return;
}

void iauRx(double phi, double r[3][3])
{
	double s, c, a10, a11, a12, a20, a21, a22;
	
	
	s = sin(phi);
	c = cos(phi);
	
	a10 =   c*r[1][0] + s*r[2][0];
	a11 =   c*r[1][1] + s*r[2][1];
	a12 =   c*r[1][2] + s*r[2][2];
	a20 = - s*r[1][0] + c*r[2][0];
	a21 = - s*r[1][1] + c*r[2][1];
	a22 = - s*r[1][2] + c*r[2][2];
	
	r[1][0] = a10;
	r[1][1] = a11;
	r[1][2] = a12;
	r[2][0] = a20;
	r[2][1] = a21;
	r[2][2] = a22;
	
	return;
}

void iauRy(double theta, double r[3][3])
{
	double s, c, a00, a01, a02, a20, a21, a22;
	
	
	s = sin(theta);
	c = cos(theta);
	
	a00 = c*r[0][0] - s*r[2][0];
	a01 = c*r[0][1] - s*r[2][1];
	a02 = c*r[0][2] - s*r[2][2];
	a20 = s*r[0][0] + c*r[2][0];
	a21 = s*r[0][1] + c*r[2][1];
	a22 = s*r[0][2] + c*r[2][2];
	
	r[0][0] = a00;
	r[0][1] = a01;
	r[0][2] = a02;
	r[2][0] = a20;
	r[2][1] = a21;
	r[2][2] = a22;
	
	return;
}

void iauRz(double psi, double r[3][3])
{
	double s, c, a00, a01, a02, a10, a11, a12;
	
	
	s = sin(psi);
	c = cos(psi);
	
	a00 =   c*r[0][0] + s*r[1][0];
	a01 =   c*r[0][1] + s*r[1][1];
	a02 =   c*r[0][2] + s*r[1][2];
	a10 = - s*r[0][0] + c*r[1][0];
	a11 = - s*r[0][1] + c*r[1][1];
	a12 = - s*r[0][2] + c*r[1][2];
	
	r[0][0] = a00;
	r[0][1] = a01;
	r[0][2] = a02;
	r[1][0] = a10;
	r[1][1] = a11;
	r[1][2] = a12;
	
	return;
}	

void iauPom00(double xp, double yp, double sp, double rpom[3][3])
{
	
	/* Construct the matrix. */
	iauIr(rpom);
	iauRz(sp, rpom);
	iauRy(-xp, rpom);
	iauRx(-yp, rpom);
	
	return;
}

void iauCp(double p[3], double c[3])
{
	c[0] = p[0];
	c[1] = p[1];
	c[2] = p[2];
	
	return;
}

void iauCr(double r[3][3], double c[3][3])
{
	iauCp(r[0], c[0]);
	iauCp(r[1], c[1]);
	iauCp(r[2], c[2]);
	
	return;
}	

void iauTr(double r[3][3], double rt[3][3])
{
	double wm[3][3];
	int i, j;
	
	
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			wm[i][j] = r[j][i];
		}
	}
	iauCr(wm, rt);
	
	return;
}

void iauRxp(double r[3][3], double p[3], double rp[3])
{
	double w, wrp[3];
	int i, j;
	
	
	/* Matrix r * vector p. */
	for (j = 0; j < 3; j++) {
		w = 0.0;
		for (i = 0; i < 3; i++) {
			w += r[j][i] * p[i];
		}
		wrp[j] = w;
	}
	
	/* Return the result. */
	iauCp(wrp, rp);
	
	return;
}	

void iauTrxp(double r[3][3], double p[3], double trp[3])
{
	double tr[3][3];
	
	
	/* Transpose of matrix r. */
	iauTr(r, tr);
	
	/* Matrix tr * vector p -> vector trp. */
	iauRxp(tr, p, trp);
	
	return;
}	

void iauPvtob(double elong, double phi, double hm,
              double xp, double yp, double sp, double theta,
              double pv[2][3])
/*
**  - - - - - - - - -
**   i a u P v t o b
**  - - - - - - - - -
**
**  Position and velocity of a terrestrial observing station.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     elong   double       longitude (radians, east +ve, Note 1)
**     phi     double       latitude (geodetic, radians, Note 1)
**     hm      double       height above ref. ellipsoid (geodetic, m)
**     xp,yp   double       coordinates of the pole (radians, Note 2)
**     sp      double       the TIO locator s' (radians, Note 2)
**     theta   double       Earth rotation angle (radians, Note 3)
**
**  Returned:
**     pv      double[2][3] position/velocity vector (m, m/s, CIRS)
**
**  Notes:
**
**  1) The terrestrial coordinates are with respect to the WGS84
**     reference ellipsoid.
**
**  2) xp and yp are the coordinates (in radians) of the Celestial
**     Intermediate Pole with respect to the International Terrestrial
**     Reference System (see IERS Conventions), measured along the
**     meridians 0 and 90 deg west respectively.  sp is the TIO locator
**     s', in radians, which positions the Terrestrial Intermediate
**     Origin on the equator.  For many applications, xp, yp and
**     (especially) sp can be set to zero.
**
**  3) If theta is Greenwich apparent sidereal time instead of Earth
**     rotation angle, the result is with respect to the true equator
**     and equinox of date, i.e. with the x-axis at the equinox rather
**     than the celestial intermediate origin.
**
**  4) The velocity units are meters per UT1 second, not per SI second.
**     This is unlikely to have any practical consequences in the modern
**     era.
**
**  5) No validation is performed on the arguments.  Error cases that
**     could lead to arithmetic exceptions are trapped by the iauGd2gc
**     function, and the result set to zeros.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
**     the Astronomical Almanac, 3rd ed., University Science Books
**     (2013), Section 7.4.3.3.
**
**  Called:
**     iauGd2gc     geodetic to geocentric transformation
**     iauPom00     polar motion matrix
**     iauTrxp      product of transpose of r-matrix and p-vector
**
**  This revision:   2013 October 9
**
**  SOFA release 2016-05-03
**
**  Copyright (C) 2016 IAU SOFA Board.  See notes at end.
*/
{
/* Earth rotation rate in radians per UT1 second */
   const double OM = 1.00273781191135448 * D2PI / DAYSEC;

   double xyzm[3], rpm[3][3], xyz[3], x, y, z, s, c;


/* Geodetic to geocentric transformation (WGS84). */
   (void) iauGd2gc(1, elong, phi, hm, xyzm);

/* Polar motion and TIO position. */
   iauPom00(xp, yp, sp, rpm);
   iauTrxp(rpm, xyzm, xyz);
   x = xyz[0];
   y = xyz[1];
   z = xyz[2];

/* Functions of ERA. */
   s = sin(theta);
   c = cos(theta);

/* Position. */
   pv[0][0] = c*x - s*y;
   pv[0][1] = s*x + c*y;
   pv[0][2] = z;

/* Velocity. */
   pv[1][0] = OM * ( -s*x - c*y );
   pv[1][1] = OM * (  c*x - s*y );
   pv[1][2] = 0.0;

/* Finished. */
}
