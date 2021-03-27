#include "math.h"

#define DJM0 (2400000.5)/* Julian Date of Modified Julian Date zero */

int iauCal2jd(int iy, int im, int id, double *djm0, double *djm)
/*
**  - - - - - - - - - -
**   i a u C a l 2 j d
**  - - - - - - - - - -
**
**  Gregorian Calendar to Julian Date.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     iy,im,id  int     year, month, day in Gregorian calendar (Note 1)
**
**  Returned:
**     djm0      double  MJD zero-point: always 2400000.5
**     djm       double  Modified Julian Date for 0 hrs
**
**  Returned (function value):
**               int     status:
**                           0 = OK
**                          -1 = bad year   (Note 3: JD not computed)
**                          -2 = bad month  (JD not computed)
**                          -3 = bad day    (JD computed)
**
**  Notes:
**
**  1) The algorithm used is valid from -4800 March 1, but this
**     implementation rejects dates before -4799 January 1.
**
**  2) The Julian Date is returned in two pieces, in the usual SOFA
**     manner, which is designed to preserve time resolution.  The
**     Julian Date is available as a single number by adding djm0 and
**     djm.
**
**  3) In early eras the conversion is from the "Proleptic Gregorian
**     Calendar";  no account is taken of the date(s) of adoption of
**     the Gregorian Calendar, nor is the AD/BC numbering convention
**     observed.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 12.92 (p604).
**
**  This revision:  2013 August 7
**
**  SOFA release 2016-05-03
**
**  Copyright (C) 2016 IAU SOFA Board.  See notes at end.
*/
{
   int j, ly, my;
   long iypmy;

/* Earliest year allowed (4800BC) */
   const int IYMIN = -4799;

/* Month lengths in days */
   static const int mtab[]
                     = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};


/* Preset status. */
   j = 0;

/* Validate year and month. */
   if (iy < IYMIN) return -1;
   if (im < 1 || im > 12) return -2;

/* If February in a leap year, 1, otherwise 0. */
   ly = ((im == 2) && !(iy%4) && (iy%100 || !(iy%400)));

/* Validate day, taking into account leap years. */
   if ( (id < 1) || (id > (mtab[im-1] + ly))) j = -3;

/* Return result. */
   my = (im - 14) / 12;
   iypmy = (long) (iy + my);
   *djm0 = DJM0;
   *djm = (double)((1461L * (iypmy + 4800L)) / 4L
                 + (367L * (long) (im - 2 - 12 * my)) / 12L
                 - (3L * ((iypmy + 4900L) / 100L)) / 4L
                 + (long) id - 2432076L);

/* Return status. */
   return j;
}

int iauDat(int iy, int im, int id, double fd, double *deltat )
/*
 **  - - - - - - -
 **   i a u D a t
 **  - - - - - - -
 **
 **  For a given UTC date, calculate delta(AT) = TAI-UTC.
 **
 **     :------------------------------------------:
 **     :                                          :
 **     :                 IMPORTANT                :
 **     :                                          :
 **     :  A new version of this function must be  :
 **     :  produced whenever a new leap second is  :
 **     :  announced.  There are four items to     :
 **     :  change on each such occasion:           :
 **     :                                          :
 **     :  1) A new line must be added to the set  :
 **     :     of statements that initialize the    :
 **     :     array "changes".                     :
 **     :                                          :
 **     :  2) The constant IYV must be set to the  :
 **     :     current year.                        :
 **     :                                          :
 **     :  3) The "Latest leap second" comment     :
 **     :     below must be set to the new leap    :
 **     :     second date.                         :
 **     :                                          :
 **     :  4) The "This revision" comment, later,  :
 **     :     must be set to the current date.     :
 **     :                                          :
 **     :  Change (2) must also be carried out     :
 **     :  whenever the function is re-issued,     :
 **     :  even if no leap seconds have been       :
 **     :  added.                                  :
 **     :                                          :
 **     :  Latest leap second:  2015 June 30       :
 **     :                                          :
 **     :__________________________________________:
 **
 **  This function is part of the International Astronomical Union's
 **  SOFA (Standards Of Fundamental Astronomy) software collection.
 **
 **  Status:  support function.
 **
 **  Given:
 **     iy     int      UTC:  year (Notes 1 and 2)
 **     im     int            month (Note 2)
 **     id     int            day (Notes 2 and 3)
 **     fd     double         fraction of day (Note 4)
 **
 **  Returned:
 **     deltat double   TAI minus UTC, seconds
 **
 **  Returned (function value):
 **            int      status (Note 5):
 **                       1 = dubious year (Note 1)
 **                       0 = OK
 **                      -1 = bad year
 **                      -2 = bad month
 **                      -3 = bad day (Note 3)
 **                      -4 = bad fraction (Note 4)
 **                      -5 = internal error (Note 5)
 **
 **  Notes:
 **
 **  1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
 **     to call the function with an earlier date.  If this is attempted,
 **     zero is returned together with a warning status.
 **
 **     Because leap seconds cannot, in principle, be predicted in
 **     advance, a reliable check for dates beyond the valid range is
 **     impossible.  To guard against gross errors, a year five or more
 **     after the release year of the present function (see the constant
 **     IYV) is considered dubious.  In this case a warning status is
 **     returned but the result is computed in the normal way.
 **
 **     For both too-early and too-late years, the warning status is +1.
 **     This is distinct from the error status -1, which signifies a year
 **     so early that JD could not be computed.
 **
 **  2) If the specified date is for a day which ends with a leap second,
 **     the UTC-TAI value returned is for the period leading up to the
 **     leap second.  If the date is for a day which begins as a leap
 **     second ends, the UTC-TAI returned is for the period following the
 **     leap second.
 **
 **  3) The day number must be in the normal calendar range, for example
 **     1 through 30 for April.  The "almanac" convention of allowing
 **     such dates as January 0 and December 32 is not supported in this
 **     function, in order to avoid confusion near leap seconds.
 **
 **  4) The fraction of day is used only for dates before the
 **     introduction of leap seconds, the first of which occurred at the
 **     end of 1971.  It is tested for validity (0 to 1 is the valid
 **     range) even if not used;  if invalid, zero is used and status -4
 **     is returned.  For many applications, setting fd to zero is
 **     acceptable;  the resulting error is always less than 3 ms (and
 **     occurs only pre-1972).
 **
 **  5) The status value returned in the case where there are multiple
 **     errors refers to the first error detected.  For example, if the
 **     month and day are 13 and 32 respectively, status -2 (bad month)
 **     will be returned.  The "internal error" status refers to a
 **     case that is impossible but causes some compilers to issue a
 **     warning.
 **
 **  6) In cases where a valid result is not available, zero is returned.
 **
 **  References:
 **
 **  1) For dates from 1961 January 1 onwards, the expressions from the
 **     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
 **
 **  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
 **     the 1992 Explanatory Supplement.
 **
 **  Called:
 **     iauCal2jd    Gregorian calendar to JD
 **
 **  This revision:  2015 February 27
 **
 **  SOFA release 2016-05-03
 **
 **  Copyright (C) 2016 IAU SOFA Board.  See notes at end.
 */
{
	/* Release year for this version of iauDat */
	enum { IYV = 2015};
	
	/* Reference dates (MJD) and drift rates (s/day), pre leap seconds */
	static const double drift[][2] = {
		{ 37300.0, 0.0012960 },
		{ 37300.0, 0.0012960 },
		{ 37300.0, 0.0012960 },
		{ 37665.0, 0.0011232 },
		{ 37665.0, 0.0011232 },
		{ 38761.0, 0.0012960 },
		{ 38761.0, 0.0012960 },
		{ 38761.0, 0.0012960 },
		{ 38761.0, 0.0012960 },
		{ 38761.0, 0.0012960 },
		{ 38761.0, 0.0012960 },
		{ 38761.0, 0.0012960 },
		{ 39126.0, 0.0025920 },
		{ 39126.0, 0.0025920 }
	};
	
	/* Number of Delta(AT) expressions before leap seconds were introduced */
	enum { NERA1 = (int) (sizeof drift / sizeof (double) / 2) };
	
	/* Dates and Delta(AT)s */
	static const struct {
		int iyear, month;
		double delat;
	} changes[] = {
		{ 1960,  1,  1.4178180 },
		{ 1961,  1,  1.4228180 },
		{ 1961,  8,  1.3728180 },
		{ 1962,  1,  1.8458580 },
		{ 1963, 11,  1.9458580 },
		{ 1964,  1,  3.2401300 },
		{ 1964,  4,  3.3401300 },
		{ 1964,  9,  3.4401300 },
		{ 1965,  1,  3.5401300 },
		{ 1965,  3,  3.6401300 },
		{ 1965,  7,  3.7401300 },
		{ 1965,  9,  3.8401300 },
		{ 1966,  1,  4.3131700 },
		{ 1968,  2,  4.2131700 },
		{ 1972,  1, 10.0       },
		{ 1972,  7, 11.0       },
		{ 1973,  1, 12.0       },
		{ 1974,  1, 13.0       },
		{ 1975,  1, 14.0       },
		{ 1976,  1, 15.0       },
		{ 1977,  1, 16.0       },
		{ 1978,  1, 17.0       },
		{ 1979,  1, 18.0       },
		{ 1980,  1, 19.0       },
		{ 1981,  7, 20.0       },
		{ 1982,  7, 21.0       },
		{ 1983,  7, 22.0       },
		{ 1985,  7, 23.0       },
		{ 1988,  1, 24.0       },
		{ 1990,  1, 25.0       },
		{ 1991,  1, 26.0       },
		{ 1992,  7, 27.0       },
		{ 1993,  7, 28.0       },
		{ 1994,  7, 29.0       },
		{ 1996,  1, 30.0       },
		{ 1997,  7, 31.0       },
		{ 1999,  1, 32.0       },
		{ 2006,  1, 33.0       },
		{ 2009,  1, 34.0       },
		{ 2012,  7, 35.0       },
		{ 2015,  7, 36.0       }
	};
	
	/* Number of Delta(AT) changes */
	enum { NDAT = (int) (sizeof changes / sizeof changes[0]) };
	
	/* Miscellaneous local variables */
	int j, i, m;
	double da, djm0, djm;
	
	
	/* Initialize the result to zero. */
	*deltat = da = 0.0;
	
	/* If invalid fraction of a day, set error status and give up. */
	if (fd < 0.0 || fd > 1.0) return -4;
	
	/* Convert the date into an MJD. */
	j = iauCal2jd(iy, im, id, &djm0, &djm);
	
	/* If invalid year, month, or day, give up. */
	if (j < 0) return j;
	
	/* If pre-UTC year, set warning status and give up. */
	if (iy < changes[0].iyear) return 1;
	
	/* If suspiciously late year, set warning status but proceed. */
	if (iy > IYV + 5) j = 1;
	
	/* Combine year and month to form a date-ordered integer... */
	m = 12*iy + im;
	
	/* ...and use it to find the preceding table entry. */
	for (i = NDAT-1; i >=0; i--) {
		if (m >= (12 * changes[i].iyear + changes[i].month)) break;
	}
	
	/* Prevent underflow warnings. */
	if (i < 0) return -5;
	
	/* Get the Delta(AT). */
	da = changes[i].delat;
	
	/* If pre-1972, adjust for drift. */
	if (i < NERA1) da += (djm + fd - drift[i][0]) * drift[i][1];
	
	/* Return the Delta(AT) value. */
	*deltat = da;
	
	/* Return the status. */
	return j;
}

int iauJd2cal(double dj1, double dj2,
              int *iy, int *im, int *id, double *fd)
/*
 **  - - - - - - - - - -
 **   i a u J d 2 c a l
 **  - - - - - - - - - -
 **
 **  Julian Date to Gregorian year, month, day, and fraction of a day.
 **
 **  This function is part of the International Astronomical Union's
 **  SOFA (Standards Of Fundamental Astronomy) software collection.
 **
 **  Status:  support function.
 **
 **  Given:
 **     dj1,dj2   double   Julian Date (Notes 1, 2)
 **
 **  Returned (arguments):
 **     iy        int      year
 **     im        int      month
 **     id        int      day
 **     fd        double   fraction of day
 **
 **  Returned (function value):
 **               int      status:
 **                           0 = OK
 **                          -1 = unacceptable date (Note 3)
 **
 **  Notes:
 **
 **  1) The earliest valid date is -68569.5 (-4900 March 1).  The
 **     largest value accepted is 1e9.
 **
 **  2) The Julian Date is apportioned in any convenient way between
 **     the arguments dj1 and dj2.  For example, JD=2450123.7 could
 **     be expressed in any of these ways, among others:
 **
 **            dj1             dj2
 **
 **         2450123.7           0.0       (JD method)
 **         2451545.0       -1421.3       (J2000 method)
 **         2400000.5       50123.2       (MJD method)
 **         2450123.5           0.2       (date & time method)
 **
 **  3) In early eras the conversion is from the "proleptic Gregorian
 **     calendar";  no account is taken of the date(s) of adoption of
 **     the Gregorian calendar, nor is the AD/BC numbering convention
 **     observed.
 **
 **  Reference:
 **
 **     Explanatory Supplement to the Astronomical Almanac,
 **     P. Kenneth Seidelmann (ed), University Science Books (1992),
 **     Section 12.92 (p604).
 **
 **  This revision:  2013 August 7
 **
 **  SOFA release 2016-05-03
 **
 **  Copyright (C) 2016 IAU SOFA Board.  See notes at end.
 */
{
	/* Minimum and maximum allowed JD */
	const double DJMIN = -68569.5;
	const double DJMAX = 1e9;
	
	long jd, l, n, i, k;
	double dj, d1, d2, f1, f2, f, d;
	
	
	/* Verify date is acceptable. */
	dj = dj1 + dj2;
	if (dj < DJMIN || dj > DJMAX) return -1;
	
	/* Copy the date, big then small, and re-align to midnight. */
	if (dj1 >= dj2) {
		d1 = dj1;
		d2 = dj2;
	} else {
		d1 = dj2;
		d2 = dj1;
	}
	d2 -= 0.5;
	
	/* Separate day and fraction. */
	f1 = fmod(d1, 1.0);
	f2 = fmod(d2, 1.0);
	f = fmod(f1 + f2, 1.0);
	if (f < 0.0) f += 1.0;
	d = floor(d1 - f1) + floor(d2 - f2) + floor(f1 + f2 - f);
	jd = (long) floor(d) + 1L;
	
	/* Express day in Gregorian calendar. */
	l = jd + 68569L;
	n = (4L * l) / 146097L;
	l -= (146097L * n + 3L) / 4L;
	i = (4000L * (l + 1L)) / 1461001L;
	l -= (1461L * i) / 4L - 31L;
	k = (80L * l) / 2447L;
	*id = (int) (l - (2447L * k) / 80L);
	l = k / 11L;
	*im = (int) (k + 2L - 12L * l);
	*iy = (int) (100L * (n - 49L) + i + l);
	*fd = f;
	
	return 0;
}
