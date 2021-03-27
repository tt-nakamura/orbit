#define PI 3.14159265358979
#define PI2 6.28318530717959
#define GRAV 0.01720209895// gaussian grav const
#define C_LIGHT 173.14463331135// velocity of light / AU/day
#define DEGREE 0.0174532925199433// radian / degree
#define GM 0.000295912208285591// mu = GRAV^2
#define TT_TAI 32.184// (dynamical time) - (atomic time) / sec

#define DAU (149597870e3)// Astronomical unit (m)
#define DJM0 (2400000.5)// Julian Date of Modified Julian Date zero
#define DAYSEC (86400.0)// Seconds per day

#define	hms2h(h,m,s) ((((s)/60.+(m))/60.)+(h))
#define OBLIQUITY 0.409092804222329// angle between ecliptic and equator
#define cos_o 0.917482062069182// cos(OBLIQUITY)
#define sin_o 0.397777155931914// sin(OBLIQUITY)