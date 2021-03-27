#include<cmath>

void CSZ(double& C, double& S, double z)
// universal functions C(z),S(z)
// reference:
//  R.R.Bate, D.D.Mueller and J.E.White
//   "Fundamentals of Astrodynamics" p196
{
	static double z0(5.e-1);
	if(fabs(z) < z0) {
		double dc, ds, C1(0);
		dc=C=0.5;
		ds=S=1./6;
		for(int i=4; C1!=C; i+=2) {
			dc *= -z/i/(i-1);
			ds *= -z/i/(i+1);
			C1 = C;
			C += dc;
			S += ds;
		}
	}
	else if(z>0) {
		double x = sqrt(z);
		C = (1 - cos(x))/z;
		S = (x - sin(x))/(x*z);
	}
	else {
        double x = sqrt(-z);
		C = (1 - cosh(x))/z;
		S = (x - sinh(x))/(x*z);
	}
}

void CSZ(double& C, double& S, double& dC, double& dS, double z)
// universal functions C(z),S(z),C'(z),S'(z)
// reference:
//  R.R.Bate, D.D.Mueller and J.E.White
//   "Fundamentals of Astrodynamics" p196
{
	static double z0(5.e-1);
	if(fabs(z) < z0) {
		double dc, ds, C1(0);
		dc = dC = -1./24;
		ds = dS = -1./120;
		C = 0.5;
		S = 1./6;
		for(int i=6, j=2; C1!=C; i+=2, j++) {
			dc *= z;
			ds *= z;
			C1 = C;
			C += dc;
			S += ds;
			dc/=-i; dc/=i-1;
			ds/=-i; ds/=i+1;
			dC += j*dc;
			dS += j*ds;
		}
	}
	else {
        if(z>0) {
			double x = sqrt(z);
			C = (1 - cos(x))/z;
			S = (x - sin(x))/(x*z);
		}
		else {
            double x = sqrt(-z);
			C = (1 - cosh(x))/z;
			S = (x - sinh(x))/(x*z);
		}
		dC = 0.5*(1 - z*S - 2*C)/z;
		dS = 0.5*(C - 3*S)/z;
	}
}

