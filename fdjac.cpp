// W.H.Press et al "Numerical Recipes in C++"
//  http://numerical.recipes

#include<cmath>
#include "Vec.h"

void fdjac(const Vec_DP &x, const Vec_DP &y, Mat_DP &df,
	      void vecfunc(const Vec_DP &, Vec_DP &))
{
	const double EPS=1.0e-8;
	int i,j;
	double h;
	
	int n=x.size();
	int m=y.size();
	Vec_DP f(m),xh(x);
	for (j=0;j<n;j++) {
		h=EPS*fabs(x[j]);
		if (h == 0.0) h=EPS;
		xh[j]+=h;
		h=xh[j]-x[j];
		vecfunc(xh,f);
		xh[j]=x[j];
		for (i=0;i<m;i++)
			df[i][j]=(f[i]-y[i])/h;
	}
}
