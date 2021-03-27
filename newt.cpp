// W.H.Press et al "Numerical Recipes in C++"
//  http://numerical.recipes

#include <cmath>
#include <limits>
#include "Vec.h"
#include "util.h"
using namespace std;

Vec_DP *fvec_p;
void (*nrfuncv)(const Vec_DP&, Vec_DP&);

double fmin(const Vec_DP &x)
{
	int i;
	double sum;
	
	Vec_DP &fvec=*fvec_p;
	nrfuncv(x,fvec);
	int n=x.size();
	for (sum=0.0,i=0;i<n;i++) sum += SQR(fvec[i]);
	return 0.5*sum;
}

void lnsrch(const Vec_DP &xold, const double fold, const Vec_DP &g, Vec_DP &p,
		Vec_DP &x, double &f, const double stpmax, bool &check,
	       double func(const Vec_DP &))
{
	const double ALF=1.0e-4, TOLX=numeric_limits<double>::epsilon();
	const int MAXITS=200;
	int i,its;
	double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
	double rhs1,rhs2,slope,sum,temp,test,tmplam;
	
	int n=xold.size();
	check=false;
	sum=0.0;
	for (i=0;i<n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=0;i<n;i++) p[i] *= stpmax/sum;
	slope=0.0;
	for (i=0;i<n;i++)
		slope += g[i]*p[i];
	if (slope >= 0.0) error("Roundoff problem in lnsrch.");
	test=0.0;
	for (i=0;i<n;i++) {
		temp=fabs(p[i])/MAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (its=0;its<MAXITS;its++) {
		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
		f=func(x);
		if (alam < alamin) {
			for (i=0;i<n;i++) x[i]=xold[i];
			check=true;
			return;
		} else if (f <= fold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(f-fold-slope));
			else {
				rhs1=f-fold-alam*slope;
				rhs2=f2-fold-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc < 0.0) tmplam=0.5*alam;
					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
					else tmplam=-slope/(b+sqrt(disc));
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = f;
		alam=MAX(tmplam,0.1*alam);
	}
}

void newt(Vec_DP &x, bool &check, void vecfunc(const Vec_DP &, Vec_DP &))
{
	void ludcmp(Mat_DP&, Vec<int>&, double&);
	void lubksb(const Mat_DP&, const Vec<int>&, Vec_DP&);
	void fdjac(const Vec_DP&, const Vec_DP&, Mat_DP&,
		   void(const Vec_DP &, Vec_DP &));
	const int MAXITS=200;
	const double TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=100.0;
	const double TOLX=numeric_limits<double>::epsilon();
	int i,j,its;
	double d,den,f,fold,stpmax,sum,temp,test;
	
	int n=x.size();
	Vec<int> indx(n);
	Vec_DP g(n),p(n),xold(n);
	Mat_DP fjac(n,n);
	fvec_p=new Vec_DP(n);
	nrfuncv=vecfunc;
	Vec_DP &fvec=*fvec_p;
	f=fmin(x);
	test=0.0;
	for (i=0;i<n;i++)
		if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
	if (test < 0.01*TOLF) {
		check=false;
		delete fvec_p;
		return;
	}
	sum=0.0;
	for (i=0;i<n;i++) sum += SQR(x[i]);
	stpmax=STPMX*MAX(sqrt(sum),double(n));
	for (its=0;its<MAXITS;its++) {
		fdjac(x,fvec,fjac,vecfunc);
		for (i=0;i<n;i++) {
			sum=0.0;
			for (j=0;j<n;j++) sum += fjac[j][i]*fvec[j];
			g[i]=sum;
		}
		for (i=0;i<n;i++) xold[i]=x[i];
		fold=f;
		for (i=0;i<n;i++) p[i] = -fvec[i];
		ludcmp(fjac,indx,d);
		lubksb(fjac,indx,p);
		lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin);
		test=0.0;
		for (i=0;i<n;i++)
			if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
		if (test < TOLF) {
			check=false;
			delete fvec_p;
			return;
		}
		if (check) {
			test=0.0;
			den=MAX(f,0.5*n);
			for (i=0;i<n;i++) {
				temp=fabs(g[i])*MAX(fabs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			check=(test < TOLMIN);
			delete fvec_p;
			return;
		}
		test=0.0;
		for (i=0;i<n;i++) {
			temp=(fabs(x[i]-xold[i]))/MAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) {
			delete fvec_p;
			return;
		}
	}
//	error("MAXITS exceeded in newt");
	check=true;
}

bool newt(Vec_DP &x, void vecfunc(const Vec_DP &, Vec_DP &))
{
    bool check;
    newt(x,check,vecfunc);
    return check;
}