// W.H.Press et al "Numerical Recipes in C++"
//  http://numerical.recipes

#include "Vec.h"
#include "util.h"
#include<cmath>

void fdjac(const Vec_DP&, const Vec_DP&, Mat_DP&,
		   void(const Vec_DP &, Vec_DP &));
void gaussj(Mat_DP&, Vec_DP&);

void covsrt(Mat_DP &covar, const Vec<bool> &ia, const int mfit)
{
	int i,j,k;
	
	int ma=ia.size();
	for (i=mfit;i<ma;i++)
		for (j=0;j<i+1;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit-1;
	for (j=ma-1;j>=0;j--) {
		if (ia[j]) {
			for (i=0;i<ma;i++) SWAP(covar[i][k],covar[i][j]);
			for (i=0;i<ma;i++) SWAP(covar[k][i],covar[j][i]);
			k--;
		}
	}
}

void mrqcof(const Vec_DP &sig, const Vec_DP &a,
		const Vec<bool> &ia, Mat_DP &alpha, Vec_DP &beta, double &chisq,
		void funcs(const Vec_DP &, Vec_DP &))
{
	int i,j,k,l,m,mfit=0;
	double wt,sig2i,dy;
	
	int ndata=sig.size();
	int ma=a.size();
	Vec_DP ymod(ndata);
	Mat_DP dyda(ndata,ma);
	for (j=0;j<ma;j++)
		if (ia[j]) mfit++;
	for (j=0;j<mfit;j++) {
		for (k=0;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	chisq=0.0;
	funcs(a,ymod);
	fdjac(a,ymod,dyda,funcs);
	for (i=0;i<ndata;i++) {
		sig2i=1.0/(sig[i]*sig[i]);
		dy=-ymod[i];
		for (j=0,l=0;l<ma;l++) {
			if (ia[l]) {
				wt=dyda[i][l]*sig2i;
				for (k=0,m=0;m<l+1;m++)
					if (ia[m]) alpha[j][k++] += wt*dyda[i][m];
				beta[j++] += dy*wt;
			}
		}
		chisq += dy*dy*sig2i;
	}
	for (j=1;j<mfit;j++)
		for (k=0;k<j;k++) alpha[k][j]=alpha[j][k];
}

void mrqmin(const Vec_DP &sig, Vec_DP &a,
	const Vec<bool> &ia, Mat_DP &covar, Mat_DP &alpha, double &chisq,
	void funcs(const Vec_DP &, Vec_DP &), double &alamda)
{
	static int mfit;
	static double ochisq;
	int j,k,l;

	int ma=a.size();
	static Vec_DP atry(ma),beta(ma),da(ma);
	if (alamda < 0.0) {
		mfit=0;
		for (j=0;j<ma;j++)
			if (ia[j]) mfit++;
		alamda=0.001;
		mrqcof(sig,a,ia,alpha,beta,chisq,funcs);
		ochisq=chisq;
		for (j=0;j<ma;j++) atry[j]=a[j];
	}
	Mat_DP temp(mfit,mfit);
	for (j=0;j<mfit;j++) {
		for (k=0;k<mfit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]*=(1.0+alamda);
		for (k=0;k<mfit;k++) temp[j][k]=covar[j][k];
		da[j]=beta[j];
	}
	gaussj(temp,da);
	for (j=0;j<mfit;j++)
		for (k=0;k<mfit;k++) covar[j][k]=temp[j][k];
	if (alamda == 0.0) {
		covsrt(covar,ia,mfit);
		covsrt(alpha,ia,mfit);
		return;
	}
	for (j=0,l=0;l<ma;l++)
		if (ia[l]) atry[l]=a[l]+da[j++];
	mrqcof(sig,atry,ia,covar,da,chisq,funcs);
	if (chisq < ochisq) {
		alamda *= 0.1;
		ochisq=chisq;
		for (j=0;j<mfit;j++) {
			for (k=0;k<mfit;k++) alpha[j][k]=covar[j][k];
				beta[j]=da[j];
		}
		for (l=0;l<ma;l++) a[l]=atry[l];
	} else {
		alamda *= 10.0;
		chisq=ochisq;
	}
}

void mrqmin(const Vec_DP &sig, Vec_DP &a, const Vec<bool>& ia,
	    Mat_DP &covar, double &chisq,
	    void funcs(const Vec_DP&, Vec_DP&))
{
	static double eps(1.e-3);
	static int IMAX(100), JMAX(10);
	double	alamda(-1), ochisq;
	int	i,j;
	Mat_DP alpha(a.size(), a.size());
	mrqmin(sig, a, ia, covar, alpha, chisq, funcs, alamda);
	for(i=j=0; i<IMAX && chisq>0; i++) {	
		ochisq = chisq;
		mrqmin(sig, a, ia, covar, alpha, chisq, funcs, alamda);
		if(ochisq-chisq <= eps*chisq) j++;
		else j=0;
		if(j==JMAX) break;
	}
	alamda=0.0;
	mrqmin(sig, a, ia, covar, alpha, chisq, funcs, alamda);
}
