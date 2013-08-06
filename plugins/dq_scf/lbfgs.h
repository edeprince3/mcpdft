#ifndef LBFGS_H
#define LBFGS_H
#include "ftoc.h"

extern "C" {
  void F77NAME(lbfgs)(integers &n,
		      integers &m,
		      double *x,
		      double &f,
		      double *g,
		      logicals &diagco,
		      double *diag,
		      integers *iprint,
		      double &eps,
		      double &xtol,
		      double &ftol,
		      double &stpmax,
		      int &maxfev,
		      double *w,
		      integers &iflag,
		      integers &info);
};

inline void LBFGS(int &n, int &m, double *x, double &f, double *g, int &diagco,
		  double *diag, int *iprint, double &eps, double &xtol,
		  double &ftol, double &stpmax, int &maxfev, double *w, 
		  int &iflag, int &info)
{
  F77NAME(lbfgs)(n,m,x,f,g,diagco,diag,iprint,eps,xtol,ftol,stpmax,maxfev,w,
		 iflag,info);
};

#endif
