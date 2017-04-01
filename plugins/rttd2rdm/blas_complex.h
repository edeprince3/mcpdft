/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef BLAS_COMPLEX_H
#define BLAS_COMPLEX_H

/**
 * fortran-ordered blas routines
 */

#ifndef FC_SYMBOL
#define FC_SYMBOL 2
#endif

#if   FC_SYMBOL==1
#define F77NAME(x) x
#elif FC_SYMBOL==2
#define F77NAME(x) x##_
#endif

#include<stdlib.h>
//#include<complex.h>

typedef long int integer;
typedef double doublereal;
struct mycomplex {
    double val[2];
    //long double val;
};
//typedef double _Complex complexdoublereal;
//typedef struct mycomplex complexdoublereal;
typedef double complexdoublereal;
//typedef void complexdoublereal;

namespace psi{ namespace fnocc{

/**
 * diagonalize a hermitian matrix
 */
void DiagonalizeHermitianMatrix(integer N,doublereal*re,doublereal*im,doublereal*W);

/**
 * name mangling fortran complex number wrapper
 */
extern "C" {
    void F77NAME(wrap)(doublereal*re,doublereal*im,doublereal*w,integer&n);
};
inline void WRAP(doublereal*re,doublereal*im,doublereal*w,integer&n){
    F77NAME(wrap)(re,im,w,n);
}

/**
 * name mangling cheev
 */
extern "C" {
    void F77NAME(zheev)(char&JOBZ,char&UPLO,integer&N,complexdoublereal*A,integer&LDA,doublereal*W,complexdoublereal*WORK,integer&LWORK,doublereal*RWORK,integer&INFO);
};
inline void ZHEEV(char&JOBZ,char&UPLO,integer&N,complexdoublereal*A,integer&LDA,doublereal*W,complexdoublereal*WORK,integer&LWORK,doublereal*RWORK,integer&INFO){
    F77NAME(zheev)(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO);
}

}}

#endif
