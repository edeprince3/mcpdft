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

#include"blas.h"
#include<stdlib.h>
#include<stdio.h>

namespace psi { namespace fnocc{

// position in a symmetric packed matrix
long int Position(long int i,long int j){
  if (i<j){
    return ((j*(j+1))>>1)+i;
  }
  return ((i*(i+1))>>1)+j;
}

/**
 *  Diagonalize a hermitian matrix
 */
void DiagonalizeHermitianMatrix(long int N,double*re,double*im,double*W){
    WRAP(re,im,W,N);

    /*char JOBZ = 'V';
    char UPLO = 'U';
    integer LDA = N;
    integer LWORK = 2*N-1;
    complexdoublereal*WORK=(complexdoublereal*)malloc(1000*LWORK*sizeof(complexdoublereal)); 
    complexdoublereal*A=(complexdoublereal*)malloc(1000*N*sizeof(complexdoublereal)); 
    for (int i = 0; i < N; i++) {
        A[2*i]   = re[i];
        A[2*i+1] = im[i];
    }

    integer INFO=0;
    doublereal*RWORK=(doublereal*)malloc(1000*(3*N-2)*sizeof(doublereal)); 
//for (int i = 0; i < N*N; i++) printf("%5i %20.12lf %20.12lf\n",i,A[i].val[0],A[i].val[1]);
    ZHEEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO);
//for (int i = 0; i < N*N; i++) printf("%5i %20.12lf %20.12lf\n",i,A[i].val[0],A[i].val[1]);
    //free(WORK);
    //free(RWORK);*/
}

}}
