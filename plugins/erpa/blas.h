#ifndef ERPA_BLAS_H
#define ERPA_BLAS_H

#include"fortran.h"

namespace psi{ namespace erpa{

/**
 * name mangling dggev
 */
extern "C" {
    void dggev(char&JOBVL, char&JOBVR, long int &N, double*A,long int&LDA,double*B,long int &LDB,
        double * ALPHAR, double * ALPHAI, double * BETA, double * VL, long int &LDVL, double * VR, long int &LDVR,
        double*WORK,long int&LWORK,long int&INFO);
};
inline void DGGEV(char&JOBVL, char&JOBVR, long int &N, double*A,long int&LDA,double*B,long int &LDB,
        double * ALPHAR, double * ALPHAI, double * BETA, double * VL, long int &LDVL, double * VR, long int &LDVR,
        double*WORK,long int&LWORK,long int&INFO){
    dggev(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO);
}
/**
 * name manging dgeev
 */
extern "C" {
void F77NAME(dgeev)(char &jobvl,char &jobvr,integer &n,doublereal *a,integer &lda,
                    doublereal *wr,doublereal *wi, doublereal *vl,integer &ldvl,doublereal *vr,
                    integer &ldvr,doublereal * work,integer &lwork,integer &info);
};
inline void DGEEV(char &jobvl,char &jobvr,integer &n,doublereal*a,integer &lda,
                  doublereal*wr,doublereal*wi, doublereal*vl,integer&ldvl,doublereal*vr,
                  integer &ldvr,doublereal * work,integer &lwork,integer&info)
{
  F77NAME(dgeev)(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info);
};



/**
  * diagonalize general matrix and keep eigenvectors
  */
void NonSymmetricEigenvalueEigenvector(long int dim, double * M, double * eigval, double * el, double * er) {
  integer info;
  char vl = 'V';
  char vr = 'V';
  integer n = dim;
  integer lwork = 4*n;
  double * work  = (doublereal*)malloc(lwork*sizeof(doublereal));
  double * wi = (double*)malloc(n*sizeof(double));

  DGEEV(vl,vr,n,M,n,eigval,wi,el,n,er,n,work,lwork,info);

  // sort eigenvalues and eigenvectors
  //int count = 0;
  //int*skip=(int*)malloc(n*sizeof(int));
  //for (int i = 0; i < n; i++) skip[i] = 0;
  //for (int i = 0; i < n; i++) {
  //    double min = 1e99;
  //    int mini   = -1;
  //    for (int j = 0; j < n; j++) {
  //        if (skip[j]) continue;
  //        if (eigval[j] < min) {
  //            min  = eigval[j];
  //            mini = j;
  //        }
  //    }
  //    skip[mini] = 1;
  //    wi[i] = min;
  //}
  //C_DCOPY(n,wi,1,eigval,1);

  //free(skip);

  free(wi);
  free(work);
}
/**
  * diagonalize general matrix, don't compute eigenvectors
  */
void NonSymmetricEigenvalue(long int dim, double * M, double * eigval) {
  integer info;
  char vl = 'N';
  char vr = 'N';
  integer n = dim;
  integer lwork = 4*n;
  double * work  = (doublereal*)malloc(lwork*sizeof(doublereal));
  double * wi = (double*)malloc(n*sizeof(double));

  double * eigvec, * leigvec;

  DGEEV(vl,vr,n,M,n,eigval,wi,leigvec,n,eigvec,n,work,lwork,info);

  // sort eigenvalues and eigenvectors
  int count = 0;
  int*skip=(int*)malloc(n*sizeof(int));
  for (int i = 0; i < n; i++) skip[i] = 0;
  for (int i = 0; i < n; i++) {
      double min = 1e99;
      int mini   = -1;
      for (int j = 0; j < n; j++) {
          if (skip[j]) continue;
          if (eigval[j] < min) {
              min  = eigval[j];
              mini = j;
          }
      }
      skip[mini] = 1;
      wi[i] = min;
  }
  C_DCOPY(n,wi,1,eigval,1);
  //free(eigvecl);
  free(skip);
  free(wi);
  free(work);
}
void GeneralizedEigenvalueProblem(long int N, double * A, double * B, double * c, double * eig){
    char JOBVR = 'N';
    char JOBVL = 'N';

    long int LDA  = N;
    long int LDB  = N;
    long int LDVR = N;
    long int LDVL = N;
    long int LWORK = 8*N;

    double*WORK=(double*)malloc(LWORK*sizeof(double));
    double*ALPHAR=(double*)malloc(N*sizeof(double));
    double*ALPHAI=(double*)malloc(N*sizeof(double));
    double*BETA=(double*)malloc(N*sizeof(double));
    double*VR=(double*)malloc(N*N*sizeof(double));
    double*VL=(double*)malloc(N*N*sizeof(double));

    integer INFO=0;
    DGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO);

    //printf("hey info %5i, dimension %5i\n",INFO,N);

    outfile->Printf("\n");
    outfile->Printf("    ==> positive excitation energies <==\n");
    outfile->Printf("\n");
    for (long int i = 0; i < N; i++) {
        eig[i] = ALPHAR[i] / BETA[i];
        if ( eig[i] > 0.0 && fabs(ALPHAI[i]/BETA[i]) < 1e-6) {
            outfile->Printf("%20.12lf %20.12lf %20.12lf\n",eig[i], eig[i] * 27.21138, ALPHAI[i]/BETA[i]);
            //for (long int j = 0; j < N; j++) { 
            //    printf("        %20.12lf %20.12lf\n",VR[i*N+j],VL[i*N+j]);
            //    //printf("        %20.12lf %20.12lf\n",VR[j*N+i],VL[j*N+i]);
            //}
        }else {
            //eig[i] = 0.0;
        }
    }
    outfile->Printf("\n");
    C_DCOPY(N*N,VR,1,B,1);
    C_DCOPY(N*N,VL,1,A,1);

    free(VR);
    free(VL);
    free(WORK);
    free(ALPHAR);
    free(ALPHAI);
    free(BETA);
}
/**
 * name mangling dggev
 */
extern "C" {
    void dsygv(long int &ITYPE,char&JOBZ,char&UPLO, long int &N, double*A,long int&LDA,double*B,long int &LDB,
        double *W, double*WORK,long int&LWORK,long int&INFO);
};
inline void DSYGV(long int &ITYPE,char&JOBZ,char&UPLO, long int &N, double*A,long int&LDA,double*B,long int &LDB,
        double *W, double*WORK,long int&LWORK,long int&INFO){
    dsygv(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO);
}

int SymmetricGeneralizedEigenvalueProblem(long int N, double * A, double * B, double * c, double * eig){

    long int ITYPE = 1; // Ax = l Bx (actually solve Bx = 1/l Ax)
    char JOBZ = 'V';
    char UPLO = 'U';
    long int LDA  = N;
    long int LDB  = N;

    long int LWORK = 3*N-1;
    double*WORK=(double*)malloc(LWORK*sizeof(double));

    integer INFO=0;
    //printf("hey itype %5i\n",ITYPE);
    DSYGV(ITYPE,JOBZ,UPLO,N,B,LDA,A,LDB,eig,WORK,LWORK,INFO);

    return (int)INFO;

    free(WORK);
}

// diagonalize real, nonsymmetric matrix
void NonsymmetricEigenvalue(long int N, double * A, double * VL, double * VR, double * WR, double *WI){

    char JOBVL = 'V';
    char JOBVR = 'V';
    long int LDA  = N;
    long int LDVL = N;
    long int LDVR = N;
    long int LWORK = 4*N;
    double * WORK = (double*)malloc(LWORK*sizeof(double));
    long int INFO;

    DGEEV(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO);

    // kill complex eigenvalues
    //for (int i = 0; i < N; i++) {
    //    if ( fabs(WI[i]) > 1e-6 ) {
    //        WR[i] = 0.0;
    //        WI[i] = 0.0;
    //    }
    //}

    free(WORK);
}

}} // end of namespace

#endif
