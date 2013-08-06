#include"sdp.h"
#include<libpsio/psio.hpp>
#include<libmints/mints.h>
#include<../bin/fnocc/blas.h>
#include"lib/ftoc.h"
#include"lib/lbfgs.h"
#include<time.h>

using namespace psi;
using namespace fnocc;
namespace psi{ namespace dq_scf{

SDPSolver::SDPSolver(boost::shared_ptr<Wavefunction> reference_wavefunction,Options & options):
  Wavefunction(options,_default_psio_lib_)
{
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}


SDPSolver::~SDPSolver()
{
    free(V2);
    free(V1);
    free(D1);
    free(Q1);
}

double SDPSolver::compute_energy(){
    fflush(outfile);
    fprintf(outfile,"\n\n");
    fprintf(outfile, "        *******************************************************\n");
    fprintf(outfile, "        *                                                     *\n");
    fprintf(outfile, "        *                       DQ-SCF                        *\n");
    fprintf(outfile, "        *             Variational 1-RDM Optimization          *\n");
    fprintf(outfile, "        *                                                     *\n");
    fprintf(outfile, "        *                   Eugene DePrince                   *\n");
    fprintf(outfile, "        *                                                     *\n");
    fprintf(outfile, "        *******************************************************\n");
    fprintf(outfile,"\n\n");
    fflush(outfile);

    // get 1- and 2-electron integrals (TODO: these are held in memory currently)
    K2();

    double enuc = Process::environment.molecule()->nuclear_repulsion_energy();

    // lbfgs:
    int iprint[2];
    iprint[0]       = -1;
    iprint[1]       = 0;
    int mbfgs       = 5;
    int n           = 2*nmo*nmo;
    double xtol     = 1.0e-7;
    double ftol     = 1.0e-7;
    int maxfev      = 100;
    double stpmax   = 1.0e2;
    int diagco      = 0;
    double eps      = 1.0e-7;
    double * hess_f = (double*)malloc(n*sizeof(double));
    double * w      = (double*)malloc((n*(2*mbfgs+1)+2*mbfgs)*sizeof(double));
    int iflag=0;
    int info=0;

    double *  vars = (double*)malloc( n * sizeof(double));
    double * dvars = (double*)malloc( n * sizeof(double));

    R1d   = vars;
    R1q   = vars + nmo*nmo;
    dR1d  = dvars;
    dR1q  = dvars + nmo*nmo;

    memset((void*)vars,'\0',n*sizeof(double));
    memset((void*)dvars,'\0',n*sizeof(double));

    // random guess
    Guess();
    double dum = K2D1();

    int oiter = 0;
    // outer iterations:
    double max = 1000.0;
    double enow = 0.0;
    double elast = 1.0;
    fprintf(outfile,"\n");
    fprintf(outfile,"    reference energy: %20.12lf\n",dum+enuc);
    fprintf(outfile,"\n");
    fprintf(outfile,"    oiter");
    fprintf(outfile,"  iiter");
    fprintf(outfile,"           f");
    fprintf(outfile,"           e");
    fprintf(outfile,"      Tr(D1)");
    fprintf(outfile,"   max error");
    fprintf(outfile,"  error norm\n");

    double conv = 1e9;
    do {
        elast = enow;
        info = 0;
        for (int i = 0; i < n; i++) dvars[i] = 0.0;

        double f = eval();
        int iiter = 0;

        // inner iterations:
        do {
            LBFGS(n,mbfgs,vars,f,dvars,diagco,hess_f,iprint,eps,xtol,ftol,stpmax,maxfev,w,iflag,info);
            if (iflag==1){
                memset((void*)dvars,'\0',n*sizeof(double));
                f = eval();
            }
            else break;

            if (iflag == 1){
                iiter++;
            }
        }while(iiter<100000);
        iflag=info=0;
        
        double newmax = 0.0;
        int imax = -1;
        for (int i = 0; i < nconstraints; i++) {
            if ( fabs(cerror[i]) > newmax ) {
                newmax = fabs(cerror[i]);
                imax = i;
             }
        }
        if ( newmax < 0.25 * max ){
            f = eval();
            for (int i = 0; i < nconstraints; i++) {
                lambda[i] -= cerror[i] / mu;
            }
        }else{
            mu /= 10.0;
        }
        max = newmax;
        oiter++;
        enow = energy;
        conv = F_DNRM2(nconstraints,cerror,1);
        fprintf(outfile,"      %3i %6i %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf\n",oiter,iiter,f+enuc,enow+enuc,cval[0],max,conv);fflush(outfile);

    }while(conv > 1e-7 );

    memset((void*)dvars,'\0',n*sizeof(double));
    eval();

    fprintf(outfile,"\n");
    fprintf(outfile,"    optimized energy: %20.12lf\n",energy+enuc);
    fprintf(outfile,"\n");
    Process::environment.globals["CURRENT ENERGY"] = energy+enuc;

    free(vars);
    free(dvars);
    free(hess_f);
    free(w);
}

void SDPSolver::BuildD1(){
    F_DGEMM('n','t',nmo,nmo,nmo,1.0,R1d,nmo,R1d,nmo,0.0,D1,nmo);
}

void SDPSolver::BuildQ1(){

    int offset = 1;

    F_DGEMM('n','t',nmo,nmo,nmo,1.0,R1q,nmo,R1q,nmo,0.0,Q1,nmo);
    for (int i = 0; i < nmo; i++) {
        for (int j = 0; j < nmo; j++) {
            cval[offset+i*nmo+j] = Q1[i*nmo+j] + D1[i*nmo+j] - (i==j);
            cerror[offset+i*nmo+j] = cval[offset+i*nmo+j] - c[offset+i*nmo+j];
        }
    }

    // derivative of constraints: d1ij + q1ij - dij = 0
    F_DGEMM('n','n',nmo,nmo,nmo,-2.0,lambda+offset,nmo,R1d,nmo,1.0,dR1d,nmo);
    F_DGEMM('n','n',nmo,nmo,nmo,4.0/mu,cerror+offset,nmo,R1d,nmo,1.0,dR1d,nmo);
    F_DGEMM('n','n',nmo,nmo,nmo,-2.0,lambda+offset,nmo,R1q,nmo,0.0,dR1q,nmo);
    F_DGEMM('n','n',nmo,nmo,nmo,4.0/mu,cerror+offset,nmo,R1q,nmo,1.0,dR1q,nmo);

    /*for (int k = 0; k < nmo; k++) {
        for (int beta = 0; beta < nmo; beta++) {
            dR1q[beta*nmo+k] = 0.0;
            for (int i = 0; i < nmo; i++) {
                dR1d[beta*nmo+k] -= 2.0 * R1d[beta*nmo+i] * lambda[offset+i*nmo+k];
                dR1d[beta*nmo+k] += 2.0 * 2.0 * R1d[beta*nmo+i] * 1.0/mu * cerror[offset+i*nmo+k];
                dR1q[beta*nmo+k] -= 2.0 * R1q[beta*nmo+i] * lambda[offset+i*nmo+k];
                dR1q[beta*nmo+k] += 2.0 * 2.0 * R1q[beta*nmo+i] * 1.0/mu * cerror[offset+i*nmo+k];
                //for (int j = 0; j < nmo; j++) {
                //    dR1d[beta*nmo+k] -= 2.0 * R1d[beta*nmo+i] * lambda[offset+i*nmo+j] * (j==k);
                //    dR1d[beta*nmo+k] += 2.0 * 2.0 * R1d[beta*nmo+i] * 1.0/mu * cerror[offset+i*nmo+j] * (j==k);
                //    //dR1d[beta*nmo+k] -= R1d[beta*nmo+j] * lambda[offset+i*nmo+j] * (i==k);
                //    //dR1d[beta*nmo+k] += 2.0 * R1d[beta*nmo+j] * 1.0/mu * cerror[offset+i*nmo+j] * (i==k);

                //    dR1q[beta*nmo+k] -= 2.0 * R1q[beta*nmo+i] * lambda[offset+i*nmo+j] * (j==k);
                //    dR1q[beta*nmo+k] += 2.0 * 2.0 * R1q[beta*nmo+i] * 1.0/mu * cerror[offset+i*nmo+j] * (j==k);
                //    //dR1q[beta*nmo+k] -= R1q[beta*nmo+j] * lambda[offset+i*nmo+j] * (i==k);
                //    //dR1q[beta*nmo+k] += 2.0 * R1q[beta*nmo+j] * 1.0/mu * cerror[offset+i*nmo+j] * (i==k);
                //}
            }
        }
    }*/
}

double SDPSolver::eval(){
    BuildD1();
    energy = K2D1();
    nrm = TraceD1();
    BuildQ1();

    for (int i = 0; i < nconstraints; i++) {
        cerror[i] = cval[i] - c[i];
    }
    double f = energy;
    for (int i = 0; i < nconstraints; i++) {
         f += - lambda[i] * cerror[i] + 1.0 / mu * cerror[i]*cerror[i];
    }
    return f;
}

double SDPSolver::TraceD1(){
    double dum = 0.0;
    for (int i = 0; i < nmo; i++) {
        dum += D1[i*nmo+i];
    }
    cval[0] = dum;
    nrm = dum;

    cerror[0] = cval[0] - c[0];

    // derivative of contraints:
    //for (int p = 0; p < nmo; p++) {
    //    for (int q = 0; q < nmo; q++) {
    //        // trace of d1
    //        dR1d[q*nmo+p] += 2.0 * (2.0/mu*cerror[0] - lambda[0] ) * R1d[q*nmo+p];
    //    }
    //}
    F_DAXPY(nmo*nmo,2.0 * (2.0/mu*cerror[0] - lambda[0] ),R1d,1,dR1d,1);
    return nrm;
}

void SDPSolver::Guess(){
    srand(time(NULL));
    memset((void*)R1q,'\0',nmo*nmo*sizeof(double));
    memset((void*)R1d,'\0',nmo*nmo*sizeof(double));
    for (int i = nalpha_; i < nmo; i++){
        R1q[i*nmo+i] = 1.0;
    }
    for (int i = 0; i < nalpha_; i++){
        R1d[i*nmo+i] = 1.0;
    }
    for (int i = 0; i < nmo*nmo; i++) {
        R1d[i] += ((double)rand()/RAND_MAX - 0.5 ) * 2.0 / 1000.0;
    }

    BuildD1();
    nrm    = TraceD1();
    nrm = c[0] / cval[0];
    for (int i = 0; i < nmo*nmo; i++) {
        R1d[i] *= nrm;
    }
    for (int i = 0; i < nmo; i++) {
        for (int j = 0; j < nmo; j++) {
            R1q[i*nmo+j] = (i==j) - R1d[i*nmo+j];
        }
    }
}

double SDPSolver::K2D1(){

    double en = F_DDOT(nmo*nmo,D1,1,V1,1);
    for (int i = 0; i < nmo; i++) {
        for (int k = 0; k < nmo; k++) {
            double dum = 0.0;
            for (int j = 0; j < nmo; j++) {
                for (int l = 0; l < nmo; l++) {
                    dum += D1[j*nmo+l] * ( V2[i*nmo*nmo*nmo+k*nmo*nmo+j*nmo+l] - 0.5 * V2[i*nmo*nmo*nmo+l*nmo*nmo+k*nmo+j]);
                }
            }
            en += dum * D1[i*nmo+k];
        }
    }
    //for (int i = 0; i < nmo; i++) {
    //    for (int j = 0; j < nmo; j++) {
    //        en += D1[i*nmo+j] * V1[i*nmo+j];
    //        for (int k = 0; k < nmo; k++) {
    //            for (int l = 0; l < nmo; l++) {
    //                en += 1.0 * D1[i*nmo+k]*D1[j*nmo+l] * V2[i*nmo*nmo*nmo+k*nmo*nmo+j*nmo+l];
    //                en -= 0.5 * D1[i*nmo+k]*D1[j*nmo+l] * V2[i*nmo*nmo*nmo+l*nmo*nmo+k*nmo+j];
    //            }
    //        }
    //    }
    //}
    en *= 2.0;

    // derivative of energy:
    /*for (int i = 0; i < nmo; i++) {
        for (int j = 0; j < nmo; j++) {
            for (int alpha = 0; alpha < nmo; alpha++) {
                dR1d[alpha*nmo+i] += 2.0 * R1d[alpha*nmo+j] * V1[i*nmo+j];
                dR1d[alpha*nmo+j] += 2.0 * R1d[alpha*nmo+i] * V1[i*nmo+j];
            }
            for (int k = 0; k < nmo; k++) {
                for (int l = 0; l < nmo; l++) {
                    for (int alpha = 0; alpha < nmo; alpha++) {
                        for (int beta = 0; beta < nmo; beta++) {
                            dR1d[alpha*nmo+i] += 2.0 * 1.0 * R1d[alpha*nmo+k]*R1d[beta*nmo+j]*R1d[beta*nmo+l] * V2[i*nmo*nmo*nmo+k*nmo*nmo+j*nmo+l];
                            dR1d[alpha*nmo+k] += 2.0 * 1.0 * R1d[alpha*nmo+i]*R1d[beta*nmo+j]*R1d[beta*nmo+l] * V2[i*nmo*nmo*nmo+k*nmo*nmo+j*nmo+l];
                            dR1d[beta*nmo+j]  += 2.0 * 1.0 * R1d[alpha*nmo+i]*R1d[alpha*nmo+k]*R1d[beta*nmo+l] * V2[i*nmo*nmo*nmo+k*nmo*nmo+j*nmo+l];
                            dR1d[beta*nmo+l]  += 2.0 * 1.0 * R1d[alpha*nmo+i]*R1d[alpha*nmo+k]*R1d[beta*nmo+j] * V2[i*nmo*nmo*nmo+k*nmo*nmo+j*nmo+l];

                            dR1d[alpha*nmo+i] -= 2.0 * 0.5 * R1d[alpha*nmo+k]*R1d[beta*nmo+j]*R1d[beta*nmo+l] * V2[i*nmo*nmo*nmo+l*nmo*nmo+k*nmo+j];
                            dR1d[alpha*nmo+k] -= 2.0 * 0.5 * R1d[alpha*nmo+i]*R1d[beta*nmo+j]*R1d[beta*nmo+l] * V2[i*nmo*nmo*nmo+l*nmo*nmo+k*nmo+j];
                            dR1d[beta*nmo+j]  -= 2.0 * 0.5 * R1d[alpha*nmo+i]*R1d[alpha*nmo+k]*R1d[beta*nmo+l] * V2[i*nmo*nmo*nmo+l*nmo*nmo+k*nmo+j];
                            dR1d[beta*nmo+l]  -= 2.0 * 0.5 * R1d[alpha*nmo+i]*R1d[alpha*nmo+k]*R1d[beta*nmo+j] * V2[i*nmo*nmo*nmo+l*nmo*nmo+k*nmo+j];
                        }
                    }
                }
            }
        }
    }*/
    double * I = (double*)malloc(nmo*nmo*sizeof(double));
    for (int i = 0; i < nmo; i++) {
        for (int m = 0; m < nmo; m++) {
            double dum = 0.0;
            for (int k = 0; k < nmo; k++) {
                for (int l = 0; l < nmo; l++) {
                    dum += D1[k*nmo+l] * ( 4.0 * V2[m*nmo*nmo*nmo+i*nmo*nmo+k*nmo+l]
                                               - V2[m*nmo*nmo*nmo+l*nmo*nmo+i*nmo+k]
                                               - V2[i*nmo*nmo*nmo+l*nmo*nmo+m*nmo+k] );
                }
                //for (int l = k+1; l < nmo; l++) {
                //    dum += 2.0 * (D1[k*nmo+l] * ( 4.0 * V2[m*nmo*nmo*nmo+i*nmo*nmo+k*nmo+l]
                //                               - V2[m*nmo*nmo*nmo+l*nmo*nmo+i*nmo+k]
                //                               - V2[i*nmo*nmo*nmo+l*nmo*nmo+m*nmo+k] ));
                //}
                //dum += D1[k*nmo+k] * ( 4.0 * V2[m*nmo*nmo*nmo+i*nmo*nmo+k*nmo+k]
                //                           - V2[m*nmo*nmo*nmo+k*nmo*nmo+i*nmo+k]
                //                           - V2[i*nmo*nmo*nmo+k*nmo*nmo+m*nmo+k] );
            }
            I[i*nmo+m] = dum + 2.0 * V1[i*nmo+m];
        }
    }
    /*for (int m = 0; m < nmo; m++) {
        for (int alpha = 0; alpha < nmo; alpha++) {
            double dum = 0.0;
            for (int i = 0; i < nmo; i++) {
                dum += R1d[alpha*nmo+i] * I[i*nmo+m];//( 2.0 * V1[i*nmo+m] + I[i*nmo+m] );
            }
            dR1d[alpha*nmo+m] = 2.0 * dum;
        }
    }*/
    F_DGEMM('n','n',nmo,nmo,nmo,2.0,I,nmo,R1d,nmo,0.0,dR1d,nmo);
    free (I);
    
    return en;
}

void SDPSolver::K2(){
    V2 = (double*)malloc(nmo*nmo*nmo*nmo*sizeof(double));
    V1 = (double*)malloc(nmo*nmo*sizeof(double));
    D1 = (double*)malloc(nmo*nmo*sizeof(double));
    Q1 = (double*)malloc(nmo*nmo*sizeof(double));

    boost::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));
    boost::shared_ptr<Matrix> tei = mints->mo_eri(Ca_,Ca_);
    double ** teip = tei->pointer();
    F_DCOPY(nmo*nmo*nmo*nmo,&teip[0][0],1,V2,1);

    // one-electron integrals in ao basis
    boost::shared_ptr<Matrix> T = mints->so_kinetic();
    boost::shared_ptr<Matrix> V = mints->so_potential();
    V->add(T);
    double**trans = Ca_->pointer();

    double ** oei = V->pointer();
    int full = Ca_->colspi()[0];
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],oei,trans);
    F_DCOPY(full*full,&oei[0][0],1,V1,1);
}

// so->mo transformation for 1-body matrix
void SDPSolver::SoToMo(int nsotemp,int nmotemp,double**mat,double**trans){
  double*tmp = (double*)malloc(sizeof(double)*nsotemp*nsotemp);
  int i,j,k;
  double dum;
  for (i=0; i<nsotemp; i++){
      for (j=0; j<nmotemp; j++){
          dum = 0.;
          for (k=0; k<nsotemp; k++){
              dum += mat[i][k]*trans[k][j];
          }
          tmp[i*nmotemp+j]=dum;
      }
  }
  for (i=0; i<nmotemp; i++){
      for (j=0; j<nmotemp; j++){
          dum = 0.;
          for (k=0; k<nsotemp; k++){
              dum += tmp[k*nmotemp+j]*trans[k][i];
          }
          mat[i][j]=dum;
      }
  }
  //F_DGEMM('n','n',nmotemp,nsotemp,nsotemp,1.0,&trans[0][0],nmotemp,&mat[0][0],nsotemp,0.0,tmp,nmotemp);
  //F_DGEMM('n','t',nmotemp,nmotemp,nsotemp,1.0,tmp,nmotemp,&trans[0][0],nmotemp,0.0,&mat[0][0],nsotemp);
  free(tmp);
}

void SDPSolver::common_init(){

    escf    = reference_wavefunction_->reference_energy();
    doccpi_ = reference_wavefunction_->doccpi();
    soccpi_ = reference_wavefunction_->soccpi();
    frzcpi_ = reference_wavefunction_->frzcpi();
    frzvpi_ = reference_wavefunction_->frzvpi();
    nmopi_  = reference_wavefunction_->nmopi();

    Da_ = SharedMatrix(reference_wavefunction_->Da());
    Ca_ = SharedMatrix(reference_wavefunction_->Ca());
    Fa_ = SharedMatrix(reference_wavefunction_->Fa());
    Db_ = SharedMatrix(reference_wavefunction_->Db());
    Cb_ = SharedMatrix(reference_wavefunction_->Cb());
    Fb_ = SharedMatrix(reference_wavefunction_->Fb());

    epsilon_a_= boost::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_= boost::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());
    nalpha_ = reference_wavefunction_->nalpha();
    nbeta_  = reference_wavefunction_->nbeta();

    nso = nmo = ndocc = nvirt = nfzc = nfzv = 0;
    for (int h=0; h<nirrep_; h++){
        nfzc   += frzcpi_[h];
        nfzv   += frzvpi_[h];
        nso    += nsopi_[h];
        nmo    += nmopi_[h]-frzcpi_[h]-frzvpi_[h];
        ndocc  += doccpi_[h];
    }
    ndoccact = ndocc - nfzc;
    nvirt    = nmo - ndoccact;

    // memory is from process::environment
    memory = Process::environment.get_memory();

    // set the wavefunction name
    name_ = "DQ-SCF";

    // orbital energies
    // qt
    eps = (double*)malloc(2*nmo*sizeof(double));
    int count=0;
    for (int h=0; h<nirrep_; h++){
        for (int norb = frzcpi_[h]; norb<doccpi_[h]; norb++){
            eps[count++] = epsilon_a_->get(h,norb);
        }
    }
    for (int h=0; h<nirrep_; h++){
        for (int norb = doccpi_[h]; norb<nmopi_[h]-frzvpi_[h]; norb++){
            eps[count++] = epsilon_a_->get(h,norb);
        }
    }

    trd1a = nalpha_;
    trd1b = nalpha_;

    //
    mu = 0.1;
    nconstraints = 1 + nmo*nmo;

    lambda = (double*)malloc(nconstraints*sizeof(double));
    c      = (double*)malloc(nconstraints*sizeof(double));
    cval   = (double*)malloc(nconstraints*sizeof(double));
    cerror = (double*)malloc(nconstraints*sizeof(double));

    for (int i = 0; i < nconstraints; i++) lambda[i] = 0.0;
    for (int i = 0; i < nconstraints; i++) c[i]      = 0.0;
    for (int i = 0; i < nconstraints; i++) cval[i]   = 0.0;
    for (int i = 0; i < nconstraints; i++) cerror[i] = 0.0;

    // trace condition:
    c[0] = nalpha_;
}

}}
