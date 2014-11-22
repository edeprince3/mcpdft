#include"psi4-dec.h"
#include<liboptions/liboptions.h>
#include<libmints/mints.h>
#include<libpsio/psio.h>
#include<physconst.h>
#include<libqt/qt.h>
#include<../bin/fnocc/blas.h>
#include<../bin/fnocc/frozen_natural_orbitals.h>
#include<psifiles.h>

#include"tdhf.h"

// boost numerical integrators live here:
#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>

void rk4_call( state_type &x , state_type &dxdt , double t ){
    psi::tdhf_cqed::MyTDHF->rk4_call_gah(x,dxdt,t);
}

#ifdef _OPENMP
    #include<omp.h>
#endif

using namespace psi;
using namespace boost;
using namespace fnocc;

namespace psi{ namespace tdhf_cqed {

TDHF::TDHF(boost::shared_ptr<Wavefunction> reference_wavefunction,Options & options):
  Wavefunction(options,_default_psio_lib_)
{
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}
TDHF::~TDHF(){
}
void TDHF::common_init(){

    escf    = reference_wavefunction_->reference_energy();
    doccpi_ = reference_wavefunction_->doccpi();
    soccpi_ = reference_wavefunction_->soccpi();
    frzcpi_ = reference_wavefunction_->frzcpi();
    frzvpi_ = reference_wavefunction_->frzvpi();
    nmopi_  = reference_wavefunction_->nmopi();

    Da_ = SharedMatrix(reference_wavefunction_->Da());
    Ca_ = SharedMatrix(reference_wavefunction_->Ca());
    Fa_ = SharedMatrix(reference_wavefunction_->Fa());

    epsilon_a_= boost::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
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

    if ( nfzc > 0 ) {
        throw PsiException("TDHF does not work with frozen core (yet).",__FILE__,__LINE__);
    }
    // memory is from process::environment
    memory = Process::environment.get_memory();

    // set the wavefunction name
    name_ = "TDHF";

    // orbital energies
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
    
    boost::shared_ptr<MintsHelper> mints (new MintsHelper());
    T   = mints->so_kinetic();
    V   = mints->so_potential();

    TransformIntegrals();

    //eri = mints->mo_eri(Ca_,Ca_);
    //for (int i = 0; i < nmo; i++) {
    //    for (int j = 0; j < nmo; j++) {
    //        for (int k = 0; k < nmo; k++) {
    //            for (int l = 0; l < nmo; l++) {
    //                double dum = 0.0;
    //                for (int q = 0; q < nQ; q++) {
    //                    dum += Qmo[q*nmo*nmo+i*nmo+j] * Qmo[q*nmo*nmo+k*nmo+l];
    //                }
    //                if (fabs(dum - eri->pointer()[i*nmo+j][k*nmo+l]) > 1e-6) {
    //                    printf("%5i %5i %5i %5i %20.12lf %20.12lf %20.12lf\n",i,j,k,l,dum , eri->pointer()[i*nmo+j][k*nmo+l],dum - eri->pointer()[i*nmo+j][k*nmo+l]);
    //                }
    //            }
    //        }
    //    }
    //}

    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],T->pointer(),Ca_->pointer());
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],V->pointer(),Ca_->pointer());

    Dre_init   = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Dim_init   = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Dre        = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Dim        = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Vext       = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Fre        = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Fim        = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Fre_init   = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Fim_init   = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    for (int i = 0; i < ndocc; i++) {
        Dre->pointer()[i][i] = 1.0;
    }

    // get dipole integrals:
    dipole = mints->so_dipole();
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],dipole[0]->pointer(),Ca_->pointer());
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],dipole[1]->pointer(),Ca_->pointer());
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],dipole[2]->pointer(),Ca_->pointer());

    // get polarization:
    polarization = (double*)malloc(sizeof(double)*3);
    if (options_["POLARIZATION"].has_changed()){
       if (options_["POLARIZATION"].size() != 3)
          throw PsiException("The POLARIZATION array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) polarization[i] = options_["POLARIZATION"][i].to_double();
    }else{
       polarization[0] = 1.0;
       polarization[1] = 1.0;
       polarization[2] = 1.0;
    }

    BuildFock(&(Dre->pointer())[0][0],&(Dim->pointer())[0][0],0.0);
    Dre_init->copy(Dre);
    Dim_init->copy(Dim);
    Fre_init->copy(Fre);
    Fim_init->copy(Fim);

    total_time = options_.get_double("TOTAL_TIME");
    time_step  = options_.get_double("TIME_STEP");
    laser_amp  = options_.get_double("LASER_AMP");
    laser_freq = options_.get_double("LASER_FREQ");
    transition_dpm = options_.get_double("LASER_TDPM");
    laser_time = options_.get_double("LASER_TIME");
    total_iter = total_time / time_step + 1;

    // which pulse shape do we want?
    if (options_.get_str("LASER_SHAPE") == "SIN_SQUARED") {
        // from prl:
        pulse_shape_ = 0;
    }else if (options_.get_str("LASER_SHAPE") == "TRAPEZOID") {
        // from 2007 schlegel paper (jcp 126, 244110 (2007))
        pulse_shape_ = 1;
    }else if (options_.get_str("LASER_SHAPE") == "PI_PULSE") {
        // pi pulse from licn paper
        pulse_shape_ = 2;
    }else if (options_.get_str("LASER_SHAPE") == "CONTINUOUS") {
        // continuous wave for rabi flopping
        pulse_shape_ = 3;
    }

    // if doing linear response, act on reference with dipole operator
    linear_response = false;
    if ( options_.get_bool("LINEAR_RESPONSE") ){
       c1 = (double*)malloc(ndocc*nvirt*sizeof(double));
       linear_response = true;
       throw PsiException ("Linear response not implemented for tdhf cqed",__FILE__, __LINE__);
       //for (int i = 0; i < ndocc; i++) {
       //    for (int a = 0; a < nvirt; a++) {
       //        Dre->pointer()[a+ndocc][i]  = dipole[0]->pointer()[i][a+ndocc];
       //        Dre->pointer()[a+ndocc][i] += dipole[1]->pointer()[i][a+ndocc];
       //        Dre->pointer()[a+ndocc][i] += dipole[2]->pointer()[i][a+ndocc];
       //    }
       //}
       double c0 = 0.0;
       for (int i = 0; i < ndocc; i++) {
           c0 += dipole[0]->pointer()[i][i];
           c0 += dipole[1]->pointer()[i][i];
           c0 += dipole[2]->pointer()[i][i];
       }
       for (int i = 0; i < ndocc; i++) {
           for (int a = 0; a < nvirt; a++) {
               c1[a*ndocc+i]  = dipole[0]->pointer()[i][a+ndocc];
               c1[a*ndocc+i] += dipole[1]->pointer()[i][a+ndocc];
               c1[a*ndocc+i] += dipole[2]->pointer()[i][a+ndocc];
           }
       }
       for (int i = 0; i < ndocc; i++) {
           for (int j = 0; j < ndocc; j++) {
               double dum = 0.0;
               for (int c = 0; c < nvirt; c++) {
                   dum += c1[c*ndocc+i] * c1[c*ndocc+j];
               }
               Dre->pointer()[i][j] -= dum;
           }
       }
       for (int a = 0; a < nvirt; a++) {
           for (int b = 0; b < nvirt; b++) {
               double dum = 0.0;
               for (int k = 0; k < ndocc; k++) {
                   dum += c1[a*ndocc+k] * c1[b*ndocc+k];
               }
               Dre->pointer()[a+ndocc][b+ndocc] += dum;
           }
       }
       for (int i = 0; i < ndocc; i++) {
           for (int a = 0; a < nvirt; a++) {
               Dre->pointer()[i][a+ndocc] = c0 * c1[a*ndocc+i];
               Dre->pointer()[a+ndocc][i] = c0 * c1[a*ndocc+i];
           }
       }
       Dre_init->copy(Dre);
       free(c1);
    }

    // pad the correlation function with zeros just to get more output points
    //extra_pts = 4*(ttot/time_step+2);
    extra_pts = 1000000;
    // correlation function or dipole acceleration (fourier transformed)
    midpt = total_time/time_step+extra_pts + 1;
    corr_func = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+2*total_time/time_step+2+2*extra_pts));
  
    // maximum frequency to output (eV)
    max_freq = 200.0;

    // stencil for second derivative of dipole moment
    stencil = (double*)malloc(sizeof(double)*5);
    memset((void*)stencil,'\0',5*sizeof(double));

    // dipole potential integrals:
    DipolePotentialIntegrals();

}

void TDHF::BuildFock(double * Dre, double * Dim, double curtime) {

    //for (int i = 0; i < nmo; i++) {
    //    for (int j = 0; j < nmo; j++) {
    //        double dumre = 0.0;
    //        double dumim = 0.0;
    //        for (int a = 0; a < nmo; a++) {
    //            for (int b = 0; b < nmo; b++) {
    //                dumre += (2.0 * eri->pointer()[i*nmo+j][a*nmo+b] - eri->pointer()[i*nmo+a][j*nmo+b]) * (Dre[a*nmo+b]);
    //                dumim += (2.0 * eri->pointer()[i*nmo+j][a*nmo+b] - eri->pointer()[i*nmo+a][j*nmo+b]) * (Dim[a*nmo+b]);
    //            }
    //        }
    //        dumre += T->pointer()[i][j] + V->pointer()[i][j];
    //        Fre->pointer()[i][j] = dumre;
    //        Fim->pointer()[i][j] = dumim;
    //    }
    //}


    for (int q = 0; q < nQ; q++) {
        double dumre = 0.0;
        double dumim = 0.0;
        for (int a = 0; a < nmo; a++) {
            for (int b = 0; b < nmo; b++) {
                dumre += 2.0 * Qmo[q*nmo*nmo+a*nmo+b] * Dre[a*nmo+b];
                dumim += 2.0 * Qmo[q*nmo*nmo+a*nmo+b] * Dim[a*nmo+b];
            }
        }
        Ire[q] = dumre;
        Iim[q] = dumim;
    }
    for (int i = 0; i < nmo; i++) {
        for (int j = 0; j < nmo; j++) {
            double dumre = 0.0;
            double dumim = 0.0;
            for (int q = 0; q < nQ; q++) {
                dumre += Qmo[q*nmo*nmo+i*nmo+j] * Ire[q];
                dumim += Qmo[q*nmo*nmo+i*nmo+j] * Iim[q];
            }
            dumre += T->pointer()[i][j] + V->pointer()[i][j];
            Fre->pointer()[i][j] = dumre;
            Fim->pointer()[i][j] = dumim;
        }
    }

    // K(re)
    F_DGEMM('t','n',nmo,nmo*nQ,nmo,1.0,Dre,nmo,Qmo,nmo,0.0,Ire,nmo);
    for (int q = 0; q < nQ; q++) {
        for (int a = 0; a < nmo; a++) {
            for (int j = 0; j < nmo; j++) {
                Iim[a*nmo+j] = Ire[q*nmo*nmo+j*nmo+a];
            }
        }
        F_DCOPY(nmo*nmo,Iim,1,Ire+q*nmo*nmo,1);
    }
    F_DGEMM('n','t',nmo,nmo,nQ*nmo,-1.0,Ire,nmo,Qmo,nmo,1.0,&(Fre->pointer()[0][0]),nmo);
    // K(im)
    F_DGEMM('t','n',nmo,nmo*nQ,nmo,1.0,Dim,nmo,Qmo,nmo,0.0,Ire,nmo);
    for (int q = 0; q < nQ; q++) {
        for (int a = 0; a < nmo; a++) {
            for (int j = 0; j < nmo; j++) {
                Iim[a*nmo+j] = Ire[q*nmo*nmo+j*nmo+a];
            }
        }
        F_DCOPY(nmo*nmo,Iim,1,Ire+q*nmo*nmo,1);
    }
    F_DGEMM('n','t',nmo,nmo,nQ*nmo,-1.0,Ire,nmo,Qmo,nmo,1.0,&(Fim->pointer()[0][0]),nmo);

    //for (int i = 0; i < nmo; i++) {
    //    for (int j = 0; j < nmo; j++) {
    //        double dumre = 0.0;
    //        double dumim = 0.0;
    //        for (int a = 0; a < nmo; a++) {
    //            for (int b = 0; b < nmo; b++) {
    //                for (int q = 0; q < nQ; q++) {
    //                    //dumre += (2.0 * Qmo[q*nmo*nmo+a*nmo+b]*Qmo[q*nmo*nmo+i*nmo+j] - Qmo[q*nmo*nmo+a*nmo+i]*Qmo[q*nmo*nmo+b*nmo+j]) * (Dre[a*nmo+b]);
    //                    //dumim += (2.0 * Qmo[q*nmo*nmo+a*nmo+b]*Qmo[q*nmo*nmo+i*nmo+j] - Qmo[q*nmo*nmo+a*nmo+i]*Qmo[q*nmo*nmo+b*nmo+j]) * (Dim[a*nmo+b]);
    //                    dumre += (- Qmo[q*nmo*nmo+a*nmo+i]*Qmo[q*nmo*nmo+b*nmo+j]) * (Dre[a*nmo+b]);
    //                    dumim += (- Qmo[q*nmo*nmo+a*nmo+i]*Qmo[q*nmo*nmo+b*nmo+j]) * (Dim[a*nmo+b]);
    //                }
    //            }
    //        }
    //        //Fre->pointer()[i][j] += dumre;
    //        Fim->pointer()[i][j] += dumim;
    //    }
    //}

    Vext->zero();
    double sigma = laser_time*0.5;
    if (!linear_response) {

        // add external field

        double fxn = 0.0;

        if (pulse_shape_ == 0 ) {
    
            // from prl:
            if ( curtime < laser_time ) {
                fxn = sin(M_PI*curtime/(laser_time));
                fxn *= fxn*laser_amp*sin(laser_freq*curtime);
            }
    
        } else if ( pulse_shape_ == 1 ) {
    
            // from 2007 schlegel paper (jcp 126, 244110 (2007))
            if (curtime <= 2.0 * M_PI / laser_freq)      fxn = laser_freq * curtime / (2.0 * M_PI) * laser_amp;
            else if (curtime <= 4.0 * M_PI / laser_freq) fxn = laser_amp;
            else if (curtime <= 6.0 * M_PI / laser_freq) fxn = (3.0 - laser_freq * curtime / (2.0 * M_PI) ) * laser_amp;
            fxn *= sin(laser_freq*curtime);
    
        } else if ( pulse_shape_ == 2 ) {
    
            // pi pulse from licn paper
            double sigma = laser_time*0.5;
            if ( curtime < laser_time ) {
                fxn = cos(M_PI*curtime/(2.0*sigma) + 0.5*M_PI);
                fxn *= M_PI/(sigma*transition_dpm) * fxn * laser_amp * cos(laser_freq*curtime);
            }
    
        } else if ( pulse_shape_ == 3 ) {
    
            // continuous wave for rabi flopping
            fxn = laser_amp*sin(laser_freq*curtime);
    
        }

        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                Fre->pointer()[i][j] += fxn * dipole[0]->pointer()[i][j] * polarization[0];
                Fre->pointer()[i][j] += fxn * dipole[1]->pointer()[i][j] * polarization[1];
                Fre->pointer()[i][j] += fxn * dipole[2]->pointer()[i][j] * polarization[2];
                //Vext->pointer()[i][j] += fxn * dipole[0]->pointer()[i][j];
                //Vext->pointer()[i][j] += fxn * dipole[1]->pointer()[i][j];
                //Vext->pointer()[i][j] += fxn * dipole[2]->pointer()[i][j];
            }
        }
    }
}

// so->mo transformation for 1-body matrix
void TDHF::SoToMo(int nsotemp,int nmotemp,double**mat,double**trans){
  double*tmp = (double*)malloc(sizeof(double)*nsotemp*nsotemp);
  F_DGEMM('n','n',nmotemp,nsotemp,nsotemp,1.0,&trans[0][0],nmotemp,&mat[0][0],nsotemp,0.0,&tmp[0],nmotemp);
  F_DGEMM('n','t',nmotemp,nmotemp,nsotemp,1.0,&tmp[0],nmotemp,&trans[0][0],nmotemp,0.0,&mat[0][0],nsotemp);
  free(tmp);
}

double TDHF::compute_energy() {

    // rk4
    state_type rk4_buffer(nmo*nS*nmo*nS*2);
    stepper_type rk4;
    fftw_iter   = 0;
    CorrelationFunction();

    double * factorB = (double*)malloc(sizeof(double)*5);
    double * factorb = (double*)malloc(sizeof(double)*5);

    // factors for symplectic integrator.  they come from that sanz-serna paper.
    double fac,tntf = 1./3924.;
    factorB[0] = (642.+sqrt(471.))*tntf;
    factorB[1] = 121.*(12.-sqrt(471.))*tntf;
    factorB[2] = 1. - 2.*(factorB[0]+factorB[1]);
    factorB[3] = factorB[1];
    factorB[4] = factorB[0];

    factorb[0] = 6./11.;
    factorb[1] = .5 - factorb[0];
    factorb[2] = factorb[1];
    factorb[3] = factorb[0];
    factorb[4] = 0.;


    for ( int iter = 0; iter < total_iter; iter++ ) {

        // skg's simplectic integrator:
        /*for (int m = 0; m < 5; m++) {

            // update imaginary part
            for (int i = 0; i < m; i++) fac += factorb[i];
            double faketime = iter * time_step + time_step * fac;
            BuildFock(&(Dre->pointer())[0][0],&(Dim->pointer())[0][0],faketime);
            for (int i = 0; i < nmo; i++) {
                for (int j = 0; j < nmo; j++) {
                    double dumim = 0.0;
                    for (int k = 0; k < nmo; k++) {

                        // full tdhf:
                        dumim += Fre->pointer()[i][k] * Dre->pointer()[k][j];
                        dumim -= Fim->pointer()[i][k] * Dim->pointer()[k][j];
                        dumim -= Dre->pointer()[i][k] * Fre->pointer()[k][j];
                        dumim += Dim->pointer()[i][k] * Fim->pointer()[k][j];

                        // linear response tdhf:
                     //   dumim += Fre_init->pointer()[i][k] * Dre->pointer()[k][j];
                     //   dumim -= Fim_init->pointer()[i][k] * Dim->pointer()[k][j];
                     //   dumim -= Dre->pointer()[i][k] * Fre_init->pointer()[k][j];
                     //   dumim += Dim->pointer()[i][k] * Fim_init->pointer()[k][j];

                     //   dumim += Fre->pointer()[i][k] * Dre_init->pointer()[k][j];
                     //   dumim -= Fim->pointer()[i][k] * Dim_init->pointer()[k][j];
                     //   dumim -= Dre_init->pointer()[i][k] * Fre->pointer()[k][j];
                     //   dumim += Dim_init->pointer()[i][k] * Fim->pointer()[k][j];

                     //   dumim += Vext->pointer()[i][k] * Dre_init->pointer()[k][j];
                     //   dumim -= Dre_init->pointer()[i][k] * Vext->pointer()[k][j];

                    }
                    tempim[i*nmo+j] = dumim;
                }
            }
            F_DAXPY(nmo*nmo,- factorB[m] * time_step,tempim,1,&(Dim->pointer())[0][0],1);

            // update real part
            fac = 0.0;
            for (int i = 0; i <= m; i++) fac += factorB[i];
            faketime = iter * time_step + time_step * fac;
            BuildFock(&(Dre->pointer())[0][0],&(Dim->pointer())[0][0],faketime);
            for (int i = 0; i < nmo; i++) {
                for (int j = 0; j < nmo; j++) {
                    double dumre = 0.0;
                    for (int k = 0; k < nmo; k++) {

                        // full tdhf:
                        dumre += Fre->pointer()[i][k] * Dim->pointer()[k][j];
                        dumre += Fim->pointer()[i][k] * Dre->pointer()[k][j];
                        dumre -= Dre->pointer()[i][k] * Fim->pointer()[k][j];
                        dumre -= Dim->pointer()[i][k] * Fre->pointer()[k][j];

                        // linear response tdhf:
                   //     dumre += Fre_init->pointer()[i][k] * Dim->pointer()[k][j];
                   //     dumre += Fim_init->pointer()[i][k] * Dre->pointer()[k][j];
                   //     dumre -= Dre->pointer()[i][k] * Fim_init->pointer()[k][j];
                   //     dumre -= Dim->pointer()[i][k] * Fre_init->pointer()[k][j];

                   //     dumre += Fre->pointer()[i][k] * Dim_init->pointer()[k][j];
                   //     dumre += Fim->pointer()[i][k] * Dre_init->pointer()[k][j];
                   //     dumre -= Dre_init->pointer()[i][k] * Fim->pointer()[k][j];
                   //     dumre -= Dim_init->pointer()[i][k] * Fre->pointer()[k][j];

                   //     dumre += Vext->pointer()[i][k] * Dim_init->pointer()[k][j];
                   //     dumre -= Dim_init->pointer()[i][k] * Vext->pointer()[k][j];

                    }
                    tempre[i*nmo+j] = dumre;
                }
            }
            F_DAXPY(nmo*nmo,factorb[m] * time_step,tempre,1,&(Dre->pointer())[0][0],1);
        }*/

        // RK4
        // y(n+1) = y( n ) + 1/6 h ( k1 + 2k2 + 2k3 + k4 )
        // t(n+1) = t( n ) + h

        for (int Ai = 0; Ai < nmo*nS; Ai++) {
            for (int Bj = 0; Bj < nmo*nS; Bj++) {
                rk4_buffer[Ai*nmo*nS+Bj] = Dre->pointer()[Ai][Bj];
            }
        }
        for (int Ai = 0; Ai < nmo*nS; Ai++) {
            for (int Bj = 0; Bj < nmo*nS; Bj++) {
                rk4_buffer[nmo*nS*nmo*nS + Ai*nmo*nS+Bj] = Dim->pointer()[Ai][Bj];
            }
        }
        rk4.do_step( rk4_call , rk4_buffer , iter * time_step , time_step );
        for (int Ai = 0; Ai < nmo*nS; Ai++) {
            for (int Bj = 0; Bj < nmo*nS; Bj++) {
                Dre->pointer()[Ai][Bj] = rk4_buffer[Ai*nmo*nS+Bj];
            }
        }
        for (int Ai = 0; Ai < nmo*nS; Ai++) {
            for (int Bj = 0; Bj < nmo*nS; Bj++) {
                Dim->pointer()[Ai][Bj] = rk4_buffer[nmo*nS*nmo*nS + Ai*nmo*nS+Bj];
            }
        }

        /*k1re->zero();
        k1im->zero();
        RK4(k1re,k1im,k1re,k1im,iter,0.0);
        RK4(k2re,k2im,k1re,k1im,iter,0.5);
        RK4(k3re,k3im,k2re,k2im,iter,0.5);
        RK4(k4re,k4im,k3re,k3im,iter,1.0);

        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(Dre->pointer()[0][0]),nmo,&(Fre->pointer()[0][0]),nmo,0.0,k1im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(Dim->pointer()[0][0]),nmo,&(Fim->pointer()[0][0]),nmo,1.0,k1im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(Fre->pointer()[0][0]),nmo,&(Dre->pointer()[0][0]),nmo,1.0,k1im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(Fim->pointer()[0][0]),nmo,&(Dim->pointer()[0][0]),nmo,1.0,k1im,nmo);

        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                tempre[i*nmo+j] = Dre->pointer()[i][j] + k1re[i*nmo+j] * time_step / 2.0;
                tempim[i*nmo+j] = Dim->pointer()[i][j] + k1im[i*nmo+j] * time_step / 2.0;
            }
        }
        // k2     = f( t( n + 1/2 h ) , y( n ) + h/2 k1 )
        BuildFock(tempre,tempim,(iter+0.5)*time_step);

        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(tempim[0]),nmo,&(Fre->pointer()[0][0]),nmo,0.0,k2re,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(tempre[0]),nmo,&(Fim->pointer()[0][0]),nmo,1.0,k2re,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(Fim->pointer()[0][0]),nmo,&(tempre[0]),nmo,1.0,k2re,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(Fre->pointer()[0][0]),nmo,&(tempim[0]),nmo,1.0,k2re,nmo);

        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(tempre[0]),nmo,&(Fre->pointer()[0][0]),nmo,0.0,k2im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(tempim[0]),nmo,&(Fim->pointer()[0][0]),nmo,1.0,k2im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(Fre->pointer()[0][0]),nmo,&(tempre[0]),nmo,1.0,k2im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(Fim->pointer()[0][0]),nmo,&(tempim[0]),nmo,1.0,k2im,nmo);

        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                tempre[i*nmo+j] = Dre->pointer()[i][j] + k2re[i*nmo+j] * time_step / 2.0;
                tempim[i*nmo+j] = Dim->pointer()[i][j] + k2im[i*nmo+j] * time_step / 2.0;
            }
        }*/

        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(tempim[0]),nmo,&(Fre->pointer()[0][0]),nmo,0.0,k3re,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(tempre[0]),nmo,&(Fim->pointer()[0][0]),nmo,1.0,k3re,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(Fim->pointer()[0][0]),nmo,&(tempre[0]),nmo,1.0,k3re,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(Fre->pointer()[0][0]),nmo,&(tempim[0]),nmo,1.0,k3re,nmo);

        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(tempre[0]),nmo,&(Fre->pointer()[0][0]),nmo,0.0,k3im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(tempim[0]),nmo,&(Fim->pointer()[0][0]),nmo,1.0,k3im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(Fre->pointer()[0][0]),nmo,&(tempre[0]),nmo,1.0,k3im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(Fim->pointer()[0][0]),nmo,&(tempim[0]),nmo,1.0,k3im,nmo);

        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                tempre[i*nmo+j] = Dre->pointer()[i][j] + k3re[i*nmo+j] * time_step;
                tempim[i*nmo+j] = Dim->pointer()[i][j] + k3im[i*nmo+j] * time_step;
            }
        }
        // k4     = f( t( n + h )     , y( n ) + h k3 )
        BuildFock(tempre,tempim,(iter+1)*time_step);

        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(tempim[0]),nmo,&(Fre->pointer()[0][0]),nmo,0.0,k4re,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(tempre[0]),nmo,&(Fim->pointer()[0][0]),nmo,1.0,k4re,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(Fim->pointer()[0][0]),nmo,&(tempre[0]),nmo,1.0,k4re,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(Fre->pointer()[0][0]),nmo,&(tempim[0]),nmo,1.0,k4re,nmo);

        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(tempre[0]),nmo,&(Fre->pointer()[0][0]),nmo,0.0,k4im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(tempim[0]),nmo,&(Fim->pointer()[0][0]),nmo,1.0,k4im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(Fre->pointer()[0][0]),nmo,&(tempre[0]),nmo,1.0,k4im,nmo);
        F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(Fim->pointer()[0][0]),nmo,&(tempim[0]),nmo,1.0,k4im,nmo);

        // y(n+1) = y( n ) + 1/6 h ( k1 + 2k2 + 2k3 + k4 )

        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                Dre->pointer()[i][j] += 1.0 / 6.0 * time_step * ( k1re[i*nmo+j] + 2.0 * k2re[i*nmo+j] + 2.0 * k3re[i*nmo+j] + k4re[i*nmo+j] );
                Dim->pointer()[i][j] += 1.0 / 6.0 * time_step * ( k1im[i*nmo+j] + 2.0 * k2im[i*nmo+j] + 2.0 * k3im[i*nmo+j] + k4im[i*nmo+j] );
            }
        }



        // evaluate dipole moment:
        double dpre = F_DDOT(nmo*nmo,&(Dre->pointer())[0][0],1,&(dipole[0]->pointer())[0][0],1);
        dpre       += F_DDOT(nmo*nmo,&(Dre->pointer())[0][0],1,&(dipole[1]->pointer())[0][0],1);
        dpre       += F_DDOT(nmo*nmo,&(Dre->pointer())[0][0],1,&(dipole[2]->pointer())[0][0],1);
        double dpim = F_DDOT(nmo*nmo,&(Dim->pointer())[0][0],1,&(dipole[0]->pointer())[0][0],1);
        dpim       += F_DDOT(nmo*nmo,&(Dim->pointer())[0][0],1,&(dipole[1]->pointer())[0][0],1);
        dpim       += F_DDOT(nmo*nmo,&(Dim->pointer())[0][0],1,&(dipole[2]->pointer())[0][0],1);
        stencil[0] = stencil[1];
        stencil[1] = stencil[2];
        stencil[2] = stencil[3];
        stencil[3] = stencil[4];
        stencil[4] = dpre;

        double dipole_acceleration = -1.0/12.0 * (stencil[0]+stencil[4])
                                   +  4.0/3.0  * (stencil[1]+stencil[3])
                                   -  5.0/2.0  *  stencil[2];
        dipole_acceleration /= (time_step*time_step);
        // start accumulating the dipole acceleration after the pulse is over
        if (!linear_response && iter * time_step >laser_time+0*time_step){
           corr_func[fftw_iter][0] = dipole_acceleration;
           corr_func[fftw_iter][1] = 0.0;
           fftw_iter++;
        }

        CorrelationFunction();

        double en = F_DDOT(nmo*nmo,&Fre->pointer()[0][0],1,&Dre->pointer()[0][0],1)
                  + F_DDOT(nmo*nmo,&Fim->pointer()[0][0],1,&Dim->pointer()[0][0],1)
                  + F_DDOT(nmo*nmo,&Dre->pointer()[0][0],1,&T->pointer()[0][0],1)
                  + F_DDOT(nmo*nmo,&Dre->pointer()[0][0],1,&V->pointer()[0][0],1);
        en += Process::environment.molecule()->nuclear_repulsion_energy();
        //printf("%20.12lf\n",en+ Process::environment.molecule()->nuclear_repulsion_energy());

        //fprintf(outfile,"@TDCI TIME %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",iter*time_step,dpre,dipole_acceleration,cf_real,cf_imag);
        fprintf(outfile,"@TDCI TIME %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",iter*time_step,dpre,dipole_acceleration,en,cf_imag);
        //fprintf(outfile,"@TDCI TIME %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",iter*time_step,dpre,dipole_acceleration,Dre->pointer()[0][0],Dre->pointer()[1][1],en);

    }

    free(factorb);
    free(factorB);


    free(tempre);
    free(tempim);
    free(k1re);
    free(k2re);
    free(k3re);
    free(k4re);
    free(k1im);
    free(k2im);
    free(k3im);
    free(k4im);

    // fourier transform and spectrum
    FFTW();
    Spectrum();

    free(corr_func);
    

    return 0.0;
}
void TDHF::CorrelationFunction(){
  // the correlation function
  //cf_real = F_DDOT(nmo,&(Dinit->pointer())[0][0],nmo+1,&(Dre->pointer())[0][0],nmo+1);
  //cf_imag = F_DDOT(nmo,&(Dinit->pointer())[0][0],nmo+1,&(Dim->pointer())[0][0],nmo+1);
  cf_real = 0.0;
  cf_imag = 0.0;
  for (int i = 0; i < ndocc; i++) {
      for (int a = 0; a < nvirt; a++) {
          cf_real += Dre_init->pointer()[i][a+ndocc] * Dre->pointer()[i][a+ndocc];
          cf_imag += Dre_init->pointer()[i][a+ndocc] * Dim->pointer()[i][a+ndocc];
          //cf_real += c1[a*ndocc+i] * Dre->pointer()[i][a+ndocc];
          //cf_imag += c1[a*ndocc+i] * Dim->pointer()[i][a+ndocc];
          //cf_real += c1[a*ndocc+i] * Dre->pointer()[a+ndocc][i];
          //cf_imag += c1[a*ndocc+i] * Dim->pointer()[a+ndocc][i];
      }
  }
  //for (int i = 0; i < nmo; i++) {
  //    for (int j = 0; j < nmo; j++) {
  //        cf_real += Dinit->pointer()[i][j] * Dre->pointer()[i][j];
  //        cf_imag += Dinit->pointer()[i][j] * Dim->pointer()[i][j];
  //    }
  //}
  //if (linear_response){
  //   corr_func[midpt + fftw_iter][0] =  cf_real;
  //   corr_func[midpt + fftw_iter][1] =  cf_imag;
  //   corr_func[midpt - fftw_iter][0] =  cf_real;
  //   corr_func[midpt - fftw_iter][1] = -cf_imag;
  //   fftw_iter++;
  //}
}

void TDHF::rk4_call_gah( state_type &x , state_type &dxdt , double t ){
    for (int i = 0; i < nmo*nS; i++) {
        for (int j = 0; j < nmo*nS; j++) {
            tempre->pointer()[i][j] = x[i*nmo*nS+j];
            tempim->pointer()[i][j] = x[nmo*nmo*nS*nS + i*nmo*nS+j];
        }
    }
        
    // kout = f( t( n + mh ) , y( n ) + mh kin)
    ExtField(t); 
    ElectronicContribution(tempre,tempim,k1re,k1im);
    PlasmonContribution(tempre,tempim,k1re,k1im);
    InteractionContribution(tempre,tempim,k1re,k1im,dpre,dipole_p);
    BuildLindblad(tempre,tempim);
    k1re->add(Lre);
    k1im->add(Lim);

    for (int i = 0; i < nmo*nS; i++) {
        for (int j = 0; j < nmo*nS; j++) {
            dxdt[i*nmo*nS+j]                 = tempre->pointer()[i][j];
            dxdt[nmo*nmo*nS*nS + i*nmo*nS+j] = tempim->pointer()[i][j];
        }
    }
}

void TDHF::RK4(boost::shared_ptr<Matrix> koutre, boost::shared_ptr<Matrix>koutim,
               boost::shared_ptr<Matrix> kinre, boost::shared_ptr<Matrix>kinim, int iter, double step) {

    for (int i = 0; i < nmo*nS; i++) {
        for (int j = 0; j < nmo*nS; j++) {
            tempre->pointer()[i][j] = Dre->pointer()[i][j] + kinre->pointer()[i][j] * step*time_step;
            tempim->pointer()[i][j] = Dim->pointer()[i][j] + kinim->pointer()[i][j] * step*time_step;
        }
    }
    for (int i = 0; i < nmo*nS; i++) {
        for (int j = 0; j < nmo*nS; j++) {
            double val1 = tempre->pointer()[i][j];
            double val2 = tempre->pointer()[j][i];
            tempre->pointer()[i][j] = 0.5 * (val1 + val2);
            tempre->pointer()[j][i] = 0.5 * (val1 + val2);
            val1 = tempim->pointer()[i][j];
            val2 = tempim->pointer()[j][i];
            tempim->pointer()[i][j] =  0.5 * (val1 - val2);
            tempim->pointer()[j][i] = -0.5 * (val1 - val2);
        }
    }
    /*double trace = 0.0;
    for (int i = 0; i < nmo*nS; i++) {
        trace += tempre->pointer()[i][i];
    }
    for (int i = 0; i < nmo*nS; i++) {
        for (int j = 0; j < nmo*nS; j++) {
            tempre->pointer()[i][j] *= nalpha_/trace; 
            tempim->pointer()[i][j] *= nalpha_/trace; 
        }
    }*/
        
    // kout = f( t( n + mh ) , y( n ) + mh kin)
    ExtField(iter*time_step + step*time_step); 
    ElectronicContribution(tempre,tempim,koutre,koutim);
    PlasmonContribution(tempre,tempim,koutre,koutim);
    InteractionContribution(tempre,tempim,koutre,koutim,dpre,dipole_p);
    BuildLindblad(tempre,tempim);
    koutre->add(Lre);
    koutim->add(Lim);

}
void TDHF::FFTW(){
  if (linear_response) {
      fftw_plan p;
      for (long int i=0; i<extra_pts; i++){
          corr_func[midpt + fftw_iter+i][0] = 0.0;
          corr_func[midpt + fftw_iter+i][1] = 0.0;
          corr_func[midpt - fftw_iter-i][0] = 0.0;
          corr_func[midpt - fftw_iter-i][1] = 0.0;
      }
      p = fftw_plan_dft_1d((int)(2*extra_pts+2*fftw_iter+1),corr_func,corr_func,FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
  }else {
      fftw_plan p;
      for (long int i=0; i<extra_pts; i++){
          corr_func[fftw_iter+i][0] = 0.0;
          corr_func[fftw_iter+i][1] = 0.0;
      }
      p = fftw_plan_dft_1d((int)(extra_pts+fftw_iter),corr_func,corr_func,FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
  }
}

// output absorption spectrum
void TDHF::Spectrum(){
  int i;
  double val,valr,vali,twopi = 2.0*M_PI;
  double w;
  if (linear_response){
     fprintf(outfile,"\n");
     fprintf(outfile,"        ************************************************************\n");
     fprintf(outfile,"        *                                                          *\n");
     fprintf(outfile,"        *        Absorption spectrum (oscillator strength)         *\n");
     fprintf(outfile,"        *         as computed by the Fourier transform of          *\n");
     fprintf(outfile,"        *           the dipole autocorrelation function            *\n");
     fprintf(outfile,"        *                                                          *\n");
     fprintf(outfile,"        *      I(w) = 2/3 w FourierTransform ( <M(0)|M(t)> )       *\n");
     fprintf(outfile,"        *                                                          *\n");
     fprintf(outfile,"        ************************************************************\n");
     fprintf(outfile,"\n");
     fprintf(outfile,"                                w(eV)");
     fprintf(outfile,"                 I(w)\n");
//double integral = 0.0;
     for (i=0; i<(int)(1+2*fftw_iter+2*extra_pts); i++){
         w    = twopi*i/((2*extra_pts+1+2*fftw_iter)*time_step);
         if (w*pc_hartree2ev>max_freq) break;
         valr = corr_func[i][0]/(fftw_iter);
         vali = corr_func[i][1]/(fftw_iter);
         val  = 2.0/3.0 * w * sqrt(valr*valr+vali*vali);
//integral += val;
         fprintf(outfile,"      @Frequency %20.12lf %20.12lf %20.12lf %20.12lf \n",w*pc_hartree2ev,val,2./3.*w*valr,2./3.*w*vali);
            //"abs %20.12lf %20.12lf %20.12lf\n",twopi*i/(fftw_iter*time_step)*pc_hartree2ev,
            //corr_func[i][0]/(fftw_iter),corr_func[i][1]/(fftw_iter));
     }
//printf("heyhey number of electrons %20.12lf\n",integral);// * w / (1+2*fftw_iter+2.*extra_pts) );
  }else{
     fprintf(outfile,"\n");
     fprintf(outfile,"        ***********************************************************\n");
     fprintf(outfile,"        *                                                         *\n");
     fprintf(outfile,"        *                    Emission spectrum                    *\n");
     fprintf(outfile,"        *         as computed by the Fourier transform of         *\n");
     fprintf(outfile,"        *                   dipole acceleration                   *\n");
     fprintf(outfile,"        *                                                         *\n");
     fprintf(outfile,"        *     I(w) = |FourierTransform ( d^2 D(t) / dt^2 )|^2     *\n");
     fprintf(outfile,"        *                                                         *\n");
     fprintf(outfile,"        ***********************************************************\n");
     fprintf(outfile,"\n");
     fprintf(outfile,"                                w(eV)");
     fprintf(outfile,"                 I(w)\n");
     for (i=1; i<(int)(fftw_iter+extra_pts); i++){
         w    = twopi*i/((extra_pts+fftw_iter)*time_step);
         if (w*pc_hartree2ev>max_freq) break;
         valr = corr_func[i][0]/fftw_iter;
         val  = valr*valr;
         fprintf(outfile,"      @Frequency %20.12lf %20.12lf \n",w*pc_hartree2ev,val);
     }
  }
}

void TDHF::TransformIntegrals() {

    long int o = ndoccact;
    long int v = nvirt;
    long int full = o+v+nfzc+nfzv;

    nQ = (int)Process::environment.globals["NAUX (CC)"];

    long int memory = Process::environment.get_memory();
    // subtract out 20 nmo*nmo + 250 MB to be sure we have enough room
    memory -= sizeof(double)* 20L * nmo * nmo - 250L * 1024L * 1024L;

    if ( memory < sizeof(double) * (2L*nmo*nmo*nQ) ) {
        throw PsiException("TDHF::TransformIntegrals: not enough memory",__FILE__,__LINE__);
    }

    Qmo = (double*)malloc(nmo*nmo*nQ*sizeof(double));
    memset((void*)Qmo,'\0',nmo*nmo*nQ*sizeof(double));

    Ire = (double*)malloc(nmo*nmo*nQ*sizeof(double));
    Iim = (double*)malloc( ( nmo * nmo > nQ ? nmo * nmo : nQ ) *sizeof(double));

    double ** Ca = Ca_->pointer();

    // available memory:
    memory -= sizeof(double) * (2L*nmo*nmo*nQ);
    int ndoubles = memory / sizeof(double) / 2;
    if ( nso * nso * nQ < ndoubles ) ndoubles = nso*nso*nQ;

    double * buf1 = (double*)malloc(ndoubles * sizeof(double));
    double * buf2 = (double*)malloc(ndoubles * sizeof(double));

    // (Q|rs)
    boost::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    psio_address addr1  = PSIO_ZERO;
    psio_address addrvo = PSIO_ZERO;
    long int nrows = 1;
    long int rowsize = nQ;
    while ( rowsize*nso*nso > ndoubles ) {
        nrows++;
        rowsize = nQ / nrows;
        if (nrows * rowsize < nQ) rowsize++;
        if (rowsize == 1) break;
    }
    long int lastrowsize = nQ - (nrows - 1L) * rowsize;
    long int * rowdims = new long int [nrows];
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qso CC",(char*)&buf1[0],rowdims[row]*nso*nso*sizeof(double),addr1,&addr1);
        F_DGEMM('n','n',full,nso*rowdims[row],nso,1.0,&Ca[0][0],full,buf1,nso,0.0,buf2,full);
        for (int q = 0; q < rowdims[row]; q++) {
            for (int mu = 0; mu < nso; mu++) {
                F_DCOPY(full,buf2+q*nso*full+mu*full,1,buf1+q*nso*full+mu,nso);
            }
        }
        F_DGEMM('n','n',full,full*rowdims[row],nso,1.0,&Ca[0][0],full,buf1,nso,0.0,buf2,full);

        // Qmo
        #pragma omp parallel for schedule (static)
        for (int q = 0; q < rowdims[row]; q++) {
            for (int a = 0; a < nmo; a++) {
                for (int b = 0; b < nmo; b++) {
                    Qmo[(q+rowdims[0]*row)*nmo*nmo+a*nmo+b] = buf2[q*full*full+a*full+b];
                }
            }
        }
    }
    free(buf1);
    free(buf2);
    delete rowdims;
    psio->close(PSIF_DCC_QSO,1);
}

}}
