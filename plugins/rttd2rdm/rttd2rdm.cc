/*
 *@BEGIN LICENSE
 *
 * rttd2rdm by Psi4 Developer, a plugin to:
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

#include "psi4/psi4-dec.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/physconst.h"
#include "psi4/libqt/qt.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libqt/qt.h"

#include "rttd2rdm.h"
#include "blas.h"
#include "blas_complex.h"

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() 0.0
    #define omp_get_max_threads() 1
#endif

using namespace std;
using namespace psi;
using namespace fnocc;

namespace psi{ namespace rttd2rdm {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "RTTD2RDM"|| options.read_globals()) {
        /*- total time in a.u. -*/
        options.add_double("TOTAL_TIME", 100.0);
        /*- time step in a.u. -*/
        options.add_double("TIME_STEP", 0.2);
        /*- pulse shape -*/
        options.add_str("LASER_SHAPE", "SIN_SQUARED", "SIN_SQUARED TRAPEZOID PI_PULSE CONTINUOUS GAUSSIAN");
        /*- transition dipole moment for pi pulse -*/
        options.add_double("LASER_TDPM", -0.415638122584);
        /*- amplitude of pulse in a.u.-*/
        options.add_double("LASER_AMP", 0.05);
        /*- frequency of pulse in a.u. (default is the 1 fs pulse) -*/
        options.add_double("LASER_FREQ", 0.1519829846);
        /*- width of pulse in a.u. (1 fs default) -*/
        options.add_double("LASER_TIME", 41.3414);
        /*- polarization (default x+y+z). -*/
        options.add("POLARIZATION",new ArrayType());
        /*- DF basis for SCF and TD2RDM -*/
        options.add_str("DF_BASIS_SCF", "");
    }

    return true;
}

extern "C"
SharedWavefunction rttd2rdm(SharedWavefunction ref_wfn, Options& options)
{

    std::shared_ptr<TD2RDM> rdm (new TD2RDM(ref_wfn,options));
    rdm->compute_energy();

    // just return original wave function for now
    return ref_wfn;
}

TD2RDM::TD2RDM(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options):
  Wavefunction(options)
{

    reference_wavefunction_ = reference_wavefunction;
    common_init();

}

TD2RDM::~TD2RDM()
{
}

void TD2RDM::common_init() {

    shallow_copy(reference_wavefunction_);

    escf_     = reference_wavefunction_->reference_energy();
    nalpha_   = reference_wavefunction_->nalpha();
    nbeta_    = reference_wavefunction_->nbeta();
    nalphapi_ = reference_wavefunction_->nalphapi();
    nbetapi_  = reference_wavefunction_->nbetapi();
    doccpi_   = reference_wavefunction_->doccpi();
    soccpi_   = reference_wavefunction_->soccpi();
    frzcpi_   = reference_wavefunction_->frzcpi();
    frzvpi_   = reference_wavefunction_->frzvpi();
    nmopi_    = reference_wavefunction_->nmopi();
    nirrep_   = reference_wavefunction_->nirrep();
    nso_      = reference_wavefunction_->nso();
    nmo_      = reference_wavefunction_->nmo();
    nsopi_    = reference_wavefunction_->nsopi();
    molecule_ = reference_wavefunction_->molecule();
    //enuc_     = molecule_->nuclear_repulsion_energy();

    if ( nirrep_ != 1 ) {
        throw PsiException("plugin RTTD2RDM only works with symmetry c1",__FILE__,__LINE__);
    }

    epsilon_a_= std::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());

    epsilon_b_= std::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());

    nalpha_ = reference_wavefunction_->nalpha();
    nbeta_  = reference_wavefunction_->nbeta();

    nso_ = nmo_ = ndocc_ = nvirt_ = nfrzc_ = nfrzv_ = 0;
    for (int h=0; h<nirrep_; h++){
        nfrzc_  += frzcpi_[h];
        nfrzv_  += frzvpi_[h];
        nso_    += nsopi_[h];
        nmo_    += nmopi_[h]-frzcpi_[h]-frzvpi_[h];
        ndocc_  += doccpi_[h];
    }
    ndoccact_ = ndocc_ - nfrzc_;
    nvirt_    = nmo_ - ndoccact_;

    if ( nirrep_ > 1 ) {
        throw PsiException("plugin rttd2rdm only works with c1 symmetry",__FILE__,__LINE__);
    }
    if ( nfrzc_ > 0 ) {
        throw PsiException("plugin rttd2rdm does not work with frozen core orbitals",__FILE__,__LINE__);
    }
    if ( nfrzv_ > 0 ) {
        throw PsiException("plugin rttd2rdm does not work with frozen virtual orbitals",__FILE__,__LINE__);
    }

    bool is_df  = ( options_.get_str("SCF_TYPE") == "DF" );
    if ( ! is_df ) {
        throw PsiException("Error: SDP SCF only works with density fitting for now",__FILE__,__LINE__);
    }

    // grab the molecule
    std::shared_ptr<Molecule> mol = Process::environment.molecule();

    // build primary basis:
    primary_ = basisset();

    // total number of basis functions
    nso_ = primary_->nbf();

    Dre_ = (std::shared_ptr<Matrix>)(new Matrix("real density matrix",nso_,nso_));
    Dim_ = (std::shared_ptr<Matrix>)(new Matrix("imag density matrix",nso_,nso_));
    for (int i = 0; i < ndocc_; i++) {
        Dre_->pointer()[i][i] = 1.0;
    }

    tempre_ = (std::shared_ptr<Matrix>)(new Matrix("real temporary matrix",nso_,nso_));
    tempim_ = (std::shared_ptr<Matrix>)(new Matrix("imag temporary matrix",nso_,nso_));

    k1re_ = (std::shared_ptr<Matrix>)(new Matrix("real k1 matrix",nso_,nso_));
    k1im_ = (std::shared_ptr<Matrix>)(new Matrix("imag k1 matrix",nso_,nso_));

    k2re_ = (std::shared_ptr<Matrix>)(new Matrix("real k2 matrix",nso_,nso_));
    k2im_ = (std::shared_ptr<Matrix>)(new Matrix("imag k2 matrix",nso_,nso_));

    k3re_ = (std::shared_ptr<Matrix>)(new Matrix("real k3 matrix",nso_,nso_));
    k3im_ = (std::shared_ptr<Matrix>)(new Matrix("imag k3 matrix",nso_,nso_));

    k4re_ = (std::shared_ptr<Matrix>)(new Matrix("real k4 matrix",nso_,nso_));
    k4im_ = (std::shared_ptr<Matrix>)(new Matrix("imag k4 matrix",nso_,nso_));

    kre_ = (std::shared_ptr<Matrix>)(new Matrix("real k matrix",nso_,nso_));
    kim_ = (std::shared_ptr<Matrix>)(new Matrix("imag k matrix",nso_,nso_));

    Fre_ = (std::shared_ptr<Matrix>)(new Matrix("real Fock matrix",nso_,nso_));
    Fim_ = (std::shared_ptr<Matrix>)(new Matrix("imag Fock matrix",nso_,nso_));

    C2re_ = (std::shared_ptr<Matrix>)(new Matrix("real 2-cumulant matrix",nso_*nso_,nso_*nso_));
    C2im_ = (std::shared_ptr<Matrix>)(new Matrix("imag 2-cumulant matrix",nso_*nso_,nso_*nso_));

    //Jre_ = (std::shared_ptr<Matrix>)(new Matrix("real coulomb matrix",nso_,nso_));
    //Jim_ = (std::shared_ptr<Matrix>)(new Matrix("imag coulomb matrix",nso_,nso_));

    //Kre_ = (std::shared_ptr<Matrix>)(new Matrix("real exchange matrix",nso_,nso_));
    //Kim_ = (std::shared_ptr<Matrix>)(new Matrix("imag exchange matrix",nso_,nso_));

    Da_ = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Da()));
    Db_ = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Db()));
    Ca_ = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Ca()));
    Cb_ = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Cb()));
    Fa_ = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Fa()));
    Fb_ = (std::shared_ptr<Matrix>)(new Matrix(reference_wavefunction_->Fb()));

    epsilon_a_ = (std::shared_ptr<Vector>)(new Vector("Alpha orbital energies",nso_));
    epsilon_b_ = (std::shared_ptr<Vector>)(new Vector("Beta orbital energies",nso_));

    // build auxiliary basis:
    auxiliary_ = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

    // potential stuff for DFT:
    //if ( is_dft_ ) {
    //    potential_ = VBase::build_V(options_,(options_.get_str("REFERENCE") == "RKS" ? "RV" : "UV"));
    //    potential_->initialize();
    //    functional_ = potential_->functional();

    //    // Print the KS-specific stuff
    //    potential_->print_header();
    //}

    // JK object
    jk_ = (std::shared_ptr<DFJK>)(new DFJK(primary_,auxiliary_));

    // memory for jk (TODO: check exactly what I will need)
    jk_->set_memory(0.8 * memory_);

    // integral cutoff
    jk_->set_cutoff(options_.get_double("INTS_TOLERANCE"));
   
    //if ( is_dft_ ) { 
    //    // Need a temporary functional
    //    std::shared_ptr<SuperFunctional> functional = SuperFunctional::current(options_);

    //    // K matrices
    //    jk_->set_do_K(functional->is_x_hybrid());
    //    // wK matrices
    //    jk_->set_do_wK(functional->is_x_lrc());
    //    // w Value
    //    jk_->set_omega(functional->x_omega());

    //}else {
        // Do J/K, Not wK 
        jk_->set_do_J(true);
        jk_->set_do_K(true);
        jk_->set_do_wK(false);
    //}

    // Initialize calls your derived class's preiterations member
    jk_->initialize();

    Integrals();

    // TD stuff:

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
    }else if (options_.get_str("LASER_SHAPE") == "GAUSSIAN"){
        // gaussian pulse
        pulse_shape_ = 4;
    }

    total_time_     = options_.get_double("TOTAL_TIME");
    time_step_      = options_.get_double("TIME_STEP");
    laser_amp_      = options_.get_double("LASER_AMP");
    laser_freq_     = options_.get_double("LASER_FREQ");
    transition_dpm_ = options_.get_double("LASER_TDPM");
    laser_time_     = options_.get_double("LASER_TIME");

    total_iter_     = total_time_ / time_step_ + 1L;

    polarization_ = (double*)malloc(sizeof(double)*3);
    memset((void*)polarization_,'\0',3*sizeof(double));
    if (options_["POLARIZATION"].has_changed()){
       if (options_["POLARIZATION"].size() != 3)
          throw PsiException("The POLARIZATION array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) polarization_[i] = options_["POLARIZATION"][i].to_double();
       int pol = 0;
       for (int i = 0; i < 3; i++) {
           if (fabs(polarization_[i]) > 1e-9) {
               pol++;
           }
       }
       if (pol > 1) {
          throw PsiException("only plane polarized light supported",__FILE__,__LINE__);
       }
    }else{
       polarization_[0] = 0.0;
       polarization_[1] = 0.0;
       polarization_[2] = 1.0;
    }


    // fftw
    corr_func_x     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+total_time_/time_step_+2));
    corr_func_y     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+total_time_/time_step_+2));
    corr_func_z     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+total_time_/time_step_+2));

}

void TD2RDM::Integrals() {

    // grab the one-electron integrals from MintsHelper:
    std::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));

    // one-electron kinetic energy integrals
    std::shared_ptr<Matrix> T = mints->so_kinetic();

    // one-electron potential energy integrals
    std::shared_ptr<Matrix> V = mints->so_potential();

    // build the core hamiltonian
    Hcore_ = (std::shared_ptr<Matrix>)(new Matrix(T));
    Hcore_->add(V);

    // transform 1-electron integrals to MO basis:
    Hcore_->transform(Ca_);

    // get dipole integrals in MO basis:
    std::vector<std::shared_ptr<Matrix> > dipole_ = mints->so_dipole();
    mux_ = (std::shared_ptr<Matrix>)(new Matrix(dipole_[0]));
    muy_ = (std::shared_ptr<Matrix>)(new Matrix(dipole_[1]));
    muz_ = (std::shared_ptr<Matrix>)(new Matrix(dipole_[2]));
    mux_->transform(Ca_);
    muy_->transform(Ca_);
    muz_->transform(Ca_);

    // get three-index integrals
    //printf("Get three-index integrals \n"); 
    nQ_ = auxiliary_->nbf();

        // use the molecule to determine the total number of electrons
    int charge     = molecule_->molecular_charge();
    int nelectron  = 0;
    for (int i = 0; i < molecule_->natom(); i++) {
        nelectron += (int)molecule_->Z(i);
    }
    nelectron -= charge;

    // this code only works for closed shells
    if ( nelectron % 2 != 0 ) {
        throw PsiException("plugin rtcc2 only works for closed shells",__FILE__,__LINE__);
    }

    // the number of doubly occupied orbtials (or alpha electrons)
    int na = nelectron / 2;

    // construct the three-index integrals 
    // similarly, the number of active vs inactive orbitals isn't really important here.
    std::shared_ptr<DFTensor> DF (new DFTensor(primary_,auxiliary_,Ca_,na,nso_-na,na,nso_-na,options_));
    Qmo_ = DF->Qmo();
    Qp_ = Qmo_->pointer();
    //printf("Three-index integrals built sucessfully!\n"); 
}

double TD2RDM::compute_energy() {

//    ReadTPDM();
//    for (long int p = 0; p < nso_; p++){
//        for (long int q = 0; q < nso_; q++){
//            Dre_->pointer()[p][q] = Da[p*nso_+q];
//        }
//    }
    // time stepping
    for ( int iter = 0; iter < total_iter_; iter++ ) {

        // RK4
        // y(n+1) = y( n ) + 1/6 h ( k1 + 2k2 + 2k3 + k4 )
        // t(n+1) = t( n ) + h

        k1re_->zero();
        k1im_->zero();
        RK4(k1re_,k1im_,k1re_,k1im_,iter,0.0);
        RK4(k2re_,k2im_,k1re_,k1im_,iter,0.5);
        RK4(k3re_,k3im_,k2re_,k2im_,iter,0.5);
        RK4(k4re_,k4im_,k3re_,k3im_,iter,1.0);

        // y(n+1) = y( n ) + 1/6 h ( k1 + 2k2 + 2k3 + k4 )
        for (int i = 0; i < nso_; i++) {
            for (int j = 0; j < nso_; j++) {
                Dre_->pointer()[i][j] += 1.0 / 6.0 * time_step_ * ( k1re_->pointer()[i][j] 
                                                            + 2.0 * k2re_->pointer()[i][j] 
                                                            + 2.0 * k3re_->pointer()[i][j] 
                                                                  + k4re_->pointer()[i][j] );
                Dim_->pointer()[i][j] += 1.0 / 6.0 * time_step_ * ( k1im_->pointer()[i][j] 
                                                            + 2.0 * k2im_->pointer()[i][j] 
                                                            + 2.0 * k3im_->pointer()[i][j] 
                                                                  + k4im_->pointer()[i][j] );
            }
        }

        double dipole_x = C_DDOT(nso_*nso_,&(Dre_->pointer()[0][0]),1,&(mux_->pointer()[0][0]),1);
        double dipole_y = C_DDOT(nso_*nso_,&(Dre_->pointer()[0][0]),1,&(muy_->pointer()[0][0]),1);
        double dipole_z = C_DDOT(nso_*nso_,&(Dre_->pointer()[0][0]),1,&(muz_->pointer()[0][0]),1);

        // TODO: store total signal or only post pulse?
        corr_func_x[iter][0] = dipole_x;
        corr_func_y[iter][0] = dipole_y;
        corr_func_z[iter][0] = dipole_z;

        corr_func_x[iter][1] = 0.0;
        corr_func_y[iter][1] = 0.0;
        corr_func_z[iter][1] = 0.0;

        double energy = CurrentEnergy(Dre_,Dim_);

        outfile->Printf("@TIME %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le\n",iter*time_step_,
            dipole_x,dipole_y,dipole_z,
            ext_field_ * polarization_[0],ext_field_ * polarization_[1],ext_field_ * polarization_[2],energy);
    }

    // FFTW
    FFTW();

    // spectrum
    Spectrum();

    return 0.0;
}

void TD2RDM::RK4(std::shared_ptr<Matrix> koutre, std::shared_ptr<Matrix>koutim,
               std::shared_ptr<Matrix> kinre, std::shared_ptr<Matrix>kinim, int iter, double step) {

    for (int i = 0; i < nso_; i++) {
        for (int j = 0; j < nso_; j++) {
            tempre_->pointer()[i][j] = Dre_->pointer()[i][j] + kinre->pointer()[i][j] * step*time_step_;
            tempim_->pointer()[i][j] = Dim_->pointer()[i][j] + kinim->pointer()[i][j] * step*time_step_;
        }
    }
    BuildFock(tempre_,tempim_,iter*time_step_ + step*time_step_);

    C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(Fim_->pointer()[0][0]),nso_,&(tempre_->pointer()[0][0]),nso_,0.0,&(koutre->pointer()[0][0]),nso_);
    C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(Fre_->pointer()[0][0]),nso_,&(tempim_->pointer()[0][0]),nso_,1.0,&(koutre->pointer()[0][0]),nso_);
    C_DGEMM('n','n',nso_,nso_,nso_,-1.0,&(tempre_->pointer()[0][0]),nso_,&(Fim_->pointer()[0][0]),nso_,1.0,&(koutre->pointer()[0][0]),nso_);
    C_DGEMM('n','n',nso_,nso_,nso_,-1.0,&(tempim_->pointer()[0][0]),nso_,&(Fre_->pointer()[0][0]),nso_,1.0,&(koutre->pointer()[0][0]),nso_);

    C_DGEMM('n','n',nso_,nso_,nso_,-1.0,&(Fre_->pointer()[0][0]),nso_,&(tempre_->pointer()[0][0]),nso_,0.0,&(koutim->pointer()[0][0]),nso_);
    C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(Fim_->pointer()[0][0]),nso_,&(tempim_->pointer()[0][0]),nso_,1.0,&(koutim->pointer()[0][0]),nso_);
    C_DGEMM('n','n',nso_,nso_,nso_,1.0,&(tempre_->pointer()[0][0]),nso_,&(Fre_->pointer()[0][0]),nso_,1.0,&(koutim->pointer()[0][0]),nso_);
    C_DGEMM('n','n',nso_,nso_,nso_,-1.0,&(tempim_->pointer()[0][0]),nso_,&(Fim_->pointer()[0][0]),nso_,1.0,&(koutim->pointer()[0][0]),nso_);

    // 2-cumulant contribution
    //printf("compute 2-cumulant \n"); 
//    tempre_->print();
//    tempim_->print();
//    std::shared_ptr<Vector> Eval (new Vector(nso_));
//    std::shared_ptr<Matrix> Evec (new Matrix(nso_,nso_));
//    tempre_->diagonalize(Evec, Eval);
//    Eval->print();
//    Evec->print();
//    koutre->print();
//    koutim->print();
//    kre_->zero();
//    kim_->zero();
    Compute2Cumulant(Fre_,Fim_,&tempre_->pointer()[0][0],&tempim_->pointer()[0][0]);
//    koutre->add(kre_);
//    koutim->add(kim_);
    koutre->print();
    koutim->print();
    //printf("2-cumulant computed \n"); 
}

void TD2RDM::Compute2Cumulant(std::shared_ptr<Matrix> Fre,std::shared_ptr<Matrix> Fim,double * Dre0,double * Dim0) {
    
    std::shared_ptr<Vector> Eval (new Vector(nso_));
    std::shared_ptr<Matrix> Evec (new Matrix(nso_,nso_));
    std::shared_ptr<Matrix> treold (new Matrix(nso_*nso_,nso_*nso_));
    std::shared_ptr<Matrix> timold (new Matrix(nso_*nso_,nso_*nso_));
    std::shared_ptr<Matrix> tre (new Matrix(nso_*nso_,nso_*nso_));
    std::shared_ptr<Matrix> tim (new Matrix(nso_*nso_,nso_*nso_));

    double * Dre = (double*)malloc(nso_*nso_*sizeof(double));
    double * Dim = (double*)malloc(nso_*nso_*sizeof(double));
    double * tempre = (double*)malloc(nso_*nso_*sizeof(double));
    double * tempim = (double*)malloc(nso_*nso_*sizeof(double));
    double * Ftempre = (double*)malloc(nso_*nso_*sizeof(double));
    double * Ftempim = (double*)malloc(nso_*nso_*sizeof(double));
    double * Qnore = (double*)malloc(nQ_*nso_*nso_*sizeof(double));
    double * Qnoim = (double*)malloc(nQ_*nso_*nso_*sizeof(double));
    double * Eval_p = Eval->pointer();

    for (int p = 0; p < nso_*nso_; p++){
        Dre[p] = Dre0[p];
        Dim[p] = Dim0[p];
    }
   
    // find occupation numbers
    //printf("Diagonalizing 1RDM \n"); 
//    for (int p = 0; p < nso_; p++){
//        for (int q = 0; q < nso_; q++){
//            outfile->Printf(" %20.12lf ",Dre[p*nso_+q]);
//        }
//        outfile->Printf("\n");
//    }
//    outfile->Printf("\n");
//    for (int p = 0; p < nso_; p++){
//        for (int q = 0; q < nso_; q++){
//            outfile->Printf(" %20.12lf ",Dim[p*nso_+q]);
//        }
//        outfile->Printf("\n");
//    }
//    outfile->Printf("\n");

    DiagonalizeHermitianMatrix(nso_,Dre,Dim,Eval_p);
    
//    for (int p = 0; p < nso_; p++){
//        for (int q = 0; q < nso_; q++){
//            outfile->Printf(" %20.12lf ",Dre[p*nso_+q]);
//        }
//        outfile->Printf("\n");
//    }
//    outfile->Printf("\n");
//    for (int p = 0; p < nso_; p++){
//        for (int q = 0; q < nso_; q++){
//            outfile->Printf(" %20.12lf ",Dim[p*nso_+q]);
//        }
//        outfile->Printf("\n");
//    }
//    outfile->Printf("\n");
//    for (int p = 0; p < nso_; p++){
//        outfile->Printf(" %20.12lf ",Eval_p[p]);
//    }
//    outfile->Printf("\n");

    // transform Fock matrix
    //printf("Transforming Fock matrix to NO basis \n"); 
    for (int p = 0; p < nso_; p++){
        for (int q = 0; q < nso_; q++){
            tempre[p*nso_+q] = 0.0; 
            tempim[p*nso_+q] = 0.0; 
            for (int r = 0; r < nso_; r++){
                tempre[p*nso_+q] += Fre->pointer()[p][r]*Dre[r*nso_+q] - Fim->pointer()[p][r]*Dim[r*nso_+q]; 
                tempim[p*nso_+q] += Fre->pointer()[p][r]*Dim[r*nso_+q] + Fim->pointer()[p][r]*Dre[r*nso_+q]; 
            }
        }
    }
    for (int p = 0; p < nso_; p++){
        for (int q = 0; q < nso_; q++){
            Ftempre[p*nso_+q] = 0.0; 
            Ftempim[p*nso_+q] = 0.0; 
            for (int r = 0; r < nso_; r++){
                Ftempre[p*nso_+q] += Dre[r*nso_+p]*tempre[r*nso_+q]; 
                Ftempim[p*nso_+q] += Dre[r*nso_+p]*tempim[r*nso_+q]; 
            }
        }
    }
    for (int p = 0; p < nso_; p++){
        for (int q = 0; q < nso_; q++){
            tempre[p*nso_+q] = 0.0; 
            tempim[p*nso_+q] = 0.0; 
            for (int r = 0; r < nso_; r++){
                tempim[p*nso_+q] += Fre->pointer()[p][r]*Dim[r*nso_+q] + Fim->pointer()[p][r]*Dre[r*nso_+q]; 
                tempre[p*nso_+q] += Fre->pointer()[p][r]*Dre[r*nso_+q] - Fim->pointer()[p][r]*Dim[r*nso_+q]; 
            }
        }
    }
    for (int p = 0; p < nso_; p++){
        for (int q = 0; q < nso_; q++){
            for (int r = 0; r < nso_; r++){
                Ftempre[p*nso_+q] += Dim[r*nso_+p]*tempim[r*nso_+q]; 
                Ftempim[p*nso_+q] -= Dim[r*nso_+p]*tempre[r*nso_+q]; 
            }
        outfile->Printf("%20.12lf ",Ftempre[p*nso_+q]);
        }
    outfile->Printf("\n");
    }
    // make sure no negative zeros are present in the Fock matrix
    for (int p = 0; p< nso_*nso_; p++){
        if (Ftempre[p] > -1.0e-12 && Ftempre[p] < 1.0e-12){
            Ftempre[p] = 0.0;
        }
        if (Ftempim[p] > -1.0e-12 && Ftempim[p] < 1.0e-12){
            Ftempim[p] = 0.0;
        }
    }
//    outfile->Printf("\n");
    free(tempre);
    free(tempim);
    // transform three-index integrals
    //printf("Transforming 3-index integrals to NO basis \n"); 
    tempre = (double*)malloc(nQ_*nso_*nso_*sizeof(double));
    tempim = (double*)malloc(nQ_*nso_*nso_*sizeof(double));
    for (int Q = 0; Q < nQ_; Q++){
        for (int p = 0; p < nso_; p++){
            for (int q = 0; q < nso_; q++){
                int pq = Q*nso_*nso_+p*nso_+q;
                tempre[pq] = 0.0; 
                tempim[pq] = 0.0; 
                for (int r = 0; r < nso_; r++){
                    int pr = p*nso_+r;
                    int rq = r*nso_+q;
                    tempre[pq] += Qp_[Q][pr]*Dre[rq] - Qp_[Q][pr]*Dim[rq]; 
                    tempim[pq] += Qp_[Q][pr]*Dim[rq] + Qp_[Q][pr]*Dre[rq]; 
                }
            }
        }
    }
    for (int Q = 0; Q < nQ_; Q++){
        for (int p = 0; p < nso_; p++){
            for (int q = 0; q < nso_; q++){
                int pq = Q*nso_*nso_+p*nso_+q;
                Qnore[pq] = 0.0; 
                Qnoim[pq] = 0.0; 
                for (int r = 0; r < nso_; r++){
                    int rp = r*nso_+p;
                    int rq = Q*nso_*nso_+r*nso_+q;
                    Qnore[pq] += Dre[rp]*tempre[rq]; 
                    Qnoim[pq] += Dre[rp]*tempim[rq]; 
                }
            }
        }
    }
    for (int Q = 0; Q < nQ_; Q++){
        for (int p = 0; p < nso_; p++){
            for (int q = 0; q < nso_; q++){
                int pq = Q*nso_*nso_+p*nso_+q;
                tempre[pq] = 0.0; 
                tempim[pq] = 0.0; 
                for (int r = 0; r < nso_; r++){
                    int pr = p*nso_+r;
                    int rq = r*nso_+q;
                    tempim[pq] += Qp_[Q][pr]*Dim[rq] + Qp_[Q][pr]*Dre[rq]; 
                    tempre[pq] += Qp_[Q][pr]*Dre[rq] - Qp_[Q][pr]*Dim[rq]; 
                }
            }
        }
    }
    for (int Q = 0; Q < nQ_; Q++){
        for (int p = 0; p < nso_; p++){
            for (int q = 0; q < nso_; q++){
                int pq = Q*nso_*nso_+p*nso_+q;
                for (int r = 0; r < nso_; r++){
                    int rp = r*nso_+p;
                    int rq = Q*nso_*nso_+r*nso_+q;
                    Qnore[pq] += Dim[rp]*tempim[rq]; 
                    Qnoim[pq] -= Dim[rp]*tempre[rq]; 
                }
            outfile->Printf("%20.12lf ",Qnore[pq]);
            }
        }
    outfile->Printf("\n");
    }
    outfile->Printf("\n");
    free(tempre);
    free(tempim);
    for (int Q = 0; Q< nQ_; Q++){
        for (int p = 0; p< nso_*nso_; p++){
            int Qp = Q*nso_*nso_+p;
            if (Qnore[Qp] > -1.0e-12 && Qnore[Qp] < 1.0e-12){
                Qnore[Qp] = 0.0;
            }
            if (Qnoim[Qp] > -1.0e-12 && Qnoim[Qp] < 1.0e-12){
                Qnoim[Qp] = 0.0;
            }
        }
    }

//    for (int p = 0; p < nso_; p++){
//        for (int q = 0; q < nso_; q++){
//            Evec->pointer()[p][q] = Evec_p[p*nso_+q];
//        }
//    }
//    free(Eval_p);
//    Eval->print();
//    for (int p = 0; p < nso_; p++){
//        for (int q = 0; q < nso_; q++){
//            outfile->Printf("%20.12lf ",Dre[p*nso_+q]);
//        }
//        outfile->Printf("\n");
//    }
//    outfile->Printf("\n");
//    for (int p = 0; p < nso_; p++){
//        for (int q = 0; q < nso_; q++){
//            outfile->Printf("%20.12lf ",Dim[p*nso_+q]);
//        }
//        outfile->Printf("\n");
//    }
//    return;
//    exit(0);
//    Dre->diagonalize(Evec,Eval,descending);
//    Fre->transform(Evec);
//    Fim->transform(Evec);

    // iterative procedure to find T
    for (long int p = 0; p < nso_; p++){
        double Ep = Eval_p[p]; 
        for (long int q = 0; q < nso_; q++){
            double Eq = Eval_p[q]; 
            for (long int r = 0; r < nso_; r++){
                double Er = Eval_p[r]; 
                for (long int s = 0; s < nso_; s++){
                    double Es = Eval_p[s]; 
                    for (long int Q = 0; Q < nQ_; Q++){
                        int pr = Q*nso_*nso_+r*nso_+p;
                        int qs = Q*nso_*nso_+s*nso_+q;
                        tre->pointer()[p*nso_+q][r*nso_+s] += Qnore[pr]*Qnore[qs] * (1.0-Ep)*(1.0-Eq)*(Er)*(Es); 
                        tre->pointer()[p*nso_+q][r*nso_+s] += Qnoim[pr]*Qnoim[qs] * (1.0-Ep)*(1.0-Eq)*(Er)*(Es); 

                        tim->pointer()[p*nso_+q][r*nso_+s] += Qnore[pr]*Qnoim[qs] * (1.0-Ep)*(1.0-Eq)*(Er)*(Es); 
                        tim->pointer()[p*nso_+q][r*nso_+s] -= Qnoim[pr]*Qnore[qs] * (1.0-Ep)*(1.0-Eq)*(Er)*(Es); 
                    }
                }
            }
        }
    }
//    tre->print();
//    tim->print();

    //printf("begining iterative procedure \n"); 
    for (int iter = 0; iter < 10000; iter++){
        treold->copy(tre);
        timold->copy(tim);
        tre->zero();
        tim->zero();
        double deltaT = 0.0;
        for (long int p = 0; p < nso_; p++){
            double Ep = Eval_p[p]; 
            for (long int q = 0; q < nso_; q++){
                double Eq = Eval_p[q]; 
                for (long int r = 0; r < nso_; r++){
                    double Er = Eval_p[r]; 
                    for (long int s = 0; s < nso_; s++){
                        double Es = Eval_p[s]; 
                        for (long int Q = 0; Q < nQ_; Q++){
                            int pr = Q*nso_*nso_+r*nso_+p;
                            int qs = Q*nso_*nso_+s*nso_+q;
                            tre->pointer()[p*nso_+q][r*nso_+s] += Qnore[pr]*Qnore[qs] * (1.0-Ep)*(1.0-Eq)*(Er)*(Es); 
                            tre->pointer()[p*nso_+q][r*nso_+s] += Qnoim[pr]*Qnoim[qs] * (1.0-Ep)*(1.0-Eq)*(Er)*(Es); 

                            tim->pointer()[p*nso_+q][r*nso_+s] += Qnore[pr]*Qnoim[qs] * (1.0-Ep)*(1.0-Eq)*(Er)*(Es); 
                            tim->pointer()[p*nso_+q][r*nso_+s] -= Qnoim[pr]*Qnore[qs] * (1.0-Ep)*(1.0-Eq)*(Er)*(Es); 
                        }
                        for (long int t = 0; t < nso_; t++){
                            if (t != p){
                                int tp = t*nso_+p;
                                tre->pointer()[p*nso_+q][r*nso_+s] -= Ftempre[tp] * treold->pointer()[t*nso_+q][r*nso_+s];
                                tre->pointer()[p*nso_+q][r*nso_+s] += Ftempim[tp] * timold->pointer()[t*nso_+q][r*nso_+s];

                                tim->pointer()[p*nso_+q][r*nso_+s] -= Ftempre[tp] * timold->pointer()[t*nso_+q][r*nso_+s];
                                tim->pointer()[p*nso_+q][r*nso_+s] -= Ftempim[tp] * treold->pointer()[t*nso_+q][r*nso_+s];
                            }
                            if (t != q){
                                int tq = t*nso_+q;
                                tre->pointer()[p*nso_+q][r*nso_+s] -= Ftempre[tq] * treold->pointer()[p*nso_+t][r*nso_+s];
                                tre->pointer()[p*nso_+q][r*nso_+s] += Ftempim[tq] * timold->pointer()[p*nso_+t][r*nso_+s];

                                tim->pointer()[p*nso_+q][r*nso_+s] -= Ftempre[tq] * timold->pointer()[p*nso_+t][r*nso_+s];
                                tim->pointer()[p*nso_+q][r*nso_+s] -= Ftempim[tq] * treold->pointer()[p*nso_+t][r*nso_+s];
                            }
                            if (t != r){
                                int rt = r*nso_+t;
                                tre->pointer()[p*nso_+q][r*nso_+s] += Ftempre[rt] * treold->pointer()[p*nso_+q][t*nso_+s];
                                tre->pointer()[p*nso_+q][r*nso_+s] -= Ftempim[rt] * timold->pointer()[p*nso_+q][t*nso_+s];

                                tim->pointer()[p*nso_+q][r*nso_+s] += Ftempre[rt] * timold->pointer()[p*nso_+q][t*nso_+s];
                                tim->pointer()[p*nso_+q][r*nso_+s] += Ftempim[rt] * treold->pointer()[p*nso_+q][t*nso_+s];
                            }
                            if (t != s){
                                int st = s*nso_+t;
                                tre->pointer()[p*nso_+q][r*nso_+s] += Ftempre[st] * treold->pointer()[p*nso_+q][r*nso_+t];
                                tre->pointer()[p*nso_+q][r*nso_+s] -= Ftempim[st] * timold->pointer()[p*nso_+q][r*nso_+t];

                                tim->pointer()[p*nso_+q][r*nso_+s] += Ftempre[st] * timold->pointer()[p*nso_+q][r*nso_+t];
                                tim->pointer()[p*nso_+q][r*nso_+s] += Ftempim[st] * treold->pointer()[p*nso_+q][r*nso_+t];
                            }
                        }
                        int pp = p*nso_+p;
                        int qq = q*nso_+q;
                        int rr = r*nso_+r;
                        int ss = s*nso_+s;
                        double denom = (Ftempre[pp] + Ftempre[qq] - Ftempre[rr] - Ftempre[ss])
                                     * (Ftempre[pp] + Ftempre[qq] - Ftempre[rr] - Ftempre[ss])
                                     + (Ftempim[pp] + Ftempim[qq] - Ftempim[rr] - Ftempim[ss])
                                     * (Ftempim[pp] + Ftempim[qq] - Ftempim[rr] - Ftempim[ss]); 
                        double tempTre = 0.0;
                        double tempTim = 0.0;
                        if (denom > 1e-12) {
                           tempTre = (tre->pointer()[p*nso_+q][r*nso_+s]
                                   * (Ftempre[pp] + Ftempre[qq] - Ftempre[rr] - Ftempre[ss])
                                   + tim->pointer()[p*nso_+q][r*nso_+s] 
                                   * (Ftempim[pp] + Ftempim[qq] - Ftempim[rr] - Ftempim[ss]))/denom;
                           tempTim = (tim->pointer()[p*nso_+q][r*nso_+s] 
                                   * (Ftempre[pp] + Ftempre[qq] - Ftempre[rr] - Ftempre[ss])
                                   - tre->pointer()[p*nso_+q][r*nso_+s]
                                   * (Ftempim[pp] + Ftempim[qq] - Ftempim[rr] - Ftempim[ss]))/denom;
                        }
                        deltaT += pow((tempTre - treold->pointer()[p*nso_+q][r*nso_+s]),2); 
                        deltaT += pow((tempTim - timold->pointer()[p*nso_+q][r*nso_+s]),2);
                        if (tempTre > 1.0e-12 || tempTre < -1.0e-12){
                            tre->pointer()[p*nso_+q][r*nso_+s] = tempTre;
                        } else tre->pointer()[p*nso_+q][r*nso_+s] = 0.0;
                        if (tempTim > 1.0e-12 || tempTim < -1.0e-12){
                            tim->pointer()[p*nso_+q][r*nso_+s] = tempTim;
                        } else tim->pointer()[p*nso_+q][r*nso_+s] = 0.0;
                    }
                }
            }
        }
//        tre->print();
//        tim->print();
//        exit(0);
        outfile->Printf(" |Delta T| =  %20.12lf  %20.12lf \n",sqrt(deltaT), deltaT);
        if (sqrt(deltaT) < 1e-12){break;}
    }
    //printf("exiting iterative procedure \n"); 
//    free(Ftempre);
//    free(Ftempim);
    free(Qnore);
    free(Qnoim);
    //printf("computing 2-cumulant \n");
    // building 2-cumulant  
    for (long int p = 0; p < nso_; p++){
        for (long int q = 0; q < nso_; q++){
            for (long int r = 0; r < nso_; r++){
                for (long int s = 0; s < nso_; s++){
                    C2re_->pointer()[p*nso_+q][r*nso_+s] = tre->pointer()[p*nso_+q][r*nso_+s] + tre->pointer()[r*nso_+s][p*nso_+q]; 
                    C2im_->pointer()[p*nso_+q][r*nso_+s] = tim->pointer()[p*nso_+q][r*nso_+s] - tim->pointer()[r*nso_+s][p*nso_+q]; 
                    if (p == r && q == s){
//                        C2re_->pointer()[p*nso_+q][r*nso_+s] -= Eval_p[p]*Eval_p[q]; 
                    } 
                }
            }
        }
    }  
    C2re_->print();
    C2im_->print();
  
    for (long int p = 0; p < nso_; p++){
        for (long int q = 0; q < nso_; q++){
            for (long int r = 0; r < nso_; r++){
                for (long int s = 0; s < nso_; s++){
                    tre->pointer()[p*nso_+q][r*nso_+s] = 1.0*C2re_->pointer()[p*nso_+q][r*nso_+s] - 0.0*C2re_->pointer()[p*nso_+q][s*nso_+r];
                    tim->pointer()[p*nso_+q][r*nso_+s] = 1.0*C2im_->pointer()[p*nso_+q][r*nso_+s] - 0.0*C2im_->pointer()[p*nso_+q][s*nso_+r];
                }
            }
        }
    }  
//    tre->copy(C2re_);
//    tim->copy(C2im_);
    tre->print();
    tim->print();
    // I[pj] = Delta[pqjs]*Qqs
    tempre = (double*)malloc(nso_*nso_*sizeof(double));
    tempim = (double*)malloc(nso_*nso_*sizeof(double));
    for (long int p = 0; p < nso_; p++){
        for (long int j = 0; j < nso_; j++){
            tempre[p*nso_+j] = 0.0; 
            tempim[p*nso_+j] = 0.0;
            for (long int q = 0; q < nso_; q++){
                double Eq = Eval_p[q];
                for (long int s = 0; s < nso_; s++){
                    double Es = Eval_p[s];
                    for (long int Q = 0; Q < nQ_; Q++){
                        long int pqjs = p*nso_*nso_*nso_+q*nso_*nso_+j*nso_+s;
                        tempre[p*nso_+j] -= tre->pointer()[p*nso_+q][j*nso_+s]*Qnore[Q*nso_*nso_+q*nso_+s] * (1.0-Eq)*(Es);
                        tempre[p*nso_+j] += tim->pointer()[p*nso_+q][j*nso_+s]*Qnoim[Q*nso_*nso_+q*nso_+s] * (1.0-Eq)*(Es);

                        tempim[p*nso_+j] -= tre->pointer()[p*nso_+q][j*nso_+s]*Qnoim[Q*nso_*nso_+q*nso_+s] * (1.0-Eq)*(Es);
                        tempim[p*nso_+j] -= tim->pointer()[p*nso_+q][j*nso_+s]*Qnore[Q*nso_*nso_+q*nso_+s] * (1.0-Eq)*(Es);
                    }
                }
            }
        }
    }  
    for (long int i = 0; i < nso_; i++){
        double Ei = Eval_p[i];
        for (long int j = 0; j < nso_; j++){
            Ftempre[i*nso_+j] = 0.0; 
            Ftempim[i*nso_+j] = 0.0; 
            for (long int p = 0; p < nso_; p++){
                double Ep = Eval_p[p];
                for (long int Q = 0; Q < nQ_; Q++){
                    Ftempre[i*nso_+j] += 0.5*tempre[p*nso_+j]*Qnore[Q*nso_*nso_+p*nso_+i] * (1.0-Ep)*(Ei);
                    Ftempre[i*nso_+j] += 0.5*tempim[p*nso_+j]*Qnoim[Q*nso_*nso_+p*nso_+i] * (1.0-Ep)*(Ei);

//                    Ftempim[i*nso_+j] += 0.5*tempre[p*nso_+j]*Qnoim[Q*nso_*nso_+p*nso_+i];
//                    Ftempim[i*nso_+j] += 0.5*tempim[p*nso_+j]*Qnore[Q*nso_*nso_+p*nso_+i];
                }
            }
        }
    }  
    // I[ip] = Delta[ispq]*Qqs
    for (long int i = 0; i < nso_; i++){
        for (long int p = 0; p < nso_; p++){
            tempre[i*nso_+p] = 0.0;
            tempim[i*nso_+p] = 0.0;
            for (long int q = 0; q < nso_; q++){
                double Eq = Eval_p[q];
                for (long int s = 0; s < nso_; s++){
                    double Es = Eval_p[s];
                    for (long int Q = 0; Q < nQ_; Q++){
                        tempre[i*nso_+p] -= tre->pointer()[i*nso_+s][p*nso_+q]*Qnore[Q*nso_*nso_+s*nso_+q] * (1.0-Es)*(Eq);
                        tempre[i*nso_+p] += tim->pointer()[i*nso_+s][p*nso_+q]*Qnoim[Q*nso_*nso_+s*nso_+q] * (1.0-Es)*(Eq);

                        tempim[i*nso_+p] -= tre->pointer()[i*nso_+s][p*nso_+q]*Qnoim[Q*nso_*nso_+s*nso_+q] * (1.0-Es)*(Eq);
                        tempim[i*nso_+p] -= tim->pointer()[i*nso_+s][p*nso_+q]*Qnore[Q*nso_*nso_+s*nso_+q] * (1.0-Es)*(Eq);
                    }
                }
            }
        }
    }
    for (long int i = 0; i < nso_; i++){
        for (long int j = 0; j < nso_; j++){
            double Ej = Eval_p[j];
            for (long int p = 0; p < nso_; p++){
                double Ep = Eval_p[p];
                for (long int Q = 0; Q < nQ_; Q++){
                    Ftempre[i*nso_+j] -= 0.5*tempre[i*nso_+p]*Qnore[Q*nso_*nso_+j*nso_+p] * (1.0-Ej)*(Ep);
                    Ftempre[i*nso_+j] -= 0.5*tempim[i*nso_+p]*Qnoim[Q*nso_*nso_+j*nso_+p] * (1.0-Ej)*(Ep);

//                    Ftempim[i*nso_+j] += 0.5*tempre[i*nso_+p]*Qnoim[Q*nso_*nso_+j*nso_+p];
//                    Ftempim[i*nso_+j] += 0.5*tempim[i*nso_+p]*Qnore[Q*nso_*nso_+j*nso_+p];
                }
            }
            outfile->Printf("%20.12lf (%20.12lf) ",Ftempre[i*nso_+j],Ftempim[i*nso_+j]);
        }
        outfile->Printf("\n");
    }
    outfile->Printf("\n");

    // back transform to MO basis
    kre_->zero();
    kim_->zero();
    for (int p = 0; p < nso_; p++){
        for (int q = 0; q < nso_; q++){
            tempre[p*nso_+q] = 0.0;
            tempim[p*nso_+q] = 0.0;
            for (int r = 0; r < nso_; r++){
                tempre[p*nso_+q] += Ftempre[p*nso_+r]*Dre[r*nso_+q];// + Ftempim[p*nso_+r]*Dim[q*nso_+r];
                tempim[p*nso_+q] += -Ftempre[p*nso_+r]*Dim[r*nso_+q];// + Ftempim[p*nso_+r]*Dre[q*nso_+r];
            }
        }
    }
    for (int p = 0; p < nso_; p++){
        for (int q = 0; q < nso_; q++){
            for (int r = 0; r < nso_; r++){
                kim_->pointer()[p][q] += Dre[p*nso_+r]*tempre[r*nso_+q] - Dim[p*nso_+r]*tempim[r*nso_+q];
                kre_->pointer()[p][q] += Dre[p*nso_+r]*tempim[r*nso_+q] + Dim[p*nso_+r]*tempre[r*nso_+q];
            }
        }
    }
//    for (int p = 0; p < nso_; p++){
//        for (int q = 0; q < nso_; q++){
//            tempre[p*nso_+q] = 0.0;
//            tempim[p*nso_+q] = 0.0;
//            for (int r = 0; r < nso_; r++){
//                tempim[p*nso_+q] += Fre->pointer()[p][r]*Dim[q*nso_+r] - Fim->pointer()[p][r]*Dre[q*nso_+r];
//                tempre[p*nso_+q] += Fre->pointer()[p][r]*Dre[q*nso_+r] + Fim->pointer()[p][r]*Dim[q*nso_+r];
//            }
//        }
//    }
//    for (int p = 0; p < nso_; p++){
//        for (int q = 0; q < nso_; q++){
//            for (int r = 0; r < nso_; r++){
//                kim_->pointer()[p][q] -= Dim[p*nso_+r]*tempim[r*nso_+q];
//                kre_->pointer()[p][q] -= Dim[p*nso_+r]*tempre[r*nso_+q];
//            }
//        }
//    }

    for (int p = 0; p < nso_; p++){
        for (int q = 0; q < nso_; q++){
            if (kre_->pointer()[p][q] < 1.0e-12 && kre_->pointer()[p][q] > -1.0e-12){
                kre_->pointer()[p][q] = 0.0;
            }
            if (kim_->pointer()[p][q] < 1.0e-12 && kim_->pointer()[p][q] > -1.0e-12){
                kim_->pointer()[p][q] = 0.0;
            }
        }
    }

    kre_->print();
    kim_->print();
//    exit(0);
    
    // computing contribution to the EOM

//    C2re_->print();
//    C2im_->print();

    //free(Eval_p);
}

double TD2RDM::CurrentEnergy(std::shared_ptr<Matrix> Dre,std::shared_ptr<Matrix> Dim) {

    BuildJK(Dre,Dim);

    Jre_->scale(2.0);
    Fre_->copy(Hcore_);
    Fre_->add(Hcore_);
    Fre_->add(Jre_);
    Fre_->subtract(Kre_);

    double energy = C_DDOT(nso_*nso_,Dre->pointer()[0],1,Fre_->pointer()[0],1);

    return energy;
}
void TD2RDM::BuildFock(std::shared_ptr<Matrix> Dre,std::shared_ptr<Matrix> Dim, double curtime) {

    BuildJK(Dre,Dim);

    Jre_->scale(2.0);
    Fre_->copy(Hcore_);
    Fre_->add(Jre_);

    ExtField(curtime);

    Jim_->scale(2.0);
    Fim_->zero();
    Fim_->add(Jim_);

    //if ( is_dft_ ) {
    //    BuildV(Dre,Dim);

    //    Fre_->add(Vre_);
    //    //Fim_->add(Vim_);

    //    double alpha = functional_->x_alpha();
    //    double beta = 1.0 - alpha;

    //    if (alpha != 0.0) {

    //        Kre_->scale(alpha);
    //        Fre_->subtract(Kre_);
    //        Kre_->scale(1.0/alpha);
    //
    //        Kim_->scale(alpha);
    //        Fim_->subtract(Kim_);
    //        Kim_->scale(1.0/alpha);
    //
    //    } else {
    //
    //        //Kre_->zero();
    //        //Kim_->zero();
    //
    //    }
    //
    //    if (functional_->is_x_lrc()) {
    //        wKre_->scale(beta);
    //        Fre_->subtract(wKre_);
    //        wKre_->scale(1.0/beta);
    //
    //        wKim_->scale(beta);
    //        Fim_->subtract(wKim_);
    //        wKim_->scale(1.0/beta);
    //    } else {
    //        //wKre_->zero();
    //        //wKim_->zero();
    //    }

    //}else {
        Fre_->subtract(Kre_);
        Fim_->subtract(Kim_);
    //}

    //Fre_->print();
    //Fim_->print();
    //exit(0);

}

void TD2RDM::ExtField(double curtime) {
    // add external field

    ext_field_ = 0.0;

    if (pulse_shape_ == 0 ) {

        // from prl:
        if ( curtime < laser_time_ ) {
            ext_field_ = sin(M_PI*curtime/(laser_time_));
            ext_field_ *= ext_field_*laser_amp_*sin(laser_freq_*curtime);
        }

    } else if ( pulse_shape_ == 1 ) {

        // from 2007 schlegel paper (jcp 126, 244110 (2007))
        if (curtime <= 2.0 * M_PI / laser_freq_)      ext_field_ = laser_freq_ * curtime / (2.0 * M_PI) * laser_amp_;
        else if (curtime <= 4.0 * M_PI / laser_freq_) ext_field_ = laser_amp_;
        else if (curtime <= 6.0 * M_PI / laser_freq_) ext_field_ = (3.0 - laser_freq_ * curtime / (2.0 * M_PI) ) * laser_amp_;
        ext_field_ *= sin(laser_freq_*curtime);

    } else if ( pulse_shape_ == 2 ) {

        // pi pulse from licn paper
        double sigma = laser_time_*0.5;
        if ( curtime < laser_time_ ) {
            ext_field_ = cos(M_PI*curtime/(2.0*sigma) + 0.5*M_PI);
            ext_field_ *= M_PI/(sigma*transition_dpm_) * ext_field_ * laser_amp_ * cos(laser_freq_*curtime);
        }

    } else if ( pulse_shape_ == 3 ) {

        // continuous wave for rabi flopping
        ext_field_ = laser_amp_*sin(laser_freq_*curtime);
    } else if ( pulse_shape_ == 4 ) {

        // Gaussian pulse
        ext_field_ = laser_amp_ * exp(-((curtime-1.5*laser_time_)*(curtime-1.5*laser_time_))
                  / (0.3606738*laser_time_*laser_time_)) * sin(laser_freq_*curtime);

    }

    for (int i = 0; i < nso_; i++) {
        for (int j = 0; j < nso_; j++) {
            Fre_->pointer()[i][j] -= ext_field_ * mux_->pointer()[i][j] * polarization_[0];
            Fre_->pointer()[i][j] -= ext_field_ * muy_->pointer()[i][j] * polarization_[1];
            Fre_->pointer()[i][j] -= ext_field_ * muz_->pointer()[i][j] * polarization_[2];
        }
    }

}

void TD2RDM::BuildV(std::shared_ptr<Matrix> Dre,std::shared_ptr<Matrix> Dim) {

    // TODO: no such thing as C_left/C_right ... need to decompose density
    // D  = R^TR
    // C' = CR
    // Alternatively, the potential could be evaluated in the natural orbital basis ...
    // According to Lopata and Govind JCTC 7, 1334-1355 (2011), the XC potential
    // only depends on the real part of the density.  So, I'll just only need to
    // factorize (or diagonalize) the real part.

    // diagonalize real part of the density:
    std::shared_ptr<Matrix> eigvec (new Matrix(Dre));
    std::shared_ptr<Vector> eigval (new Vector(nso_));
    Dre->diagonalize(eigvec,eigval,descending);

    std::vector<SharedMatrix>& C  = potential_->C();

    // use real and imaginary C matrices

    std::shared_ptr<Matrix> Cre (new Matrix(Ca_) );
    //std::shared_ptr<Matrix> Cim (new Matrix(Ca_) );

    Cre->zero();
    //Cim->zero();

    C.clear();

    // TODO: make sure im doing this right ...
    for (int i = 0; i < nso_; i++) {
        for (int j = 0; j < nso_; j++) {
            double dum = 0.0;
            for (int k = 0; k < nso_; k++) {
                dum += Ca_->pointer()[i][k] * eigvec->pointer()[k][j];
            }
            Cre->pointer()[i][j] = dum * sqrt(fabs(eigval->pointer()[j]));
        }
    }

    //Cre->gemm('t','n',1.0,Dre,Ca_,0.0);
    //Cim->gemm('t','n',1.0,Dim,Ca_,0.0);

    //Cre->transpose_this();
    //Cim->transpose_this();

    C.push_back(Cre);
    //C.push_back(reference_wavefunction_->Ca_subset("SO", "OCC"));
    //C.push_back(Cim);

    // run the potential object
    potential_->compute();

    Vre_ = potential_->V()[0];
    //Vim_ = potential_->V()[1];

    // transform V matrix to the MO basis:

    Vre_->transform(Ca_);
    //Vim_->transform(Ca_);

}

void TD2RDM::BuildJK(std::shared_ptr<Matrix> Dre,std::shared_ptr<Matrix> Dim) {

    std::vector<SharedMatrix>& C_left  = jk_->C_left();
    std::vector<SharedMatrix>& C_right = jk_->C_right();

    // use real and imaginary C matrices

    std::shared_ptr<Matrix> Cre (new Matrix(Ca_) );
    std::shared_ptr<Matrix> Cim (new Matrix(Ca_) );

    Cre->zero();
    Cim->zero();

    C_left.clear();
    C_right.clear();

    Cre->gemm('t','n',1.0,Dre,Ca_,0.0);
    Cim->gemm('t','n',1.0,Dim,Ca_,0.0);

    Cre->transpose_this();
    Cim->transpose_this();

    C_left.push_back(Cre);
    C_left.push_back(Cim);

    C_right.push_back(Ca_);
    C_right.push_back(Ca_);

    // Let jk compute for the given C_left/C_right

    jk_->compute();

    Jre_ = jk_->J()[0];
    Jim_ = jk_->J()[1];


    // transform J and K matrices to the MO basis:

    Jre_->transform(Ca_);
    Jim_->transform(Ca_);

    //if ( !is_dft_ ) {
        Kre_ = jk_->K()[0];
        Kim_ = jk_->K()[1];
        Kre_->transform(Ca_);
        Kim_->transform(Ca_);
    //}else if ( is_dft_ && functional_->is_x_hybrid() ) {
    //    Kre_ = jk_->K()[0];
    //    Kim_ = jk_->K()[1];
    //    Kre_->transform(Ca_);
    //    Kim_->transform(Ca_);
    //}

}

void TD2RDM::FFTW(){
    fftw_plan p_x;
    p_x = fftw_plan_dft_1d((int)(total_time_ / time_step_),corr_func_x,corr_func_x,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p_x);
    fftw_destroy_plan(p_x);

    fftw_plan p_y;
    p_y = fftw_plan_dft_1d((int)(total_time_ / time_step_),corr_func_y,corr_func_y,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p_y);
    fftw_destroy_plan(p_y);

    fftw_plan p_z;
    p_z = fftw_plan_dft_1d((int)(total_time_ / time_step_),corr_func_z,corr_func_z,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p_z);
    fftw_destroy_plan(p_z);
}
void TD2RDM::Spectrum(){
     outfile->Printf("\n");
     outfile->Printf("        ***********************************************************\n");
     outfile->Printf("        *                                                         *\n");
     outfile->Printf("        *                    Emission spectrum                    *\n");
     outfile->Printf("        *         as computed by the Fourier transform of         *\n");
     outfile->Printf("        *                   dipole acceleration, electric fields  *\n");
     outfile->Printf("        *                                                         *\n");
     outfile->Printf("        *     I(w) = |FourierTransform ( d^2 D(t) / dt^2 )|^2     *\n");
     outfile->Printf("        ***********************************************************\n");
     outfile->Printf("\n");
     outfile->Printf("                           w(eV)");
     outfile->Printf("                 Ix(w)");
     outfile->Printf("                 Iy(w)");
     outfile->Printf("                 Iz(w)");
     outfile->Printf("\n\n");

     // components of external and induced fields
     double max_freq = 200.0;
     double twopi = 2.0 * M_PI;
     double fftw_iter = total_time_/time_step_;
     for (long int i=1; i<fftw_iter; i++){
         double w = twopi*i/(fftw_iter*time_step_);
         if (w*pc_hartree2ev>max_freq) break;

         double valr = corr_func_x[i][0]/fftw_iter;
         double vali = corr_func_x[i][1]/fftw_iter;
         double valx  = valr*valr + vali*vali;

         valr = corr_func_y[i][0]/fftw_iter;
         vali = corr_func_y[i][1]/fftw_iter;
         double valy  = valr*valr + vali*vali;

         valr = corr_func_z[i][0]/fftw_iter;
         vali = corr_func_z[i][1]/fftw_iter;
         double valz  = valr*valr + vali*vali;

         outfile->Printf("      @Frequency %20.12le %20.12le %20.12le %20.12le \n", w*pc_hartree2ev, valx,valy, valz);
         
     }
}
void TD2RDM::DiagonalizeHermitianMatrix(long int N,double*re,double*im,double*W){
    WRAP(re,im,W,N);
}
struct tpdm {
    int i;
    int j;
    int k;
    int l;
    double val;
};

void TD2RDM::ReadTPDM(){

//    int o = ndoccact_;
//    int v = nvirt_;
    int nmo_ = nso_;

    std::shared_ptr<PSIO> psio (new PSIO());
    if ( !psio->exists(PSIF_V2RDM_D2AB) ) return;
    if ( !psio->exists(PSIF_V2RDM_D2AA) ) return;
    if ( !psio->exists(PSIF_V2RDM_D2BB) ) return;

    //Ca_->print();

    D2aa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    D2bb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    psio_address addr_aa = PSIO_ZERO;
    psio_address addr_bb = PSIO_ZERO;
    psio_address addr_ab = PSIO_ZERO;

    // ab
    psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

    long int nab;
    psio->read_entry(PSIF_V2RDM_D2AB,"length",(char*)&nab,sizeof(long int));

    for (int n = 0; n < nab; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2ab[id] = d2.val;
    }
    psio->close(PSIF_V2RDM_D2AB,1);
    // aa
    psio->open(PSIF_V2RDM_D2AA,PSIO_OPEN_OLD);

    long int naa;
    psio->read_entry(PSIF_V2RDM_D2AA,"length",(char*)&naa,sizeof(long int));

    for (int n = 0; n < naa; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2aa[id] = d2.val;
    }
    psio->close(PSIF_V2RDM_D2AA,1);

    // bb
    psio->open(PSIF_V2RDM_D2BB,PSIO_OPEN_OLD);

    long int nbb;
    psio->read_entry(PSIF_V2RDM_D2BB,"length",(char*)&nbb,sizeof(long int));

    for (int n = 0; n < nbb; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2bb[id] = d2.val;
    }
    psio->close(PSIF_V2RDM_D2BB,1);

   // check traces:
    double traa = 0.0;
    double trbb = 0.0;
    double trab = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            traa += D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
            trbb += D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
            trab += D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
        }
    }

    Da = (double*)malloc(nmo_*nmo_*sizeof(double));
    Db = (double*)malloc(nmo_*nmo_*sizeof(double));

    memset((void*)Da,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)Db,'\0',nmo_*nmo_*sizeof(double));

    double tra = 0.0;
    double trb = 0.0;

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {

            double duma = 0.0;
            double dumb = 0.0;
            for (int k = 0; k < nmo_; k++) {
                duma += D2ab[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];
                duma += D2aa[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];

                dumb += D2ab[k*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+j];
                dumb += D2bb[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];
            }
            Da[i*nmo_+j] = 1.0/(nalpha_+nbeta_-1.0) * duma;
            Db[i*nmo_+j] = 1.0/(nalpha_+nbeta_-1.0) * dumb;

            if ( i == j ) {
                tra += Da[i*nmo_+j];
                trb += Db[i*nmo_+j];
            }

        }
    }
}




}} // End namespaces

