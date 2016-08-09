#include"psi4-dec.h"
#include<liboptions/liboptions.h>
#include<libmints/mints.h>
#include<libpsio/psio.h>
#include<physconst.h>
#include<libqt/qt.h>
#include<../bin/fnocc/blas.h>
#include"frozen_natural_orbitals.h"
#include<psifiles.h>

#include"tdhf.h"

// boost numerical integrators live here:
#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>


#ifdef _OPENMP
    #include<omp.h>
#endif

using namespace psi;
using namespace fnocc;

namespace psi{ namespace tdhf_cqed {

boost::shared_ptr<TDHF> MyTDHF;

void is_this_necessary(boost::shared_ptr<Wavefunction> wfn, Options & options) {
    MyTDHF = (boost::shared_ptr<TDHF>) (new TDHF(wfn,options) );
    MyTDHF->compute_energy();
}

void rk4_call( state_type &x , state_type &dxdt , double t ){
    MyTDHF->rk4_call_gah(x,dxdt,t);
}

TDHF::TDHF(boost::shared_ptr<Wavefunction> reference_wavefunction,Options & options):
  Wavefunction(options)
{
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}
TDHF::~TDHF(){
}
void TDHF::common_init(){

//    printf("Initializing ...\n");
    escf    = reference_wavefunction_->reference_energy();
    doccpi_ = reference_wavefunction_->doccpi();
    soccpi_ = reference_wavefunction_->soccpi();
    frzcpi_ = reference_wavefunction_->frzcpi();
    frzvpi_ = reference_wavefunction_->frzvpi();
    nmopi_  = reference_wavefunction_->nmopi();
    nsopi_  = reference_wavefunction_->nsopi();
    molecule_ = reference_wavefunction_->molecule();
    nirrep_ = reference_wavefunction_->nirrep();

//    printf("Getting coefficients ...\n");

    Da_ = SharedMatrix(reference_wavefunction_->Da());
    Ca_ = SharedMatrix(reference_wavefunction_->Ca());
    Fa_ = SharedMatrix(reference_wavefunction_->Fa());

//    printf("Getting orbital parameters ... %i	%i \n",nirrep_, nsopi_);

    epsilon_a_= boost::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
//    printf("here1 ...\n");
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
//    printf("here2 ...\n");
    nalpha_ = reference_wavefunction_->nalpha();
//    printf("here3 ...\n");
    nbeta_  = reference_wavefunction_->nbeta();
    nso = nmo = ndocc = nvirt = nfzc = nfzv = 0;
//    printf("Assigning values ...\n");
    for (int h=0; h<nirrep_; h++){
        nfzc   += frzcpi_[h];
        nfzv   += frzvpi_[h];
        nso    += nsopi_[h];
        nmo    += nmopi_[h]-frzcpi_[h]-frzvpi_[h];
        ndocc  += doccpi_[h];
    }
//    printf("Exiting loop ...\n");
    ndoccact = ndocc - nfzc;
    nvirt    = nmo - ndoccact;

    if ( nfzc > 0 ) {
        throw PsiException("TDHF does not work with frozen core (yet).",__FILE__,__LINE__);
    }
    if ( nso != nmo ) {
        throw PsiException("TDHF does not work with nmo != nso (yet).",__FILE__,__LINE__);
    }
//    printf("Getting memory ...\n");

    // memory is from process::environment
    memory = Process::environment.get_memory();
    // set the wavefunction name
    name_ = "TDHF";

    // orbital energies
    eps = (double*)malloc(nmo*sizeof(double));
    memset((void*)eps,'\0',nmo*sizeof(double));
    int count=0;
    for (int h=0; h<nirrep_; h++){
        for (int norb = frzcpi_[h]; norb<doccpi_[h]; norb++){
            eps[count++] = epsilon_a_->get(h,norb);
        }
    }
    for (int h=0; h<nirrep_; h++){
        for (int norb = doccpi_[h]; norb<nmopi_[h]-frzvpi_[h]; norb++){
            eps[count++] = epsilon_a_->get(h,norb);
        //    printf("epsilon %i  = %lf \n",count, eps[count]);
        }
    }
    
    //printf("Computing the Kinetic and Potential energies ...\n");
    boost::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));
    T   = mints->so_kinetic();
    V   = mints->so_potential();

    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],T->pointer(),Ca_->pointer());
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],V->pointer(),Ca_->pointer());

    // if freezing the core, need to add frozen core contributions to the one-electron integrals:
    //TransformIntegralsFull();
    TransformIntegrals();

    // testing 4-index integrals:
    tei = (double*)malloc(nmo*nmo*nmo*nmo*sizeof(double));
    memset((void*)tei,'\0',nmo*nmo*nmo*nmo*sizeof(double));
    F_DGEMM('n','t',nmo*nmo,nmo*nmo,nQ,2.0,Qmo,nmo*nmo,Qmo,nmo*nmo,0.0,tei,nmo*nmo);
    #pragma omp parallel for schedule (static)
    for (int p = 0; p < nmo; p++) {
        for (int q = 0; q < nmo; q++) {
            for (int r = 0; r < nmo; r++) {
                for (int s = 0; s < nmo; s++) {
                    //for (int Q = 0; Q < nQ; Q++) {
                        //tei[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] += 2.0 * Qmo[Q*nmo*nmo+p*nmo+q]*Qmo[Q*nmo*nmo+r*nmo+s];
                        //tei[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] -=       Qmo[Q*nmo*nmo+p*nmo+r]*Qmo[Q*nmo*nmo+s*nmo+q];
                    //}
                    tei[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] -=       C_DDOT(nQ,&Qmo[p*nmo+r],nmo*nmo,&Qmo[s*nmo+q],nmo*nmo);
                }
            }
        }
    }
   
    nS = options_.get_int("N_PLASMON_STATES");

    int offset = 0;
    offset_dre           = offset; offset += nmo*nmo;
    offset_dim           = offset; offset += nmo*nmo;
    offset_dre_plasmon_x = offset; offset += nS*nS;
    offset_dim_plasmon_x = offset; offset += nS*nS;
    offset_dre_plasmon_y = offset; offset += nS*nS;
    offset_dim_plasmon_y = offset; offset += nS*nS;
    offset_dre_plasmon_z = offset; offset += nS*nS;
    offset_dim_plasmon_z = offset; offset += nS*nS;

    Dre        = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Dim        = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));

    Dre_plasmon_x = (boost::shared_ptr<Matrix>) (new Matrix(nS,nS));
    Dim_plasmon_x = (boost::shared_ptr<Matrix>) (new Matrix(nS,nS));
       
    Dre_plasmon_y = (boost::shared_ptr<Matrix>) (new Matrix(nS,nS));
    Dim_plasmon_y = (boost::shared_ptr<Matrix>) (new Matrix(nS,nS));
       
    Dre_plasmon_z = (boost::shared_ptr<Matrix>) (new Matrix(nS,nS));
    Dim_plasmon_z = (boost::shared_ptr<Matrix>) (new Matrix(nS,nS));

    F1re       = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    F1im       = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Fre        = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Fim        = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Fre_temp   = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    Fim_temp   = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));

    D1_e_re    = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    D1_e_im    = (boost::shared_ptr<Matrix>) (new Matrix(nmo,nmo));
    D1_p_re    = (boost::shared_ptr<Matrix>) (new Matrix(nS,nS));
    D1_p_im    = (boost::shared_ptr<Matrix>) (new Matrix(nS,nS));


    // new hamiltonians for plasmon part
    long int nS_scf  = options_.get_int("N_SCF_PLASMON_STATES");
    Dip_x              = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Dip_y              = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Dip_z              = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Hp_x               = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Hp_y               = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Hp_z               = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Hp_int_x           = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Hp_int_y           = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Hp_int_z           = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Hplasmon_total_x   = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Hplasmon_total_y   = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Hplasmon_total_z   = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Eigvec_x           = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Eigvec_y           = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Eigvec_z           = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    htemp_x            = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    htemp_y            = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    htemp_z            = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    htemp_int_x        = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    htemp_int_y        = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    htemp_int_z        = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    htemp_dip_x        = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    htemp_dip_y        = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    htemp_dip_z        = (boost::shared_ptr<Matrix>) (new Matrix(nS_scf,nS_scf));
    Eigval_x           = (boost::shared_ptr<Vector>) (new Vector(nS_scf));
    Eigval_y           = (boost::shared_ptr<Vector>) (new Vector(nS_scf));
    Eigval_z           = (boost::shared_ptr<Vector>) (new Vector(nS_scf));
    temp_x = (boost::shared_ptr<Matrix>) (new Matrix(nS,nS));
    temp_y = (boost::shared_ptr<Matrix>) (new Matrix(nS,nS));
    temp_z = (boost::shared_ptr<Matrix>) (new Matrix(nS,nS));

    for (int i = 0; i < ndoccact; i++) {
        Dre->pointer()[i][i] = 1.0;
    }
    for (int i = 0; i < nS; i++) {
        Dre_plasmon_x->pointer()[i][i] = 1.0;
        Dre_plasmon_y->pointer()[i][i] = 1.0;
        Dre_plasmon_z->pointer()[i][i] = 1.0;
    }

    //printf("Getting dipole integrals ...\n");
    // get dipole integrals:
    dipole = mints->so_dipole();
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],dipole[0]->pointer(),Ca_->pointer());
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],dipole[1]->pointer(),Ca_->pointer());
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],dipole[2]->pointer(),Ca_->pointer());


    //printf("Getting polarization ...\n");
    // get polarization:
    polarization = (double*)malloc(sizeof(double)*3);
    memset((void*)polarization,'\0',3*sizeof(double));
    if (options_["POLARIZATION"].has_changed()){
       if (options_["POLARIZATION"].size() != 3)
          throw PsiException("The POLARIZATION array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) polarization[i] = options_["POLARIZATION"][i].to_double();
       int pol = 0;
       for (int i = 0; i < 3; i++) {
           if (fabs(polarization[i]) > 1e-9) {
               pol++;
           }
       }
       if (pol > 1) {
          throw PsiException("only plane polarized light supported",__FILE__,__LINE__);
       }
    }else{
       polarization[0] = 0.0;
       polarization[1] = 0.0;
       polarization[2] = 1.0;
    }

    // get plasmon coordinates:
    plasmon_coordinates = (double*)malloc(sizeof(double)*3);
    memset((void*)plasmon_coordinates,'\0',3*sizeof(double));
    if (options_["PLASMON_COORDINATES"].has_changed()){
       if (options_["PLASMON_COORDINATES"].size() != 3)
          throw PsiException("The PLASMON COORDINATES array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) plasmon_coordinates[i] = options_["PLASMON_COORDINATES"][i].to_double();
    }else{
       plasmon_coordinates[0] = 0.0;
       plasmon_coordinates[1] = 0.0;
       plasmon_coordinates[2] = 100.0;
    }

    // get plasmon tdm:
    plasmon_tdm = (double*)malloc(sizeof(double)*3);
    memset((void*)plasmon_tdm,'\0',3*sizeof(double));
    if (options_["PLASMON_TDM"].has_changed()){
       if (options_["PLASMON_TDM"].size() != 3)
          throw PsiException("The PLASMON TDM array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) plasmon_tdm[i] = options_["PLASMON_TDM"][i].to_double();
    }else{
       plasmon_tdm[0] = 2990.0/2.54175;
       plasmon_tdm[1] = 2990.0/2.54175;
       plasmon_tdm[2] = 2990.0/2.54175;
    }

    plasmon_tdm_x = plasmon_tdm[0]*polarization[0]; 
    plasmon_tdm_y = plasmon_tdm[1]*polarization[1]; 
    plasmon_tdm_z = plasmon_tdm[2]*polarization[2]; 

    // Getting molecular information to compute the center of mass
    boost::shared_ptr<Molecule>mol=Process::environment.molecule();
    natom_ = mol->natom();
    com_ = (double*)malloc(sizeof(double)*3);
    memset((void*)com_,'\0',3*sizeof(double));
    
    double temp_x = 0.0;
    double temp_y = 0.0;
    double temp_z = 0.0;
    double temp_m = 0.0;
    for (int i = 0; i < natom_ ; i++){
       temp_x += mol->mass(i) * mol->x(i);         
       temp_y += mol->mass(i) * mol->y(i);         
       temp_z += mol->mass(i) * mol->z(i);         
       temp_m += mol->mass(i);
    }

    com_[0] = temp_x/temp_m;
    com_[1] = temp_y/temp_m;
    com_[2] = temp_z/temp_m;

    // nuclear contribution to dipole moment:
    nuc_dip_x_ = 0.0;
    nuc_dip_y_ = 0.0;
    nuc_dip_z_ = 0.0;
    for (int i = 0; i < natom_; i++) {
        nuc_dip_x_ += mol->Z(i) * mol->x(i);
        nuc_dip_y_ += mol->Z(i) * mol->y(i);
        nuc_dip_z_ += mol->Z(i) * mol->z(i);
    }

    nuc_dip = nuc_dip_x_ + nuc_dip_y_ + nuc_dip_z_ ;
    //printf("Center of mass X = %lf	Y = %lf	Z = %lf \n",com_[0],com_[1],com_[2]);

    
    // Plasmonic dipole moment operator

    Dip_x->zero();
    Dip_y->zero();
    Dip_z->zero();
    for (int s=0; s<nS_scf-1; s++){
        Dip_x->pointer()[s+1][s] += plasmon_tdm_x*sqrt(s+1);
        Dip_x->pointer()[s][s+1] += plasmon_tdm_x*sqrt(s+1);

        Dip_y->pointer()[s+1][s] += plasmon_tdm_y*sqrt(s+1);
        Dip_y->pointer()[s][s+1] += plasmon_tdm_y*sqrt(s+1);

        Dip_z->pointer()[s+1][s] += plasmon_tdm_z*sqrt(s+1);
        Dip_z->pointer()[s][s+1] += plasmon_tdm_z*sqrt(s+1);
    }

    Dip_x->print();
    Dip_y->print();
    Dip_z->print();
   
    //exit(0); 
    // Plasmon Hamiltonian

    plasmon_e = (double*)malloc(sizeof(double)*3);
    memset((void*)plasmon_e,'\0',3*sizeof(double));
    if (options_["PLASMON_E"].has_changed()){
       if (options_["PLASMON_E"].size() != 3)
          throw PsiException("The PLASMON E array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) plasmon_e[i] = options_["PLASMON_E"][i].to_double();
    }else{
       plasmon_e[0] = 2.042/27.21138;  // Energy for the Au nanoparticle taken from Gray's paper
       plasmon_e[1] = 2.042/27.21138;
       plasmon_e[2] = 2.042/27.21138;
    }

    for (int A=0; A<nS_scf; A++){
        Hp_x->pointer()[A][A] = A*plasmon_e[0];
        Hp_y->pointer()[A][A] = A*plasmon_e[1];
        Hp_z->pointer()[A][A] = A*plasmon_e[2];
    }

    // Interaction Hamiltonian

    double r = (com_[0]-plasmon_coordinates[0])*(com_[0]-plasmon_coordinates[0])
             + (com_[1]-plasmon_coordinates[1])*(com_[1]-plasmon_coordinates[1])
             + (com_[2]-plasmon_coordinates[2])*(com_[2]-plasmon_coordinates[2]);
    r = sqrt(r);

    double delta_x = plasmon_coordinates[0] - com_[0];
    double delta_y = plasmon_coordinates[1] - com_[1];
    double delta_z = plasmon_coordinates[2] - com_[2];

    double * r_vector;

    r_vector = (double*)malloc(sizeof(double)*3);
    memset((void*)r_vector,'\0',3*sizeof(double));

    r_vector[0] = delta_x;
    r_vector[1] = delta_y;
    r_vector[2] = delta_z;

    double oer3 = 1.0 /(r*r*r);
    coupling_strength = 1.0 * oer3;

    plasmon_tdm_x = polarization[0]*plasmon_tdm[0];
    plasmon_tdm_y = polarization[1]*plasmon_tdm[1];
    plasmon_tdm_z = polarization[2]*plasmon_tdm[2];

    e_dip_x = 2.0*C_DDOT(nso*nso,&(Dre->pointer())[0][0],1,&(dipole[0]->pointer())[0][0],1);
    e_dip_y = 2.0*C_DDOT(nso*nso,&(Dre->pointer())[0][0],1,&(dipole[1]->pointer())[0][0],1);
    e_dip_z = 2.0*C_DDOT(nso*nso,&(Dre->pointer())[0][0],1,&(dipole[2]->pointer())[0][0],1);

    for (int A = 0; A < nS_scf; A++) {
        if (A < nS_scf - 1) {
            Hp_int_x->pointer()[A][A+1] += (e_dip_x + nuc_dip_x_)*plasmon_tdm[0];
            Hp_int_x->pointer()[A][A+1] -= 3.0*plasmon_tdm[0]*delta_x*(r_vector[0]*(nuc_dip_x_ + e_dip_x)+r_vector[1]*(nuc_dip_y_ + e_dip_y)+r_vector[2]*(nuc_dip_z_ + e_dip_z))/(r*r);
            Hp_int_x->pointer()[A][A+1] *= coupling_strength*sqrt(A+1);

            Hp_int_y->pointer()[A][A+1] += (e_dip_y + nuc_dip_y_)*plasmon_tdm[1];
            Hp_int_y->pointer()[A][A+1] -= 3.0*plasmon_tdm[1]*delta_y*(r_vector[0]*(nuc_dip_x_ + e_dip_x)+r_vector[1]*(nuc_dip_y_ + e_dip_y)+r_vector[2]*(nuc_dip_z_ + e_dip_z))/(r*r);
            Hp_int_y->pointer()[A][A+1] *= coupling_strength*sqrt(A+1);

            Hp_int_z->pointer()[A][A+1] += (e_dip_z + nuc_dip_z_)*plasmon_tdm[2];
            Hp_int_z->pointer()[A][A+1] -= 3.0*plasmon_tdm[2]*delta_z*(r_vector[0]*(nuc_dip_x_ + e_dip_x)+r_vector[1]*(nuc_dip_y_ + e_dip_y)+r_vector[2]*(nuc_dip_z_ + e_dip_z))/(r*r);
            Hp_int_z->pointer()[A][A+1] *= coupling_strength*sqrt(A+1);
        }
        if (A > 0) {
            Hp_int_x->pointer()[A][A-1] += (e_dip_x + nuc_dip_x_)*plasmon_tdm[0];
            Hp_int_x->pointer()[A][A-1] -= 3.0*plasmon_tdm[0]*delta_x*(r_vector[0]*(nuc_dip_x_ + e_dip_x)+r_vector[1]*(nuc_dip_y_ + e_dip_y)+r_vector[2]*(nuc_dip_z_ + e_dip_z))/(r*r);
            Hp_int_x->pointer()[A][A-1] *= coupling_strength*sqrt(A);

            Hp_int_y->pointer()[A][A-1] += (e_dip_y + nuc_dip_y_)*plasmon_tdm[1];
            Hp_int_y->pointer()[A][A-1] -= 3.0*plasmon_tdm[1]*delta_y*(r_vector[0]*(nuc_dip_x_ + e_dip_x)+r_vector[1]*(nuc_dip_y_ + e_dip_y)+r_vector[2]*(nuc_dip_z_ + e_dip_z))/(r*r);
            Hp_int_y->pointer()[A][A-1] *= coupling_strength*sqrt(A);

            Hp_int_z->pointer()[A][A-1] += (e_dip_z + nuc_dip_z_)*plasmon_tdm[2];
            Hp_int_z->pointer()[A][A-1] -= 3.0*plasmon_tdm[2]*delta_z*(r_vector[0]*(nuc_dip_x_ + e_dip_x)+r_vector[1]*(nuc_dip_y_ + e_dip_y)+r_vector[2]*(nuc_dip_z_ + e_dip_z))/(r*r);
            Hp_int_z->pointer()[A][A-1] *= coupling_strength*sqrt(A);
        } 
    }
 
    // Diagonalize total plasmon Hamiltonian

    Hplasmon_total_x -> copy(Hp_x);
    Hplasmon_total_x -> add(Hp_int_x);

    Hplasmon_total_y -> copy(Hp_y);
    Hplasmon_total_y -> add(Hp_int_y);

    Hplasmon_total_z -> copy(Hp_z);
    Hplasmon_total_z -> add(Hp_int_z);

    Hplasmon_total_x -> diagonalize(Eigvec_x, Eigval_x);
    Hplasmon_total_y -> diagonalize(Eigvec_y, Eigval_y);
    Hplasmon_total_z -> diagonalize(Eigvec_z, Eigval_z);
    
//    Eigval -> print();
//    Eigvec -> print();   

    htemp_x->zero(); 
    htemp_int_x->zero(); 
    htemp_dip_x->zero(); 
    htemp_y->zero(); 
    htemp_int_y->zero(); 
    htemp_dip_y->zero(); 
    htemp_z->zero(); 
    htemp_int_z->zero(); 
    htemp_dip_z->zero(); 
    for (int A=0; A < nS_scf; A++){
        for (int B=0; B < nS_scf; B++){
            for (int C=0; C < nS_scf; C++){
                htemp_x->pointer()[A][B] += Hp_x->pointer()[A][C]*Eigvec_x->pointer()[C][B];
                htemp_int_x->pointer()[A][B] += Hp_int_x->pointer()[A][C]*Eigvec_x->pointer()[C][B];
                htemp_dip_x->pointer()[A][B] += Dip_x->pointer()[A][C]*Eigvec_x->pointer()[C][B];

                htemp_y->pointer()[A][B] += Hp_y->pointer()[A][C]*Eigvec_y->pointer()[C][B];
                htemp_int_y->pointer()[A][B] += Hp_int_y->pointer()[A][C]*Eigvec_y->pointer()[C][B];
                htemp_dip_y->pointer()[A][B] += Dip_y->pointer()[A][C]*Eigvec_y->pointer()[C][B];

                htemp_z->pointer()[A][B] += Hp_z->pointer()[A][C]*Eigvec_z->pointer()[C][B];
                htemp_int_z->pointer()[A][B] += Hp_int_z->pointer()[A][C]*Eigvec_z->pointer()[C][B];
                htemp_dip_z->pointer()[A][B] += Dip_z->pointer()[A][C]*Eigvec_z->pointer()[C][B];
            }
        }
    }

//    htemp -> print();
    Hp_x->zero();
    Hp_int_x->zero();
    Dip_x->zero();

    Hp_y->zero();
    Hp_int_y->zero();
    Dip_y->zero();

    Hp_z->zero();
    Hp_int_z->zero();
    Dip_z->zero();

    for (int A=0; A < nS_scf; A++){
        for (int B=0; B < nS_scf; B++){
            for (int C=0; C < nS_scf; C++){
                Hp_x->pointer()[A][B] += Eigvec_x->pointer()[C][A]*htemp_x->pointer()[C][B];
                Hp_int_x->pointer()[A][B] += Eigvec_x->pointer()[C][A]*htemp_int_x->pointer()[C][B];
                Dip_x->pointer()[A][B] += Eigvec_x->pointer()[C][A]*htemp_dip_x->pointer()[C][B];

                Hp_y->pointer()[A][B] += Eigvec_y->pointer()[C][A]*htemp_y->pointer()[C][B];
                Hp_int_y->pointer()[A][B] += Eigvec_y->pointer()[C][A]*htemp_int_y->pointer()[C][B];
                Dip_y->pointer()[A][B] += Eigvec_y->pointer()[C][A]*htemp_dip_y->pointer()[C][B];

                Hp_z->pointer()[A][B] += Eigvec_z->pointer()[C][A]*htemp_z->pointer()[C][B];
                Hp_int_z->pointer()[A][B] += Eigvec_z->pointer()[C][A]*htemp_int_z->pointer()[C][B];
                Dip_z->pointer()[A][B] += Eigvec_z->pointer()[C][A]*htemp_dip_z->pointer()[C][B];
            }
        }
    }

    Hp_x->print();
    Hp_y->print();
    Hp_z->print();

    //Hplasmon_total -> print(); 
    //Hp -> print();
    //Dip -> print();
    //Hp -> print();
    //Hp_int -> print();
    //Hplasmon_total -> print();  
    //return;

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
    }else if (options_.get_str("LASER_SHAPE") == "GAUSSIAN"){
        // gaussian pulse
        pulse_shape_ = 4;
    }


    linear_response = false;

    // pad the correlation function with zeros just to get more output points
    //extra_pts = 4*(ttot/time_step+2);
    extra_pts = 0; //100000;//0;//1000000;
    // correlation function or dipole acceleration (fourier transformed)
    midpt = total_time/time_step+extra_pts + 1;
    corr_func = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+2*total_time/time_step+2+2*extra_pts));
    corr_func2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+2*total_time/time_step+2+2*extra_pts));
    corr_func3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+2*total_time/time_step+2+2*extra_pts));
    corr_func4 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+2*total_time/time_step+2+2*extra_pts));
    corr_func5 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+2*total_time/time_step+2+2*extra_pts));
    corr_func6 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+2*total_time/time_step+2+2*extra_pts));
    corr_func7 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+2*total_time/time_step+2+2*extra_pts));
    corr_func8 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+2*total_time/time_step+2+2*extra_pts));
    corr_func9 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(1+2*total_time/time_step+2+2*extra_pts));
  
    // maximum frequency to output (eV)
    max_freq = 30.0;

    // stencil for second derivative of dipole moment
    stencil = (double*)malloc(sizeof(double)*5);
    memset((void*)stencil,'\0',5*sizeof(double));

    // dipole potential integrals:
    DipolePotentialIntegrals();



}

void TDHF::DipolePotentialIntegrals() {

    boost::shared_ptr<OneBodyAOInt> efp_ints(reference_wavefunction_->integral()->ao_efp_multipole_potential());

    int nbf = reference_wavefunction_->basisset()->nbf();
    int nao = reference_wavefunction_->basisset()->nao();

    std::vector< boost::shared_ptr<Matrix> > mats;
    for(int i=0; i < 20; i++) {
        mats.push_back(boost::shared_ptr<Matrix> (new Matrix(nao, nao)));
        mats[i]->zero();
    }

    // plasmon dipole potential felt by the molecule
    Vector3 coords(plasmon_coordinates[0],plasmon_coordinates[1],plasmon_coordinates[2]);
    efp_ints->set_origin(coords);
    efp_ints->compute(mats);

    // ao/so transformation
    boost::shared_ptr<PetiteList> pet(new PetiteList(reference_wavefunction_->basisset(),reference_wavefunction_->integral(),true));
    boost::shared_ptr<Matrix> U = pet->aotoso();

    boost::shared_ptr<Matrix> Vx = Matrix::triplet(U,mats[1],U,true,false,false);
    boost::shared_ptr<Matrix> Vy = Matrix::triplet(U,mats[2],U,true,false,false);
    boost::shared_ptr<Matrix> Vz = Matrix::triplet(U,mats[3],U,true,false,false);

    // This function is important until here!! For the SCF procedure.

    // so/mo transformation
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],Vx->pointer(),Ca_->pointer());
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],Vy->pointer(),Ca_->pointer());
    SoToMo(Ca_->rowspi()[0],Ca_->colspi()[0],Vz->pointer(),Ca_->pointer());

    dipole_pot_x = (boost::shared_ptr<Matrix>) (new Matrix(Vx));
    dipole_pot_y = (boost::shared_ptr<Matrix>) (new Matrix(Vy));
    dipole_pot_z = (boost::shared_ptr<Matrix>) (new Matrix(Vz));

}

// TODO: get rid of all of the extra containers to hold F
void TDHF::BuildFockThreeIndex(double * Dre_temp, double * Dim_temp, bool use_oe_terms) {

    Fre_temp->zero();
    Fim_temp->zero();

    // J
    F_DGEMV('t',nmo*nmo,nQ,2.0,Qmo,nmo*nmo,Dre_temp,1,0.0,Ire,1);
    F_DGEMV('t',nmo*nmo,nQ,2.0,Qmo,nmo*nmo,Dim_temp,1,0.0,Iim,1);

    F_DGEMV('n',nmo*nmo,nQ,0.5,Qmo,nmo*nmo,Ire,1,0.0,&(Fre_temp->pointer()[0][0]),1);
    F_DGEMV('n',nmo*nmo,nQ,0.5,Qmo,nmo*nmo,Iim,1,0.0,&(Fim_temp->pointer()[0][0]),1);

    if ( use_oe_terms ) {
        //Fre_temp->add(T);
        //Fre_temp->add(V);
        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                Fre_temp->pointer()[i][j] += T->pointer()[i+nfzc][j+nfzc];
                Fre_temp->pointer()[i][j] += V->pointer()[i+nfzc][j+nfzc];
            }
        }
    }
    
    // K(re)
    F_DGEMM('t','n',nmo,nmo*nQ,nmo,1.0,Dre_temp,nmo,Qmo,nmo,0.0,Ire,nmo);
    for (int q = 0; q < nQ; q++) {
        for (int a = 0; a < nmo; a++) {
            for (int j = 0; j < nmo; j++) {
                Iim[a*nmo+j] = Ire[q*nmo*nmo+j*nmo+a];
            }
        }
        C_DCOPY(nmo*nmo,Iim,1,Ire+q*nmo*nmo,1);
    }
    F_DGEMM('n','t',nmo,nmo,nQ*nmo,-1.0 * 0.5,Ire,nmo,Qmo,nmo,1.0,&(Fre_temp->pointer()[0][0]),nmo);
    //F_DGEMM('n','t',nmo,nmo,nQ*nmo,-1.0,Ire,nmo,Qmo,nmo,1.0,&(Fre_temp->pointer()[0][0]),nmo);
    // K(im)
    F_DGEMM('t','n',nmo,nmo*nQ,nmo,1.0,Dim_temp,nmo,Qmo,nmo,0.0,Ire,nmo);
    for (int q = 0; q < nQ; q++) {
        for (int a = 0; a < nmo; a++) {
            for (int j = 0; j < nmo; j++) {
                Iim[a*nmo+j] = Ire[q*nmo*nmo+j*nmo+a];
            }
        }
        C_DCOPY(nmo*nmo,Iim,1,Ire+q*nmo*nmo,1);
    }
    F_DGEMM('n','t',nmo,nmo,nQ*nmo,-1.0 * 0.5,Ire,nmo,Qmo,nmo,1.0,&(Fim_temp->pointer()[0][0]),nmo);
    //F_DGEMM('n','t',nmo,nmo,nQ*nmo,-1.0,Ire,nmo,Qmo,nmo,1.0,&(Fim_temp->pointer()[0][0]),nmo);

    if ( use_oe_terms ) {
        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                Fre_temp->pointer()[i][j] -= ext_field * dipole[0]->pointer()[i+nfzc][j+nfzc] * polarization[0];
                Fre_temp->pointer()[i][j] -= ext_field * dipole[1]->pointer()[i+nfzc][j+nfzc] * polarization[1];
                Fre_temp->pointer()[i][j] -= ext_field * dipole[2]->pointer()[i+nfzc][j+nfzc] * polarization[2];
            }
        }
    }
}


// TODO: get rid of all of the extra containers to hold F
// Fij = sum_kl Dkl ( 2 (ij|kl) - (il|kj) )
// Fij = sum_kl Dkl t(kl|ij)
void TDHF::BuildFock(double * Dre_temp, double * Dim_temp, bool use_oe_terms) {

    Fre_temp->zero();
    Fim_temp->zero();

    // TODO: out-of-core version
    F_DGEMV('t',nmo*nmo,nmo*nmo,1.0,tei,nmo*nmo,Dre_temp,1,0.0,&(Fre_temp->pointer()[0][0]),1);
    F_DGEMV('t',nmo*nmo,nmo*nmo,1.0,tei,nmo*nmo,Dim_temp,1,0.0,&(Fim_temp->pointer()[0][0]),1);

    if ( use_oe_terms ) {
        //Fre_temp->add(T);
        //Fre_temp->add(V);
        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                Fre_temp->pointer()[i][j] += T->pointer()[i+nfzc][j+nfzc];
                Fre_temp->pointer()[i][j] += V->pointer()[i+nfzc][j+nfzc];
            }
        }
    }

    if ( use_oe_terms ) {
        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                Fre_temp->pointer()[i][j] -= ext_field * dipole[0]->pointer()[i+nfzc][j+nfzc] * polarization[0];
                Fre_temp->pointer()[i][j] -= ext_field * dipole[1]->pointer()[i+nfzc][j+nfzc] * polarization[1];
                Fre_temp->pointer()[i][j] -= ext_field * dipole[2]->pointer()[i+nfzc][j+nfzc] * polarization[2];
            }
        }
    }

}

void TDHF::ExtField(double curtime){ 

    //Vext->zero();
    double sigma = laser_time*0.5;
    if (!linear_response) {

        // add external field

        ext_field = 0.0;

        if (pulse_shape_ == 0 ) {

            // from prl:
            if ( curtime < laser_time ) {
                ext_field = sin(M_PI*curtime/(laser_time));
                ext_field *= ext_field*laser_amp*sin(laser_freq*curtime);
            }

        } else if ( pulse_shape_ == 1 ) {

            // from 2007 schlegel paper (jcp 126, 244110 (2007))
            if (curtime <= 2.0 * M_PI / laser_freq)      ext_field = laser_freq * curtime / (2.0 * M_PI) * laser_amp;
            else if (curtime <= 4.0 * M_PI / laser_freq) ext_field = laser_amp;
            else if (curtime <= 6.0 * M_PI / laser_freq) ext_field = (3.0 - laser_freq * curtime / (2.0 * M_PI) ) * laser_amp;
            ext_field *= sin(laser_freq*curtime);

        } else if ( pulse_shape_ == 2 ) {

            // pi pulse from licn paper
            double sigma = laser_time*0.5;
            if ( curtime < laser_time ) {
                ext_field = cos(M_PI*curtime/(2.0*sigma) + 0.5*M_PI);
                ext_field *= M_PI/(sigma*transition_dpm) * ext_field * laser_amp * cos(laser_freq*curtime);
            }

        } else if ( pulse_shape_ == 3 ) {

            // continuous wave for rabi flopping
            ext_field = laser_amp*sin(laser_freq*curtime); 
        } else if ( pulse_shape_ == 4 ) {

            // Gaussian pulse
            ext_field = laser_amp * exp(-((curtime-1.5*laser_time)*(curtime-1.5*laser_time))
                      / (0.3606738*laser_time*laser_time)) * sin(laser_freq*curtime);

        }
    }
}

void TDHF::InteractionContribution(double * tempr,
                                   double * tempi,
                                   double * kre,
                                   double * kim, 
                                   double * tempr_p,
                                   double * tempi_p,
                                   double * kre_p,
                                   double * kim_p, 
                                   boost::shared_ptr<Matrix> Hp_int,
                                   boost::shared_ptr<Matrix> dipole_pot,
                                   double mdip, 
                                   double pdip) {

    // contribution to plasmon:

    for (int A=0; A<nS; A++){
        for (int B=0; B<nS; B++){
            double dumr = 0.0;
            double dumi = 0.0;
            for (int C=0; C<nS; C++){

                 // remember, only excite one mode for now
                 kre_p[A*nS+B] -= Hp_int->pointer()[A][C]*tempi_p[C*nS+B];    
                 kre_p[A*nS+B] += Hp_int->pointer()[C][B]*tempi_p[A*nS+C];    
          
                 kim_p[A*nS+B] += Hp_int->pointer()[A][C]*tempr_p[C*nS+B];    
                 kim_p[A*nS+B] -= Hp_int->pointer()[C][B]*tempr_p[A*nS+C];    
            }
        }
    }

    // contribution to electron
    for (int i=0; i<nmo; i++){
        for (int j=0; j<nmo; j++){

            double dumr = 0.0;
            double dumi = 0.0;
            for (int q = 0; q < nmo; q++) {
                dumr += tempr[i*nmo+q] * dipole_pot->pointer()[j+nfzc][q+nfzc];
                dumr -= tempr[q*nmo+j] * dipole_pot->pointer()[q+nfzc][i+nfzc];

                dumi += tempi[i*nmo+q] * dipole_pot->pointer()[j+nfzc][q+nfzc];
                dumi -= tempi[q*nmo+j] * dipole_pot->pointer()[q+nfzc][i+nfzc];
            }

            kre[i*nmo+j] += dumi * pdip * (-1);
            kim[i*nmo+j] -= dumr * pdip * (-1);

        }
    }
}


void TDHF::PlasmonContribution(double * tempr,
                               double * tempi,
                               double * kre,
                               double * kim, 
                               boost::shared_ptr<Matrix> dip, 
                               boost::shared_ptr<Matrix> Ham, 
                               double pol) {


    // one-particle part of uncoupled plasmon hamiltonian: (TODO: rewrite as DGEMM)

    for (int A=0; A<nS; A++){
        for (int B=0; B<nS; B++){
            for (int C=0; C<nS; C++){
                 kre[A*nS+B] -= Ham->pointer()[A][C]*tempi[C*nS+B];    
                 kre[A*nS+B] += Ham->pointer()[C][B]*tempi[A*nS+C];    
          
                 kim[A*nS+B] += Ham->pointer()[A][C]*tempr[C*nS+B];
                 kim[A*nS+B] -= Ham->pointer()[C][B]*tempr[A*nS+C];
             }
        }
    }

    // external field part: (TODO: rewrite as DGEMM)

    for (int A=0; A<nS; A++){
        for (int B=0; B<nS; B++){
            for (int C=0; C<nS; C++){
                kre[A*nS+B] -= dip->pointer()[A][C]*tempi[C*nS+B]*(-ext_field)*pol;    
                kre[A*nS+B] += dip->pointer()[C][B]*tempi[A*nS+C]*(-ext_field)*pol;    
           
                kim[A*nS+B] += dip->pointer()[A][C]*tempr[C*nS+B]*(-ext_field)*pol;    
                kim[A*nS+B] -= dip->pointer()[C][B]*tempr[A*nS+C]*(-ext_field)*pol;    
            }
        }
    }

}

void TDHF::BuildLindblad(double * tempr,
                         double * tempi,
                         double * kre,
                         double * kim) {

    plasmon_dr = options_.get_double("PLASMON_DR");
    for (int s = 0; s < nS; s++){
        for (int p = 0; p < nS; p++){
            kre[s*nS+p] -= 0.5*plasmon_dr*tempr[s*nS+p]*(s+p);
            kim[s*nS+p] -= 0.5*plasmon_dr*tempi[s*nS+p]*(s+p);
            if (s < nS-1 && p < nS-1){
               kre[s*nS+p] += plasmon_dr*tempr[(s+1)*nS+(p+1)]*sqrt((s+1)*(p+1)); 
               kim[s*nS+p] += plasmon_dr*tempi[(s+1)*nS+(p+1)]*sqrt((s+1)*(p+1)); 
            }
        }
    }
}

// so->mo transformation for 1-body matrix
void TDHF::SoToMo(int nsotemp,int nmotemp,double**mat,double**trans){
  double*tmp = (double*)malloc(sizeof(double)*nsotemp*nsotemp);
  memset((void*)tmp,'\0',nsotemp*nsotemp*sizeof(double));
  F_DGEMM('n','n',nmotemp,nsotemp,nsotemp,1.0,&trans[0][0],nmotemp,&mat[0][0],nsotemp,0.0,&tmp[0],nmotemp);
  F_DGEMM('n','t',nmotemp,nmotemp,nsotemp,1.0,&tmp[0],nmotemp,&trans[0][0],nmotemp,0.0,&mat[0][0],nsotemp);
  free(tmp);
}

double TDHF::compute_energy() {

    // rk4
    state_type rk4_buffer(nmo*nmo*nS*nS*2);
    rk4_stepper rk4;
    //rk78_stepper rk4;

    fftw_iter   = 0;

    double * factorB = (double*)malloc(sizeof(double)*5);
    double * factorb = (double*)malloc(sizeof(double)*5);
    memset((void*)factorB,'\0',5*sizeof(double));
    memset((void*)factorb,'\0',5*sizeof(double));

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

        // RK4
        // y(n+1) = y( n ) + 1/6 h ( k1 + 2k2 + 2k3 + k4 )
        // t(n+1) = t( n ) + h

        //Dre-> print();
        //exit(0);

        // rk4 buffer should contain: Dre, Dim, Dpxre, Dpxim, Dpyre, Dpyim, Dpzre, Dpzim

         
        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                rk4_buffer[offset_dre + i*nmo+j] = Dre->pointer()[i][j];
                rk4_buffer[offset_dim + i*nmo+j] = Dim->pointer()[i][j];
            }
        }
        for (int i = 0; i < nS; i++) {
            for (int j = 0; j < nS; j++) {
                rk4_buffer[offset_dre_plasmon_x + i*nS+j] = Dre_plasmon_x->pointer()[i][j];
                rk4_buffer[offset_dim_plasmon_x + i*nS+j] = Dim_plasmon_x->pointer()[i][j];

                rk4_buffer[offset_dre_plasmon_y + i*nS+j] = Dre_plasmon_y->pointer()[i][j];
                rk4_buffer[offset_dim_plasmon_y + i*nS+j] = Dim_plasmon_y->pointer()[i][j];

                rk4_buffer[offset_dre_plasmon_z + i*nS+j] = Dre_plasmon_z->pointer()[i][j];
                rk4_buffer[offset_dim_plasmon_z + i*nS+j] = Dim_plasmon_z->pointer()[i][j];
            }
        }
        rk4.do_step( rk4_call , rk4_buffer , iter * time_step , time_step );
        for (int i = 0; i < nmo; i++) {
            for (int j = 0; j < nmo; j++) {
                Dre->pointer()[i][j] = rk4_buffer[offset_dre + i*nmo+j];
                Dim->pointer()[i][j] = rk4_buffer[offset_dim + i*nmo+j];
            }
        }
        for (int i = 0; i < nS; i++) {
            for (int j = 0; j < nS; j++) {
                Dre_plasmon_x->pointer()[i][j] = rk4_buffer[offset_dre_plasmon_x + i*nS+j];
                Dim_plasmon_x->pointer()[i][j] = rk4_buffer[offset_dim_plasmon_x + i*nS+j];

                Dre_plasmon_y->pointer()[i][j] = rk4_buffer[offset_dre_plasmon_y + i*nS+j];
                Dim_plasmon_y->pointer()[i][j] = rk4_buffer[offset_dim_plasmon_y + i*nS+j];

                Dre_plasmon_z->pointer()[i][j] = rk4_buffer[offset_dre_plasmon_z + i*nS+j];
                Dim_plasmon_z->pointer()[i][j] = rk4_buffer[offset_dim_plasmon_z + i*nS+j];
            }
        }

        // evaluate dipole moment:

        e_dip_x = 2.0*C_DDOT(nso*nso,&(Dre->pointer())[0][0],1,&(dipole[0]->pointer())[0][0],1);
        e_dip_y = 2.0*C_DDOT(nso*nso,&(Dre->pointer())[0][0],1,&(dipole[1]->pointer())[0][0],1);
        e_dip_z = 2.0*C_DDOT(nso*nso,&(Dre->pointer())[0][0],1,&(dipole[2]->pointer())[0][0],1);

        dipole_p_x = 0.0;
        dipole_p_y = 0.0;
        dipole_p_z = 0.0;

        // Using loops instead of DDOT
        
        temp_x->zero();
        temp_y->zero();
        temp_z->zero();

        for (int A=0; A<nS; A++){
            for (int B=0; B<nS; B++){
                for (int C=0; C<nS; C++){
                    temp_x->pointer()[A][B] += D1_p_re->pointer()[A][C]*Dip_x->pointer()[C][B];
                    temp_y->pointer()[A][B] += D1_p_re->pointer()[A][C]*Dip_y->pointer()[C][B];
                    temp_z->pointer()[A][B] += D1_p_re->pointer()[A][C]*Dip_z->pointer()[C][B];
                }
            }
        }

        for (int A=0; A<nS; A++){
            dipole_p_x += temp_x->pointer()[A][A];
            dipole_p_y += temp_y->pointer()[A][A];
            dipole_p_z += temp_z->pointer()[A][A];
        } 

        // start accumulating the dipole acceleration after the pulse is over
        //if (!linear_response && iter * time_step >2000 && iter *time_step < 2000+1.0 / laser_freq * 2.0 * M_PI)
        //if (!linear_response && iter * time_step >laser_time){
           corr_func[fftw_iter][0]  = e_dip_x;
           corr_func[fftw_iter][1]  = 0.0;

           corr_func2[fftw_iter][0] = e_dip_y;
           corr_func2[fftw_iter][1] = 0.0;

           corr_func3[fftw_iter][0] = e_dip_z;
           corr_func3[fftw_iter][1] = 0.0;
          
           corr_func4[fftw_iter][0] = dipole_p_x;
           corr_func4[fftw_iter][1] = 0.0;
       
           corr_func5[fftw_iter][0] = dipole_p_y;
           corr_func5[fftw_iter][1] = 0.0;
       
           corr_func6[fftw_iter][0] = dipole_p_z;
           corr_func6[fftw_iter][1] = 0.0;

           ExtField(iter*time_step);

           corr_func7[fftw_iter][0] = ext_field * polarization[0];
           corr_func7[fftw_iter][1] = 0.0;
       
           corr_func8[fftw_iter][0] = ext_field * polarization[1];
           corr_func8[fftw_iter][1] = 0.0;
       
           corr_func9[fftw_iter][0] = ext_field * polarization[2];
           corr_func9[fftw_iter][1] = 0.0;
       
           fftw_iter++;
        //}

        //CorrelationFunction();

        double en = Process::environment.molecule()->nuclear_repulsion_energy();

        fprintf(outfile,"@TIME %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le\n",iter*time_step,
            e_dip_x,e_dip_y,e_dip_z,
            dipole_p_x,dipole_p_y,dipole_p_z,
            ext_field * polarization[0],ext_field * polarization[1],ext_field * polarization[2]);

    }

    free(factorb);
    free(factorB);


    // fourier transform and spectrum
    FFTW();
    Spectrum();

    free(corr_func);
    free(corr_func2);
    free(corr_func3);
    free(corr_func4);
    free(corr_func5);
    free(corr_func6);
    free(corr_func7);
    free(corr_func8);
    free(corr_func9);
    

    return 0.0;
}

void TDHF::rk4_call_gah( state_type &x , state_type &dxdt , double t ){

    for (int i = 0; i < 2*(nmo*nmo+3*nS*nS); i++) {
        dxdt[i] = 0.0;
    }

    e_dip_x = 0.0;
    e_dip_y = 0.0;
    e_dip_z = 0.0;
    for (int i = 0; i < nmo; i++) {
        for (int j = 0; j < nmo; j++) {
            e_dip_x += x[offset_dre+i*nmo+j] * dipole[0]->pointer()[i+nfzc][j+nfzc];
            e_dip_y += x[offset_dre+i*nmo+j] * dipole[1]->pointer()[i+nfzc][j+nfzc];
            e_dip_z += x[offset_dre+i*nmo+j] * dipole[2]->pointer()[i+nfzc][j+nfzc];
        }
    }
    e_dip_x *= 2.0;
    e_dip_y *= 2.0;
    e_dip_z *= 2.0;

    dipole_p_x = 0.0;
    dipole_p_y = 0.0;
    dipole_p_z = 0.0;
         
    temp_x->zero();
    temp_y->zero();
    temp_z->zero();

    for (int A=0; A<nS; A++){
        for (int B=0; B<nS; B++){
            for (int C=0; C<nS; C++){
                temp_x->pointer()[A][B] += x[offset_dre_plasmon_x+A*nS+C]*Dip_x->pointer()[C][B];
                temp_y->pointer()[A][B] += x[offset_dre_plasmon_y+A*nS+C]*Dip_y->pointer()[C][B];
                temp_z->pointer()[A][B] += x[offset_dre_plasmon_z+A*nS+C]*Dip_z->pointer()[C][B];
            }
        }
    }

    for (int A=0; A<nS; A++){
        dipole_p_x += temp_x->pointer()[A][A];
        dipole_p_y += temp_y->pointer()[A][A];
        dipole_p_z += temp_z->pointer()[A][A];
    } 

    // kout = f( t( n + mh ) , y( n ) + mh kin)

    ExtField(t); 

    // electronic part
    ElectronicContribution(&x[offset_dre],&x[offset_dim],&dxdt[offset_dre],&dxdt[offset_dim]);


    // 3 components of plasmon part
    PlasmonContribution(&x[offset_dre_plasmon_x],&x[offset_dim_plasmon_x],&dxdt[offset_dre_plasmon_x],&dxdt[offset_dim_plasmon_x],Dip_x,Hp_x,polarization[0]);
    PlasmonContribution(&x[offset_dre_plasmon_x],&x[offset_dim_plasmon_y],&dxdt[offset_dre_plasmon_x],&dxdt[offset_dim_plasmon_x],Dip_y,Hp_y,polarization[1]);
    PlasmonContribution(&x[offset_dre_plasmon_x],&x[offset_dim_plasmon_z],&dxdt[offset_dre_plasmon_x],&dxdt[offset_dim_plasmon_x],Dip_z,Hp_z,polarization[2]);

    // 3 components of interaction term
    InteractionContribution(&x[offset_dre],&x[offset_dim],&dxdt[offset_dre],&dxdt[offset_dim],
                            &x[offset_dre_plasmon_x],&x[offset_dim_plasmon_x],&dxdt[offset_dre_plasmon_x],&dxdt[offset_dim_plasmon_x],
                            Hp_int_x,
                            dipole_pot_x,
                            e_dip_x,
                            dipole_p_x);

    InteractionContribution(&x[offset_dre],&x[offset_dim],&dxdt[offset_dre],&dxdt[offset_dim],
                            &x[offset_dre_plasmon_y],&x[offset_dim_plasmon_y],&dxdt[offset_dre_plasmon_y],&dxdt[offset_dim_plasmon_y],
                            Hp_int_y,
                            dipole_pot_y,
                            e_dip_y,
                            dipole_p_y);

    InteractionContribution(&x[offset_dre],&x[offset_dim],&dxdt[offset_dre],&dxdt[offset_dim],
                            &x[offset_dre_plasmon_z],&x[offset_dim_plasmon_z],&dxdt[offset_dre_plasmon_z],&dxdt[offset_dim_plasmon_z],
                            Hp_int_z,
                            dipole_pot_z,
                            e_dip_z,
                            dipole_p_z);

    BuildLindblad(&x[offset_dre_plasmon_x],&x[offset_dim_plasmon_x],&dxdt[offset_dre_plasmon_x],&dxdt[offset_dim_plasmon_x]);
    BuildLindblad(&x[offset_dre_plasmon_y],&x[offset_dim_plasmon_y],&dxdt[offset_dre_plasmon_y],&dxdt[offset_dim_plasmon_y]);
    BuildLindblad(&x[offset_dre_plasmon_z],&x[offset_dim_plasmon_z],&dxdt[offset_dre_plasmon_z],&dxdt[offset_dim_plasmon_z]);

}

void TDHF::FFTW(){

    fftw_plan p;
    for (long int i=0; i<extra_pts; i++){
        corr_func[fftw_iter+i][0] = 0.0;
        corr_func[fftw_iter+i][1] = 0.0;
    }
    p = fftw_plan_dft_1d((int)(extra_pts+fftw_iter),corr_func,corr_func,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    
    fftw_plan p2;
    for (long int i=0; i<extra_pts; i++){
        corr_func2[fftw_iter+i][0] = 0.0;
        corr_func2[fftw_iter+i][1] = 0.0;
    }
    p2 = fftw_plan_dft_1d((int)(extra_pts+fftw_iter),corr_func2,corr_func2,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p2);
    fftw_destroy_plan(p2);
    
    fftw_plan p3;
    for (long int i=0; i<extra_pts; i++){
        corr_func3[fftw_iter+i][0] = 0.0;
        corr_func3[fftw_iter+i][1] = 0.0;
    }
    p3 = fftw_plan_dft_1d((int)(extra_pts+fftw_iter),corr_func3,corr_func3,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p3);
    fftw_destroy_plan(p3);
    
    fftw_plan p4;
    for (long int i=0; i<extra_pts; i++){
        corr_func4[fftw_iter+i][0] = 0.0;
        corr_func4[fftw_iter+i][1] = 0.0;
    }
    p4 = fftw_plan_dft_1d((int)(extra_pts+fftw_iter),corr_func4,corr_func4,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p4);
    fftw_destroy_plan(p4);
    
    fftw_plan p5;
    for (long int i=0; i<extra_pts; i++){
        corr_func5[fftw_iter+i][0] = 0.0;
        corr_func5[fftw_iter+i][1] = 0.0;
    }
    p5 = fftw_plan_dft_1d((int)(extra_pts+fftw_iter),corr_func5,corr_func5,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p5);
    fftw_destroy_plan(p5);
    
    fftw_plan p6;
    for (long int i=0; i<extra_pts; i++){
        corr_func6[fftw_iter+i][0] = 0.0;
        corr_func6[fftw_iter+i][1] = 0.0;
    }
    p6 = fftw_plan_dft_1d((int)(extra_pts+fftw_iter),corr_func6,corr_func6,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p6);
    fftw_destroy_plan(p6);
    
    fftw_plan p7;
    for (long int i=0; i<extra_pts; i++){
        corr_func7[fftw_iter+i][0] = 0.0;
        corr_func7[fftw_iter+i][1] = 0.0;
    }
    p7 = fftw_plan_dft_1d((int)(extra_pts+fftw_iter),corr_func7,corr_func7,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p7);
    fftw_destroy_plan(p7);
    
    fftw_plan p8;
    for (long int i=0; i<extra_pts; i++){
        corr_func8[fftw_iter+i][0] = 0.0;
        corr_func8[fftw_iter+i][1] = 0.0;
    }
    p8 = fftw_plan_dft_1d((int)(extra_pts+fftw_iter),corr_func8,corr_func8,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p8);
    fftw_destroy_plan(p8);
    
    fftw_plan p9;
    for (long int i=0; i<extra_pts; i++){
        corr_func9[fftw_iter+i][0] = 0.0;
        corr_func9[fftw_iter+i][1] = 0.0;
    }
    p9 = fftw_plan_dft_1d((int)(extra_pts+fftw_iter),corr_func9,corr_func9,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p9);
    fftw_destroy_plan(p9);

}

// output absorption spectrum
void TDHF::Spectrum(){
  int i;
  double val,valr,vali,twopi = 2.0*M_PI;
  double w;
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


  // fourier transform:
  int nfreq = 5001;
  double maxfreq = 30*0.08188379587298;
  double df = maxfreq / (nfreq - 1);
  double eps_0 = 1.0 / (4.0 * M_PI);
  eps_med = options_.get_double("EPSILON_M");
  //double eps_med = eps_0; //2.25;
  for (int i=1; i<(int)(fftw_iter+extra_pts); i++){
      w    = twopi*i/((extra_pts+fftw_iter)*time_step);
      if (w*pc_hartree2ev>max_freq) break;
      valr = corr_func[i][0]/fftw_iter;
      vali = corr_func[i][1]/fftw_iter;
      val  = sqrt(valr*valr + vali*vali);


      double e_dip_x_r   = corr_func[i][0]/fftw_iter;
      double e_dip_x_i   = corr_func[i][1]/fftw_iter;
      double e_dip_y_r   = corr_func2[i][0]/fftw_iter;
      double e_dip_y_i   = corr_func2[i][1]/fftw_iter;
      double e_dip_z_r   = corr_func3[i][0]/fftw_iter;
      double e_dip_z_i   = corr_func3[i][1]/fftw_iter;

      double p_dip_x_r   = corr_func4[i][0]/fftw_iter;
      double p_dip_x_i   = corr_func4[i][1]/fftw_iter;
      double p_dip_y_r   = corr_func5[i][0]/fftw_iter;
      double p_dip_y_i   = corr_func5[i][1]/fftw_iter;
      double p_dip_z_r   = corr_func6[i][0]/fftw_iter;
      double p_dip_z_i   = corr_func6[i][1]/fftw_iter;

      double dx_r = e_dip_x_r + p_dip_x_r;
      double dx_i = e_dip_x_i + p_dip_x_i;

      double dy_r = e_dip_y_r + p_dip_y_r;
      double dy_i = e_dip_y_i + p_dip_y_i;

      double dz_r = e_dip_z_r + p_dip_z_r;
      double dz_i = e_dip_z_i + p_dip_z_i;

      double ave_d_r = (dx_r + dy_r + dz_r)/3.0;
      double ave_d_i = (dx_i + dy_i + dz_i)/3.0;

      double field_x_r   = corr_func7[i][0]/fftw_iter;
      double field_x_i   = corr_func7[i][1]/fftw_iter;
      double field_y_r   = corr_func8[i][0]/fftw_iter;
      double field_y_i   = corr_func8[i][1]/fftw_iter;
      double field_z_r   = corr_func9[i][0]/fftw_iter;
      double field_z_i   = corr_func9[i][1]/fftw_iter;

      double field_x_mag = sqrt(field_x_r*field_x_r + field_x_i*field_x_i);
      double field_y_mag = sqrt(field_y_r*field_y_r + field_y_i*field_y_i);
      double field_z_mag = sqrt(field_z_r*field_z_r + field_z_i*field_z_i);

      double ave_f_r = (field_x_r + field_y_r + field_z_r)/3.0;
      double ave_f_i = (field_x_i + field_y_i + field_z_i)/3.0;

      double total_field_sq = field_x_r*field_x_r
                            + field_y_r*field_y_r
                            + field_z_r*field_z_r
                            + field_x_i*field_x_i
                            + field_y_i*field_y_i
                            + field_z_i*field_z_i;

      double total_field_sq_x = field_x_r*field_x_r + field_x_i*field_x_i;
      double total_field_sq_y = field_y_r*field_y_r + field_y_i*field_y_i;
      double total_field_sq_z = field_z_r*field_z_r + field_z_i*field_z_i;


//printf("%20.12le\n",total_field_sq);

      // From Gray's paper: Phys. Rev. B 88, 075411 (2013)

      
      double sca_x = pow((w / 137.03),4) * (dx_i * dx_i + dx_r * dx_r) /(6.0*M_PI*eps_0*eps_0*(total_field_sq_x));
      double sca_y = pow((w / 137.03),4) * (dy_i * dy_i + dy_r * dy_r) /(6.0*M_PI*eps_0*eps_0*(total_field_sq_y));
      double sca_z = pow((w / 137.03),4) * (dz_i * dz_i + dz_r * dz_r) /(6.0*M_PI*eps_0*eps_0*(total_field_sq_z));

      double abs_x = w*(dx_i*field_x_r - dx_r*field_x_i)/(total_field_sq_x*sqrt(eps_med)*eps_0*137.03);
      double abs_y = w*(dy_i*field_y_r - dy_r*field_y_i)/(total_field_sq_y*sqrt(eps_med)*eps_0*137.03);
      double abs_z = w*(dz_i*field_z_r - dz_r*field_z_i)/(total_field_sq_z*sqrt(eps_med)*eps_0*137.03);

      double ext_x = sca_x + abs_x;
      double ext_y = sca_y + abs_y;
      double ext_z = sca_z + abs_z;

      double sca = (sca_x + sca_y + sca_z);
      double abs = (abs_x + abs_y + abs_z);
      double ext = (ext_x + ext_y + ext_z);

      //double sca = pow((w / 137.03),4) * (ave_d_i * ave_d_i + ave_d_r * ave_d_r) / (6.0*M_PI*eps_0*(total_field_sq));
      //double abs = w*(ave_d_i*ave_f_r - ave_d_r*ave_f_i)/(total_field_sq*sqrt(eps_med)*eps_0*137.03);
      //double ext = sca + abs; 

       

      // Taken from http://arxiv.org/abs/1312.0899v1

      //double sca = pow(w,4)*(dz_i * dz_i + dz_r * dz_r)/(6.0*M_PI*sqrt(eps_0)*137.03*(total_field_sq_z));

      //double ext = -w*(-dz_i*field_z_r + dz_r*field_z_i)/(total_field_sq_z*sqrt(eps_0));

      //double abs = ext - sca;

      //fprintf(outfile,"      @Frequency %20.12lf %20.12le %20.12le %20.12le \n",w*pc_hartree2ev,sca,abs,ext);

      //fprintf(outfile,"      @Frequency %20.12lf %20.12lf %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le \n",w*pc_hartree2ev,sca_cross,mur,mui,mumr,mumi,mupr,mupi,er,ei);
      fprintf(outfile,"      @Frequency %20.12lf %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le %20.12le \n",w*pc_hartree2ev,sca,abs,ext,sca_x,abs_x,ext_x,sca_y,abs_y,ext_y,sca_z,abs_z,ext_z,field_x_mag,field_y_mag,field_z_mag);
  }
}

void TDHF::TransformIntegralsFull() {

    long int o = ndoccact;
    long int v = nvirt;
    long int full = o+v+nfzc+nfzv;

    nQ = (int)Process::environment.globals["NAUX (CC)"];

    long int memory = Process::environment.get_memory();
    // subtract out 20 full*full + 250 MB to be sure we have enough room
    memory -= sizeof(double)* 20L * full * full - 250L * 1024L * 1024L;

    if ( memory < sizeof(double) * (2L*full*full*nQ) ) {
        throw PsiException("TDHF::TransformIntegrals: not enough memory",__FILE__,__LINE__);
    }

    double * myQmo = (double*)malloc(full*full*nQ*sizeof(double));
    memset((void*)myQmo,'\0',full*full*nQ*sizeof(double));

    double ** Ca = Ca_->pointer();

    // available memory:
    memory -= sizeof(double) * (2L*full*full*nQ);
    int ndoubles = memory / sizeof(double) / 2;
    if ( nso * nso * nQ < ndoubles ) ndoubles = nso*nso*nQ;

    double * buf1 = (double*)malloc(ndoubles * sizeof(double));
    double * buf2 = (double*)malloc(ndoubles * sizeof(double));
    memset((void*)buf1,'\0',ndoubles*sizeof(double));
    memset((void*)buf2,'\0',ndoubles*sizeof(double));

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
    long int * rowdims = (long int *)malloc(nrows*sizeof(long int));
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qso CC",(char*)&buf1[0],rowdims[row]*nso*nso*sizeof(double),addr1,&addr1);
        F_DGEMM('n','n',full,nso*rowdims[row],nso,1.0,&Ca[0][0],full,buf1,nso,0.0,buf2,full);
        for (int q = 0; q < rowdims[row]; q++) {
            for (int mu = 0; mu < nso; mu++) {
                C_DCOPY(full,buf2+q*nso*full+mu*full,1,buf1+q*nso*full+mu,nso);
            }
        }
        F_DGEMM('n','n',full,full*rowdims[row],nso,1.0,&Ca[0][0],full,buf1,nso,0.0,buf2,full);

        // Qmo
        #pragma omp parallel for schedule (static)
        for (int q = 0; q < rowdims[row]; q++) {
            for (int a = 0; a < full; a++) {
                for (int b = 0; b < full; b++) {
                    myQmo[(q+rowdims[0]*row)*full*full+a*full+b] = buf2[q*full*full+a*full+b];
                }
            }
        }
    }
    //add frozen-core contribution to oeis
    for (int i = nfzc; i < full; i++) {
        for (int j = nfzc; j < full; j++) {
            double dum = 0.0;
            for (int k = 0; k < nfzc; k++) {
                for (int q = 0; q < nQ; q++) {
                    dum += 2.0 * myQmo[q*full*full+i*full+j] * myQmo[q*full*full+k*full+k];
                    dum -=       myQmo[q*full*full+i*full+k] * myQmo[q*full*full+k*full+j];
                }
            }
            T->pointer()[i][j] += dum;
        }
    }

    free(myQmo);

    free(buf1);
    free(buf2);
    free(rowdims);
    psio->close(PSIF_DCC_QSO,1);
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

    Qmo = (double*)malloc((nmo+nfzc)*(nmo+nfzc)*nQ*sizeof(double));
    memset((void*)Qmo,'\0',(nmo+nfzc)*(nmo+nfzc)*nQ*sizeof(double));

    Ire = (double*)malloc(nmo*nmo*nQ*sizeof(double));
    Iim = (double*)malloc( ( nmo * nmo > nQ ? nmo * nmo : nQ ) *sizeof(double));
    memset((void*)Ire,'\0',nmo*nmo*nQ*sizeof(double));
    memset((void*)Iim,'\0',(nmo*nmo>nQ ? nmo*nmo : nQ)*sizeof(double));

    double ** Ca = Ca_->pointer();

    // available memory:
    memory -= sizeof(double) * (2L*nmo*nmo*nQ);
    int ndoubles = memory / sizeof(double) / 2;
    if ( nso * nso * nQ < ndoubles ) ndoubles = nso*nso*nQ;

    double * buf1 = (double*)malloc(ndoubles * sizeof(double));
    double * buf2 = (double*)malloc(ndoubles * sizeof(double));
    memset((void*)buf1,'\0',ndoubles*sizeof(double));
    memset((void*)buf2,'\0',ndoubles*sizeof(double));

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
    long int * rowdims = (long int*)malloc(nrows*sizeof(long int));
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qso CC",(char*)&buf1[0],rowdims[row]*nso*nso*sizeof(double),addr1,&addr1);
        F_DGEMM('n','n',full,nso*rowdims[row],nso,1.0,&Ca[0][0],full,buf1,nso,0.0,buf2,full);
        for (int q = 0; q < rowdims[row]; q++) {
            for (int mu = 0; mu < nso; mu++) {
                C_DCOPY(full,buf2+q*nso*full+mu*full,1,buf1+q*nso*full+mu,nso);
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
    free(rowdims);
    psio->close(PSIF_DCC_QSO,1);
}

void TDHF::ElectronicContribution(double* tempr,double* tempi,double* kre,double* kim) {

    Fre->zero();
    Fim->zero();

    // F1 = J1 - K1 + h1 - mu.E
    BuildFock(tempr,tempi,true);
    F1re->copy(Fre_temp);
    F1im->copy(Fim_temp);

    F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(F1im->pointer()[0][0]),nmo,tempr,nmo,0.0,kre,nmo);
    F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(F1re->pointer()[0][0]),nmo,tempi,nmo,1.0,kre,nmo);
    F_DGEMM('n','n',nmo,nmo,nmo,-1.0,tempr,nmo,&(F1im->pointer()[0][0]),nmo,1.0,kre,nmo);
    F_DGEMM('n','n',nmo,nmo,nmo,-1.0,tempi,nmo,&(F1re->pointer()[0][0]),nmo,1.0,kre,nmo);

    F_DGEMM('n','n',nmo,nmo,nmo,-1.0,&(F1re->pointer()[0][0]),nmo,tempr,nmo,0.0,kim,nmo);
    F_DGEMM('n','n',nmo,nmo,nmo,1.0,&(F1im->pointer()[0][0]),nmo,tempi,nmo,1.0,kim,nmo);
    F_DGEMM('n','n',nmo,nmo,nmo,1.0,tempr,nmo,&(F1re->pointer()[0][0]),nmo,1.0,kim,nmo);
    F_DGEMM('n','n',nmo,nmo,nmo,-1.0,tempi,nmo,&(F1im->pointer()[0][0]),nmo,1.0,kim,nmo);

}

}}
