#ifndef RTTD2RDM_H
#define RTTD2RDM_H
#define PSIF_V2RDM_D2AA       270
#define PSIF_V2RDM_D2AB       271
#define PSIF_V2RDM_D2BB       272

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/vector.h"
#include "fftw3.h"

#include "psi4/libfock/jk.h"
#include "psi4/libfunctional/superfunctional.h"

using namespace psi;

namespace psi{ 

class VBase;

namespace rttd2rdm {

class TD2RDM: public Wavefunction {
public:
    TD2RDM(std::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
    ~TD2RDM();

    void common_init();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }

protected:

    /// KS Potential (the heart of the algorithm)
    std::shared_ptr<VBase> potential_;

    /// Pointer to potential's functional
    std::shared_ptr<SuperFunctional> functional_;

    /// are we doing DFT?
    bool is_dft_;

    /// RK4 
    void RK4(std::shared_ptr<Matrix> koutre, std::shared_ptr<Matrix>koutim,
        std::shared_ptr<Matrix> kinre, std::shared_ptr<Matrix>kinim, int iter, double step);

    /// external field magnitude at some time
    double ext_field_;

    /// total simulation time
    double total_time_;

    /// time step
    double time_step_;

    /// strength of external field
    double laser_amp_;

    /// frequency of external field
    double laser_freq_;

    /// transition dipole moment (for pi pulse)
    double transition_dpm_;

    /// transition dipole moment (for pi pulse)
    double laser_dpm_;

    /// pulse length
    double laser_time_;

    /// total number of time steps
    long int total_iter_;

    /// polarization of external field
    double * polarization_, *D2aa, *D2ab, *D2bb, *Da, *Db;

    /// pointer to 3-index integrals
    double ** Qp_;

    /// external field shape
    int pulse_shape_;

    /// external field pulse length
    int pulse_length_;

    /// read 2-particle reduced density matrix
    void ReadTPDM();

    /// add external field term to Fock matrix
    void ExtField(double curtime);

    /// build Fock matrix
    void BuildFock(std::shared_ptr<Matrix> Dre,std::shared_ptr<Matrix> Dim,double curtime);

    /// diagonalize complex hermitian matrix
    void DiagonalizeHermitianMatrix(long int N,double*re,double*im,double*W);

    /// compute 2-cumulant
    void Compute2Cumulant(std::shared_ptr<Matrix> Fre,std::shared_ptr<Matrix> Fim,double * Dre,double * Dim);

    /// evaluate current energy
    double CurrentEnergy(std::shared_ptr<Matrix> Dre,std::shared_ptr<Matrix> Dim);

    /// build JK matrices
    void BuildJK(std::shared_ptr<Matrix> Dre,std::shared_ptr<Matrix> Dim);

    /// build dft potential matrices
    void BuildV(std::shared_ptr<Matrix> Dre,std::shared_ptr<Matrix> Dim);

    /// generate 1-electron integrals 
    void Integrals();

    /// JK object
    std::shared_ptr<DFJK> jk_;

    /// reference energy
    double escf_;

    /// number of frozen virtual orbitals
    long int nfrzv_;

    /// number of virtual orbitals
    long int nvirt_;

    /// number of doubly occupied orbitals
    long int ndocc_;

    /// number of active doubly occupied orbitals
    long int ndoccact_;

    /// number of auxiliary basis functions
    long int nQ_;

    /// basis sets
    std::shared_ptr<BasisSet> primary_;
    std::shared_ptr<BasisSet> auxiliary_;

    /// density
    std::shared_ptr<Matrix> Dre_;
    std::shared_ptr<Matrix> Dim_;

    /// two-body cumulant
    std::shared_ptr<Matrix> C2re_;
    std::shared_ptr<Matrix> C2im_;

    /// coulomb matrix
    std::shared_ptr<Matrix> Jre_;
    std::shared_ptr<Matrix> Jim_;

    /// dft potential matrix
    std::shared_ptr<Matrix> Vre_;
    std::shared_ptr<Matrix> Vim_;

    /// exchange matrix
    std::shared_ptr<Matrix> Kre_;
    std::shared_ptr<Matrix> Kim_;

    /// exchange matrix
    std::shared_ptr<Matrix> wKre_;
    std::shared_ptr<Matrix> wKim_;

    /// rk4 buffers
    std::shared_ptr<Matrix> k1re_;
    std::shared_ptr<Matrix> k1im_;
    std::shared_ptr<Matrix> k2re_;
    std::shared_ptr<Matrix> k2im_;
    std::shared_ptr<Matrix> k3re_;
    std::shared_ptr<Matrix> k3im_;
    std::shared_ptr<Matrix> k4re_;
    std::shared_ptr<Matrix> k4im_;
    std::shared_ptr<Matrix> kre_;
    std::shared_ptr<Matrix> kim_;

    /// temporary buffers
    std::shared_ptr<Matrix> tempre_;
    std::shared_ptr<Matrix> tempim_;

    std::shared_ptr<Matrix> mux_;
    std::shared_ptr<Matrix> muy_;
    std::shared_ptr<Matrix> muz_;

    /// kinetic energy integrals
    std::shared_ptr<Matrix> T_;

    /// electron-nucleus potential energy integrals
    std::shared_ptr<Matrix> V_;

    /// total Fock matrix
    std::shared_ptr<Matrix> Fre_;
    std::shared_ptr<Matrix> Fim_;

    /// core Hamiltonian
    std::shared_ptr<Matrix> Hcore_;

    /// orthogonalization matrix
    std::shared_ptr<Matrix> Shalf_;

    /// three-index integrals in the mo basis
    std::shared_ptr<Matrix> Qmo_;

    /// buffers for FFTW
    fftw_complex * corr_func_x;
    fftw_complex * corr_func_y;
    fftw_complex * corr_func_z;

    /// fftw
    void FFTW();

    /// frequency-domain spectrum
    void Spectrum();
   
};

}}

#endif
