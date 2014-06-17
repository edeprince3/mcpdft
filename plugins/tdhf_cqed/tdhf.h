#ifndef TDHF_H
#define TDHF_H
#include<libmints/wavefunction.h>
#include<libmints/vector.h>
#include"fftw3.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{ namespace tdhf_cqed {

class TDHF: public Wavefunction {
public:
    TDHF(boost::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
    ~TDHF();

    void common_init();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }
protected:
    double * eps;
    double escf;
    long int ndocc, nvirt, nfzc, nfzv, ndoccact, nmo, nso;
    long int memory;
    double total_time, time_step, laser_amp, laser_freq, laser_time, total_iter, transition_dpm;

    long int nQ;
    double *Qmo;
    double *Ire,*Iim;
    double * polarization;

    boost::shared_ptr<Matrix> T;
    boost::shared_ptr<Matrix> V;
    boost::shared_ptr<Matrix> eri;
    boost::shared_ptr<Matrix> Dre_init;
    boost::shared_ptr<Matrix> Dim_init;
    boost::shared_ptr<Matrix> Dre;
    boost::shared_ptr<Matrix> Dim;
    boost::shared_ptr<Matrix> Vext;
    boost::shared_ptr<Matrix> Fre;
    boost::shared_ptr<Matrix> Fim;
    boost::shared_ptr<Matrix> Fre_init;
    boost::shared_ptr<Matrix> Fim_init;
    std::vector<boost::shared_ptr<Matrix> > dipole;

    void TransformIntegrals();
    void SoToMo(int nsotemp,int nmotemp,double**mat,double**trans);
    void BuildFock(double * Dre, double * Dim, double curtime);

    void CorrelationFunction();

    // fourier transform
    void FFTW();
    void Spectrum();
    fftw_complex*corr_func;
    int midpt,extra_pts;
    double max_freq;
    int nfreq,fftw_iter;
    bool linear_response;
    double * stencil, cf_real, cf_imag, * c1;
    int pulse_shape_;

};

}}

#endif
