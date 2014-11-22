#ifndef TDHF_H
#define TDHF_H
#include<libmints/wavefunction.h>
#include<libmints/vector.h>
#include"fftw3.h"

// boost numerical integrators live here:
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

namespace boost {
template<class T> class shared_ptr;
}
typedef std::vector< double > state_type;

typedef runge_kutta4< state_type , double ,
                      state_type , double
                       > stepper_type;


namespace psi{ namespace tdhf_cqed {

class TDHF: public Wavefunction {
public:
    TDHF(boost::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
    ~TDHF();

    void common_init();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }
    //void operator()( state_type &x , state_type &dxdt , double t );
    void rk4_call_gah( state_type &x , state_type &dxdt , double t );
    //void rk4_call( state_type &x , state_type &dxdt , double t );
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

    // RK4 
    void RK4(boost::shared_ptr<Matrix> koutre, boost::shared_ptr<Matrix>koutim,
             boost::shared_ptr<Matrix> kinre, boost::shared_ptr<Matrix>kinim, int iter, double step);

    // boost rk4:
    double * rk4_buffer;
};


// horrible global class
extern boost::shared_ptr<TDHF> MyTDHF;

}}

#endif
