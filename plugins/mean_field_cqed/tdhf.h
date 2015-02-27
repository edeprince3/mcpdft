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
                       > rk4_stepper;
typedef runge_kutta_dopri5< state_type , double ,
                      state_type , double
                       > rk5_stepper;
typedef runge_kutta_fehlberg78< state_type , double ,
                      state_type , double
                       > rk78_stepper;


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
    double total_time, time_step, laser_amp, laser_freq, laser_time, total_iter, transition_dpm, *plasmon_e, plasmon_tdm_x, plasmon_tdm_y, plasmon_tdm_z, plasmon_dr, eps_med, ext_field, dipole_p_x, dipole_p_y, dipole_p_z, mdip, pdip, e_dip_x, e_dip_y, e_dip_z, nuc_dip;
    double transition_coupling, coupling_strength;
    int natom_;
    double * com_;

    long int nQ, nS;
    double *Qmo;
    double *tei;
    double *Ire,*Iim;
    double * polarization;
    double * plasmon_tdm;
    double * plasmon_coordinates;

    int offset_dre;
    int offset_dim;
    int offset_dre_plasmon_x;
    int offset_dim_plasmon_x;
    int offset_dre_plasmon_y;
    int offset_dim_plasmon_y;
    int offset_dre_plasmon_z;
    int offset_dim_plasmon_z;

    boost::shared_ptr<Matrix> T;
    boost::shared_ptr<Matrix> V;
    boost::shared_ptr<Matrix> eri;

    boost::shared_ptr<Matrix> Dre;
    boost::shared_ptr<Matrix> Dim;

    boost::shared_ptr<Matrix> Dre_plasmon_x;
    boost::shared_ptr<Matrix> Dim_plasmon_x;

    boost::shared_ptr<Matrix> Dre_plasmon_y;
    boost::shared_ptr<Matrix> Dim_plasmon_y;

    boost::shared_ptr<Matrix> Dre_plasmon_z;
    boost::shared_ptr<Matrix> Dim_plasmon_z;

    boost::shared_ptr<Matrix> Fre;
    boost::shared_ptr<Matrix> Fim;
    boost::shared_ptr<Matrix> F1re;
    boost::shared_ptr<Matrix> F1im;
    boost::shared_ptr<Matrix> Fre_temp;
    boost::shared_ptr<Matrix> Fim_temp;

    boost::shared_ptr<Matrix> Dip_x;
    boost::shared_ptr<Matrix> Dip_y;
    boost::shared_ptr<Matrix> Dip_z;

    boost::shared_ptr<Matrix> D1_e_re;
    boost::shared_ptr<Matrix> D1_e_im;
    boost::shared_ptr<Matrix> D1_p_re;
    boost::shared_ptr<Matrix> D1_p_im;
    boost::shared_ptr<Matrix> Hp_x;
    boost::shared_ptr<Matrix> Hp_int_x;
    boost::shared_ptr<Matrix> Eigvec_x;
    boost::shared_ptr<Matrix> htemp_x;
    boost::shared_ptr<Matrix> htemp_int_x;
    boost::shared_ptr<Matrix> htemp_dip_x;
    boost::shared_ptr<Matrix> Hplasmon_total_x;
    boost::shared_ptr<Matrix> Hp_y;
    boost::shared_ptr<Matrix> Hp_int_y;
    boost::shared_ptr<Matrix> Eigvec_y;
    boost::shared_ptr<Matrix> htemp_y;
    boost::shared_ptr<Matrix> htemp_int_y;
    boost::shared_ptr<Matrix> htemp_dip_y;
    boost::shared_ptr<Matrix> Hplasmon_total_y;
    boost::shared_ptr<Matrix> Hp_z;
    boost::shared_ptr<Matrix> Hp_int_z;
    boost::shared_ptr<Matrix> Eigvec_z;
    boost::shared_ptr<Matrix> htemp_z;
    boost::shared_ptr<Matrix> htemp_int_z;
    boost::shared_ptr<Matrix> htemp_dip_z;
    boost::shared_ptr<Matrix> Hplasmon_total_z;
    std::vector<boost::shared_ptr<Matrix> > dipole;
    boost::shared_ptr<Matrix> dipole_pot_x;
    boost::shared_ptr<Matrix> dipole_pot_y;
    boost::shared_ptr<Matrix> dipole_pot_z;

    boost::shared_ptr<Matrix> temp_x;
    boost::shared_ptr<Matrix> temp_y;
    boost::shared_ptr<Matrix> temp_z;

    boost::shared_ptr<Vector> Eigval_x;
    boost::shared_ptr<Vector> Eigval_y;
    boost::shared_ptr<Vector> Eigval_z;

    void DipolePotentialIntegrals();
    void TransformIntegrals();
    void TransformIntegralsFull();
    void SoToMo(int nsotemp,int nmotemp,double**mat,double**trans);
    void BuildFock(double * Dre_temp, double * Dim_temp,bool use_oe_terms);
    void BuildFockThreeIndex(double * Dre_temp, double * Dim_temp,bool use_oe_terms);

    void ElectronicContribution(double* tempr,double* tempi,double* kre,double* kim);

    void PlasmonContribution(double * tempr,
                             double * tempi,
                             double * kre,
                             double * kim,
                             boost::shared_ptr<Matrix> dip,
                             boost::shared_ptr<Matrix> Ham,
                             double pol);

    void InteractionContribution(double * tempr,
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
                                 double pdip);


    void ExtField(double curtime);

    void BuildLindblad(double * tempr,
                       double * tempi,
                       double * kre,
                       double * kim);


    // fourier transform
    void FFTW();
    void Spectrum();
    fftw_complex*corr_func;
    fftw_complex*corr_func2;
    fftw_complex*corr_func3;
    fftw_complex*corr_func4;
    fftw_complex*corr_func5;
    fftw_complex*corr_func6;
    fftw_complex*corr_func7;
    fftw_complex*corr_func8;
    fftw_complex*corr_func9;
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

    // nuclear component of dipole moment
    double nuc_dip_x_;
    double nuc_dip_y_;
    double nuc_dip_z_;

};

}}


#endif
