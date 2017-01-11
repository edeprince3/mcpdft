#ifndef TDHF_H
#define TDHF_H
#include"psi4/libmints/wavefunction.h"
#include"psi4/libmints/vector.h"
#include"fftw3.h"

namespace std {
template<class T> class shared_ptr;
}

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
    TDHF(std::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
    ~TDHF();

    void common_init();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }

protected:

    void RK4( int rk4_dim, double * rk4_soln , double * rk4_out, double * rk4_in , double * rk4_temp, double curtime, double step );
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

<<<<<<< HEAD
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
=======
    std::shared_ptr<Matrix> T;
    std::shared_ptr<Matrix> V;
    std::shared_ptr<Matrix> eri;

    std::shared_ptr<Matrix> Dre;
    std::shared_ptr<Matrix> Dim;

    std::shared_ptr<Matrix> Dre_plasmon_x;
    std::shared_ptr<Matrix> Dim_plasmon_x;

    std::shared_ptr<Matrix> Dre_plasmon_y;
    std::shared_ptr<Matrix> Dim_plasmon_y;

    std::shared_ptr<Matrix> Dre_plasmon_z;
    std::shared_ptr<Matrix> Dim_plasmon_z;

    std::shared_ptr<Matrix> Fre;
    std::shared_ptr<Matrix> Fim;
    std::shared_ptr<Matrix> F1re;
    std::shared_ptr<Matrix> F1im;
    std::shared_ptr<Matrix> Fre_temp;
    std::shared_ptr<Matrix> Fim_temp;

    std::shared_ptr<Matrix> Dip_x;
    std::shared_ptr<Matrix> Dip_y;
    std::shared_ptr<Matrix> Dip_z;

    std::shared_ptr<Matrix> Hp_x;
    std::shared_ptr<Matrix> Hp_int_x;
    std::shared_ptr<Matrix> Eigvec_x;
    std::shared_ptr<Matrix> htemp_x;
    std::shared_ptr<Matrix> htemp_int_x;
    std::shared_ptr<Matrix> htemp_dip_x;
    std::shared_ptr<Matrix> Hplasmon_total_x;
    std::shared_ptr<Matrix> Hp_y;
    std::shared_ptr<Matrix> Hp_int_y;
    std::shared_ptr<Matrix> Eigvec_y;
    std::shared_ptr<Matrix> htemp_y;
    std::shared_ptr<Matrix> htemp_int_y;
    std::shared_ptr<Matrix> htemp_dip_y;
    std::shared_ptr<Matrix> Hplasmon_total_y;
    std::shared_ptr<Matrix> Hp_z;
    std::shared_ptr<Matrix> Hp_int_z;
    std::shared_ptr<Matrix> Eigvec_z;
    std::shared_ptr<Matrix> htemp_z;
    std::shared_ptr<Matrix> htemp_int_z;
    std::shared_ptr<Matrix> htemp_dip_z;
    std::shared_ptr<Matrix> Hplasmon_total_z;
    std::vector<std::shared_ptr<Matrix> > dipole;
    std::shared_ptr<Matrix> dipole_pot_x;
    std::shared_ptr<Matrix> dipole_pot_y;
    std::shared_ptr<Matrix> dipole_pot_z;

    std::shared_ptr<Matrix> temp_x;
    std::shared_ptr<Matrix> temp_y;
    std::shared_ptr<Matrix> temp_z;

    std::shared_ptr<Vector> Eigval_x;
    std::shared_ptr<Vector> Eigval_y;
    std::shared_ptr<Vector> Eigval_z;
>>>>>>> b8398429... updates mean_field_cqed to be compatible with psi4-1.0 as downloaded on 12/06/16

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
                             std::shared_ptr<Matrix> dip,
                             std::shared_ptr<Matrix> Ham,
                             double pol);

    void InteractionContribution(double * tempr,
                                 double * tempi,
                                 double * kre,
                                 double * kim,
                                 double * tempr_p,
                                 double * tempi_p,
                                 double * kre_p,
                                 double * kim_p,
                                 std::shared_ptr<Matrix> Hp_int,
                                 std::shared_ptr<Matrix> dipole_pot,
                                 double mdip,   
                                 double pdip);


<<<<<<< HEAD
=======
    void PlasmonHamiltonianTransformation(std::shared_ptr<Matrix> Ham,std::shared_ptr<Matrix>Eigvec);
    void HInteraction(double * D1);


>>>>>>> b8398429... updates mean_field_cqed to be compatible with psi4-1.0 as downloaded on 12/06/16
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
    void RK4(std::shared_ptr<Matrix> koutre, std::shared_ptr<Matrix>koutim,
             std::shared_ptr<Matrix> kinre, std::shared_ptr<Matrix>kinim, int iter, double step);

    // std rk4:
    double * rk4_buffer;

    // nuclear component of dipole moment
    double nuc_dip_x_;
    double nuc_dip_y_;
    double nuc_dip_z_;

};

}}


#endif
