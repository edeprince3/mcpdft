#ifndef SDP_H
#define SDP_H
#include<libmints/wavefunction.h>
#include<libmints/vector.h>
namespace boost {
template<class T> class shared_ptr;
}

namespace psi{ namespace dq_scf{

class SDPSolver: public Wavefunction {
public:
    SDPSolver(boost::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
    ~SDPSolver();

    void common_init();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return false; }
    virtual bool same_a_b_dens() const { return false; }

protected:
    long int memory;
    long int ndocc,nmo,nvirt,ndoccact,nfzc,nso,nfzv,dim2,nconstraints;
    double * eps,escf,energy,nrm,nrmaa,nrmab,trd2,trd2aa,trd2ab,trd1a,trd1b,nrm1a,nrm1b;
    double * V1, * V2, *D1, *Q1, *R1d, *dR1d, *R1q, *dR1q;
    double * lambda, mu, * clam;
    double * c, * cval, * cerror;

    int ** ibaa, ** baa;
    int ** ibab, ** bab;
    void Maps();
    void K2();
    void Guess();
    void BuildD1();
    void BuildQ1();
    double K2D1();
    double TraceD1();
    double eval();
    void SoToMo(int nsotemp,int nmotemp,double**mat,double**trans);

    // rotate to natural orbital basis:
    void NatOrbs(int nmo,double * D1);

    // potential for ab initio dft
    void Potential();
    // back transform D2
    void BackTransform();


};

}}


#endif
