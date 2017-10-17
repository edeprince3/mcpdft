/*
 * @BEGIN LICENSE
 *
 * mydft by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <tuple>
#include "psi4/libpsi4util/libpsi4util.h"

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"

// for grid
#include "psi4/libfock/points.h"
#include "psi4/libfock/cubature.h"


#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libqt/qt.h"

// jk object
#include "psi4/libfock/jk.h"

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"

// diis solver
#include "diis.h"

#include "psi4/libscf_solver/hf.h"

#include "dft.h"

namespace psi{ namespace mydft {

// the DFTSolver class derives from the Wavefunction class and inherits its members
DFTSolver::DFTSolver(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options_):
    Wavefunction(options_){
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

DFTSolver::~DFTSolver() {
}

// initialize members of the DFTSolver class
void DFTSolver::common_init() {

    shallow_copy(reference_wavefunction_);

    // number of alpha electrons
    nalpha_   = reference_wavefunction_->nalpha();

    // number of beta electrons
    nbeta_    = reference_wavefunction_->nbeta();

    // number of alpha electrons per irrep
    nalphapi_ = reference_wavefunction_->nalphapi();

    // number of beta electrons per irrep
    nbetapi_  = reference_wavefunction_->nbetapi();

    // number of doubly occupied orbitals per irrep
    doccpi_   = reference_wavefunction_->doccpi();

    // number of singly occupied orbitals per irrep
    soccpi_   = reference_wavefunction_->soccpi();

    // number of frozen core orbitals per irrep
    frzcpi_   = reference_wavefunction_->frzcpi();

    // number of frozen virtual orbitals per irrep
    frzvpi_   = reference_wavefunction_->frzvpi();

    // number of molecular orbials per irrep
    nmopi_    = reference_wavefunction_->nmopi();

    // number of symmetry orbials per irrep
    nsopi_    = reference_wavefunction_->nsopi();

    // number of irreducible representations
    nirrep_   = reference_wavefunction_->nirrep();

    // make sure we are running in c1 symmetry
    if ( nirrep_ > 1 ) {
        throw PsiException("plugin mydft only works with symmetry c1",__FILE__,__LINE__);
    }

    // total number of symmetry orbitals
    nso_      = reference_wavefunction_->nso();

    // total number of molecular orbitals
    nmo_      = reference_wavefunction_->nmo();

    // grab the molecule from the reference wave function
    molecule_ = reference_wavefunction_->molecule();

    // grab the nuclear repulsion energy from the molecule
    enuc_     = molecule_->nuclear_repulsion_energy();

    // SO/MO transformation matrices
    Ca_ = std::shared_ptr<Matrix>(reference_wavefunction_->Ca());
    Cb_ = std::shared_ptr<Matrix>(reference_wavefunction_->Cb());

    // overlap integrals
    S_  = std::shared_ptr<Matrix>(reference_wavefunction_->S());

    // SO-basis Fock matrices
    Fa_ = std::shared_ptr<Matrix>(reference_wavefunction_->Fa());
    Fb_ = std::shared_ptr<Matrix>(reference_wavefunction_->Fb());

    // SO-basis density matrices
    Da_ = std::shared_ptr<Matrix>(reference_wavefunction_->Da());
    Db_ = std::shared_ptr<Matrix>(reference_wavefunction_->Db());

    // SO-basis xc potential matrices
    Va_ = std::shared_ptr<Matrix>(new Matrix(nso_,nso_));
    Vb_ = std::shared_ptr<Matrix>(new Matrix(nso_,nso_));

    // orbital energies
    epsilon_a_= std::shared_ptr<Vector>(new Vector(nirrep_, nmopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_= std::shared_ptr<Vector>(new Vector(nirrep_, nmopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());

    // memory is from process::environment
    memory_ = Process::environment.get_memory();

    // set the wavefunction name
    name_ = "DFT";

    // unrestricted kohn sham
    same_a_b_orbs_ = false;
    same_a_b_dens_ = false;


    // allocate memory for eigenvectors and eigenvalues of the overlap matrix
    std::shared_ptr<Matrix> Sevec ( new Matrix(nso_,nso_) );
    std::shared_ptr<Vector> Seval ( new Vector(nso_) );

    // build S^(-1/2) symmetric orthogonalization matrix
    S_->diagonalize(Sevec,Seval);

    Shalf_ = (std::shared_ptr<Matrix>)( new Matrix(nso_,nso_) );
    Shalf2 = (std::shared_ptr<Matrix>)( new Matrix(nso_,nso_) );
    for (int mu = 0; mu < nso_; mu++) {
        Shalf_->pointer()[mu][mu] = 1.0 / sqrt(Seval->pointer()[mu]);
        Shalf2->pointer()[mu][mu] = sqrt(Seval->pointer()[mu]);
    }

    // transform Seval back to nonorthogonal basis
    Shalf_->back_transform(Sevec);
    Shalf2->back_transform(Sevec);


    // obtain phi(r), del phi(r)

    // evaluate basis function values on a grid:

    // the only way I can figure out to get a properly initialized grid is to 
    // read in a reference wave function from a previous dft job
    scf::HF* scfwfn = (scf::HF*)reference_wavefunction_.get();
    std::shared_ptr<SuperFunctional> functional = scfwfn->functional();
    std::shared_ptr<VBase> potential = scfwfn->V_potential();

    // phi matrix (real-space <- AO basis mapping)
    // since grid is stored in blocks, we need to build a full phi matrix
    // from the individual blocks:
    std::shared_ptr<PointFunctions> points_func = potential->properties()[0];
    points_func->set_pointers(Da_,Db_);

    int nblocks = potential->nblocks();
    phi_points_       = 0;
    int max_functions = 0;
    int max_points    = 0;
    for (int myblock = 0; myblock < nblocks; myblock++) {
        std::shared_ptr<BlockOPoints> block = potential->get_block(myblock);
        points_func->compute_points(block);
        int npoints = block->npoints();
        phi_points_ += npoints;
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();
        if ( nlocal > max_functions ) max_functions = nlocal;
        if ( npoints > max_points )   max_points    = npoints;
    }

    super_phi_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI",phi_points_,nso_));
    temp_phi_  = std::shared_ptr<Matrix>(new Matrix("TEMP PHI",phi_points_,nso_));
    grid_x_    = std::shared_ptr<Vector>(new Vector("GRID X",phi_points_));
    grid_y_    = std::shared_ptr<Vector>(new Vector("GRID Y",phi_points_));
    grid_z_    = std::shared_ptr<Vector>(new Vector("GRID Z",phi_points_));
    grid_w_    = std::shared_ptr<Vector>(new Vector("GRID W",phi_points_));

    phi_points_ = 0;
    for (int myblock = 0; myblock < nblocks; myblock++) {
        std::shared_ptr<BlockOPoints> block = potential->get_block(myblock);
        points_func->compute_points(block);
        int npoints = block->npoints();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        double ** phi = points_func->basis_value("PHI")->pointer();

        double * x = block->x();
        double * y = block->y();
        double * z = block->z();
        double * w = block->w();

        for (int p = 0; p < npoints; p++) {

            grid_x_->pointer()[phi_points_ + p] = x[p];
            grid_y_->pointer()[phi_points_ + p] = y[p];
            grid_z_->pointer()[phi_points_ + p] = z[p];
            grid_w_->pointer()[phi_points_ + p] = w[p];

            for (int nu = 0; nu < nlocal; nu++) {
                int nug = function_map[nu];
                super_phi_->pointer()[phi_points_ + p][nug] = phi[p][nu];
            }
        }
        phi_points_ += npoints;
    }

    // grab AO->SO transformation matrix
    std::shared_ptr<Matrix> ao2so = reference_wavefunction_->aotoso();

    // transform one index of super phi matrix (AO->SO)
    std::shared_ptr<Matrix> temp (new Matrix(super_phi_));

    for (int p = 0; p < phi_points_; p++) {
        for (int sigma = 0; sigma < nso_; sigma++) {
            double dum = 0.0;
            for (int nu = 0; nu < nso_; nu++) {
                dum += super_phi_->pointer()[p][nu] * ao2so->pointer()[nu][sigma];
            }
            temp->pointer()[p][sigma] = dum;
        }
    }

    // AED: don't think we need this now.  use this to do so-mo transformed phi later
    C_DCOPY(nso_*phi_points_,temp->pointer()[0],1,super_phi_->pointer()[0],1);
/*
    // 
    // p(r) = phi . D' . phi^T, with D' in the SO basis
    // 
    // but, we're working with D in an orthonormal basis
    // 
    // D = S^{1/2}.D'.S^{1/2} 
    // 
    // p(r) = (phi.S^{-1/2}).D.(phi.S^{-1/2})^T
    // 
    // so, we need to tack on an additional S^{-1/2}
    // to the phi matrix
    // 
    for (int p = 0; p < phi_points_; p++) {
        for (int sigma = 0; sigma < nso_; sigma++) {
            double dum = 0.0;
            for (int nu = 0; nu < nso_; nu++) {
                dum += temp->pointer()[p][nu] * Shalf_->pointer()[nu][sigma];
            }
            super_phi_->pointer()[p][sigma] = dum;
        }
    }
*/




}

double DFTSolver::compute_energy() {

    // grab the one-electron integrals from MintsHelper:
    std::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));

    // one-electron kinetic energy integrals
    std::shared_ptr<Matrix> T = mints->so_kinetic();

    // one-electron potential energy integrals
    std::shared_ptr<Matrix> V = mints->so_potential();

    // build the core hamiltonian
    std::shared_ptr<Matrix> h = (std::shared_ptr<Matrix>)(new Matrix(T));
    h->add(V);

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    // get auxiliary basis:
    std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

    // total number of auxiliary basis functions
    int nQ = auxiliary->nbf();

    // determine the DFT functional and initialize the potential object
    scf::HF* scfwfn = (scf::HF*)reference_wavefunction_.get();
    std::shared_ptr<SuperFunctional> functional = scfwfn->functional();
    std::shared_ptr<VBase> potential = VBase::build_V(primary,functional,options_,(options_.get_str("REFERENCE") == "RKS" ? "RV" : "UV"));

    potential->initialize();

    // print the ks information
    potential->print_header();

    // JK object
    std::shared_ptr<DFJK> jk = (std::shared_ptr<DFJK>)(new DFJK(primary,auxiliary));

    // memory for jk (say, 80% of what is available)
    jk->set_memory(0.8 * Process::environment.get_memory());

    // integral cutoff
    jk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

    // Do J/K/wK?
    bool is_x_lrc  = functional->is_x_lrc();
    if ( options_["IP_FITTING"].has_changed() ) {
        if ( options_.get_bool("IP_FITTING") ) {
            is_x_lrc = true;
        }
    }
    double x_omega = functional->x_omega();
    if ( options_["DFT_OMEGA"].has_changed() ) {
        x_omega = options_.get_double("DFT_OMEGA");
    }

    jk->set_do_J(true);
    jk->set_do_K(functional->is_x_hybrid());
    jk->set_do_wK(is_x_lrc);
    jk->set_omega(x_omega);

    jk->initialize();

    // grab some input options_
    double e_convergence = options_.get_double("E_CONVERGENCE");
    double d_convergence = options_.get_double("D_CONVERGENCE");
    int maxiter          = options_.get_int("MAXITER");

    outfile->Printf("\n");
    outfile->Printf("    No. basis functions:            %5i\n",nso_);
    outfile->Printf("    No. auxiliary basis functions:  %5i\n",nQ);
    outfile->Printf("    No. alpha electrons:            %5i\n",nalpha_);
    outfile->Printf("    No. beta electrons:             %5i\n",nbeta_);
    outfile->Printf("    e_convergence:             %10.3le\n",e_convergence);
    outfile->Printf("    d_convergence:             %10.3le\n",d_convergence);
    outfile->Printf("    maxiter:                        %5i\n",maxiter);
    outfile->Printf("\n");
    outfile->Printf("\n");

    // form F' = ST^(-1/2) F S^(-1/2), where F = h
    Fa_->copy(h);
    Fb_->copy(h);

    std::shared_ptr<Matrix> Fprime_a ( new Matrix(Fa_) );
    std::shared_ptr<Matrix> Fprime_b ( new Matrix(Fb_) );

    Fprime_a->transform(Shalf_);
    Fprime_b->transform(Shalf_);

    // allocate memory for eigenvectors of F'
    std::shared_ptr<Matrix> Fevec_a ( new Matrix(nso_,nso_) );

    std::shared_ptr<Matrix> Fevec_b ( new Matrix(nso_,nso_) );

    // diagonalize F' to obtain C'
    Fprime_a->diagonalize(Fevec_a,epsilon_a_,ascending);
    Fprime_b->diagonalize(Fevec_b,epsilon_b_,ascending);

    // Find C = S^(-1/2)C'
    Ca_->gemm(false,false,1.0,Shalf_,Fevec_a,0.0);
    Cb_->gemm(false,false,1.0,Shalf_,Fevec_b,0.0);

    // Construct density from C
    C_DGEMM('n','t',nso_,nso_,nalpha_,1.0,&(Ca_->pointer()[0][0]),nso_,&(Ca_->pointer()[0][0]),nso_,0.0,&(Da_->pointer()[0][0]),nso_);
    C_DGEMM('n','t',nso_,nso_,nbeta_,1.0,&(Cb_->pointer()[0][0]),nso_,&(Cb_->pointer()[0][0]),nso_,0.0,&(Db_->pointer()[0][0]),nso_);

    // initial energy, E = D(H+F) + Enuc
    double e_current = enuc_;
    e_current       += 0.5 * Da_->vector_dot(h);
    e_current       += 0.5 * Db_->vector_dot(h);
    e_current       += 0.5 * Da_->vector_dot(Fa_);
    e_current       += 0.5 * Db_->vector_dot(Fb_);

    // SCF iterations

    double e_last    = 0.0;
    double dele      = 0.0;
    double gnorm_a   = 0.0;
    double gnorm_b   = 0.0;

    outfile->Printf("\n");
    outfile->Printf("    Guess energy:  %20.12lf\n",e_current);
    outfile->Printf("\n");
    outfile->Printf("    ==>  Begin SCF Iterations <==\n");
    outfile->Printf("\n");
    outfile->Printf("    ");
    outfile->Printf(" Iter ");
    outfile->Printf("              energy ");
    outfile->Printf("                  dE ");
    outfile->Printf("         RMS |[F,P]| ");
    outfile->Printf("\n");

    bool do_diis = options_.get_bool("DIIS");

    // Initialize our diis extrapolation manager.  Note that,
    // in this implemenation, the DIIS managar allocates 
    // memory for two buffers of the size of the Fock matrix.

    std::shared_ptr<DIIS> diis (new DIIS(nso_*nso_));

    int iter = 0;
    do {

        e_last = e_current;

        // grab occupied alpha orbitals (the first na)
        std::shared_ptr<Matrix> myCa (new Matrix(Ca_) );
        myCa->zero();

        // grab occupied beta orbitals (the first nb)
        std::shared_ptr<Matrix> myCb (new Matrix(Cb_) );
        myCb->zero();

        for (int mu = 0; mu < nso_; mu++) {
            for (int i = 0; i < nalpha_; i++) {
                myCa->pointer()[mu][i] = Ca_->pointer()[mu][i];
            }
            for (int i = 0; i < nbeta_; i++) {
                myCb->pointer()[mu][i] = Cb_->pointer()[mu][i];
            }
        }

        // push occupied orbitals onto JK object
        std::vector< std::shared_ptr<Matrix> >& C_left  = jk->C_left();
        C_left.clear();
        C_left.push_back(myCa);
        C_left.push_back(myCb);

        // form J/K
        jk->compute();

        // form Fa = h + Ja + Jb
        Fa_->copy(h);
        Fa_->add(jk->J()[0]);
        Fa_->add(jk->J()[1]);

        // form Fb = h + Ja + Jb
        Fb_->copy(h);
        Fb_->add(jk->J()[0]);
        Fb_->add(jk->J()[1]);

        // Construct density from C
        C_DGEMM('n','t',nso_,nso_,nalpha_,1.0,&(Ca_->pointer()[0][0]),nso_,&(Ca_->pointer()[0][0]),nso_,0.0,&(Da_->pointer()[0][0]),nso_);
        C_DGEMM('n','t',nso_,nso_,nbeta_, 1.0,&(Cb_->pointer()[0][0]),nso_,&(Cb_->pointer()[0][0]),nso_,0.0,&(Db_->pointer()[0][0]),nso_);

        if (functional->needs_xc()) {

            // set a/b densities in potential object
            potential->set_D({Da_, Db_});

            // evaluate a/b potentials
            potential->compute_V({Va_,Vb_});

            // form Fa/b = h + Ja + Jb + Va/b
            Fa_->add(Va_);
            Fb_->add(Vb_);

        }

        // exact exchange?
        if (functional->is_x_hybrid()) {
            // form F = h + 2*J + V - alpha K
            double alpha = functional->x_alpha();
            Fa_->axpy(-alpha,jk->K()[0]);
            Fb_->axpy(-alpha,jk->K()[1]); 
        }

        // LRC functional?
        if (is_x_lrc) {
            // form Fa/b = h + Ja + Jb + Va/b - alpha Ka/b - beta wKa/b
            double beta = 1.0 - functional->x_alpha();
            Fa_->axpy(-beta,jk->wK()[0]);
            Fb_->axpy(-beta,jk->wK()[1]);
        }

        // evaluate the current energy
        double one_electron_energy = Da_->vector_dot(h)
                                   + Db_->vector_dot(h);

        double two_electron_energy = 0.5 * Da_->vector_dot(jk->J()[0])
                                   + 0.5 * Da_->vector_dot(jk->J()[1])
                                   + 0.5 * Db_->vector_dot(jk->J()[0])
                                   + 0.5 * Db_->vector_dot(jk->J()[1]);

        if (functional->is_x_hybrid()) {
            double alpha = functional->x_alpha();
            two_electron_energy -= 0.5 * alpha * Da_->vector_dot(jk->K()[0]);
            two_electron_energy -= 0.5 * alpha * Db_->vector_dot(jk->K()[1]);
        }

        if (is_x_lrc) {
            double beta = 1.0 - functional->x_alpha();
            two_electron_energy -= 0.5 * beta * Da_->vector_dot(jk->wK()[0]);
            two_electron_energy -= 0.5 * beta * Db_->vector_dot(jk->wK()[1]);
        }

        double exchange_correlation_energy = 0.0;
        if (functional->needs_xc()) {
            exchange_correlation_energy = potential->quadrature_values()["FUNCTIONAL"];


            // try evaluating exc ourselves:
            //std::shared_ptr<Vector> rho (new Vector(phi_points_));
            double ** phi = super_phi_->pointer();
            double ** Dap = Da_->pointer();
            double ** Dbp = Db_->pointer();
            double exc = 0.0;
            for (int p = 0; p < phi_points_; p++) {
                double duma = 0.0;
                double dumb = 0.0;
                for (int sigma = 0; sigma < nso_; sigma++) {
                    for (int nu = 0; nu < nso_; nu++) {
                        duma += phi[p][sigma] * phi[p][nu] * Dap[sigma][nu];
                        dumb += phi[p][sigma] * phi[p][nu] * Dbp[sigma][nu];
                    }
                }
                //rho->pointer()[p] = duma + dumb;
                double rho_43 = pow(duma+dumb,4.0/3.0);
                exc += -0.75 * pow(3.0/M_PI,1.0/3.0) * rho_43 * grid_w_->pointer()[p];
                //exc += -9.0/8.0 * 2.0/3.0 * pow(3.0/M_PI,1.0/3.0) * rho_43 * grid_w_->pointer()[p];
            }
            printf("%20.12lf %20.12lf\n",exc,exchange_correlation_energy);
        }

        e_current  = enuc_;
        e_current += one_electron_energy;
        e_current += two_electron_energy;
        e_current += exchange_correlation_energy;

        //printf("%20.12lf %20.12lf %20.12lf\n",one_electron_energy,two_electron_energy,exchange_correlation_energy);

        // dele
        dele = e_current - e_last;

        // form F' = ST^(-1/2) F S^(-1/2)
        Fprime_a->copy(Fa_);
        Fprime_a->transform(Shalf_);

        Fprime_b->copy(Fb_);
        Fprime_b->transform(Shalf_);

        // Now, we add a few steps for the DIIS procedure.

        // The extrapolated parameter in DIIS for SCF
        // are the alpha and beta Fock matrices, Fa' and Fb', 
        // in the orthonormal basis. This DIIS manager will 
        // write the current Fa' and Fb' to disk to avoid 
        // storing multiple copies in main memory.

        if ( do_diis ) {
            diis->WriteVector(&(Fprime_a->pointer()[0][0]),&(Fprime_b->pointer()[0][0]));
        }

        // The error vector in DIIS for SCF is defined as 
        // the orbital gradient, in the orthonormal basis:
        // 
        // ea = ST^{-1/2} [FaDaS - SDaFa] S^{-1/2}
        // eb = ST^{-1/2} [FbDbS - SDbFb] S^{-1/2}

        std::shared_ptr<Matrix> grad_a = OrbitalGradient(Da_,Fa_,Shalf_);
        std::shared_ptr<Matrix> grad_b = OrbitalGradient(Db_,Fb_,Shalf_);
       
        // We will use the RMS of the orbital gradient 
        // to monitor convergence.
        gnorm_a = grad_a->rms();
        gnorm_b = grad_b->rms();

        // The DIIS manager will write the current error vector to disk.
        if ( do_diis ) {
            diis->WriteErrorVector(&(grad_a->pointer()[0][0]),&(grad_b->pointer()[0][0]));
        }

        // The DIIS manager extrapolates the Fock matrices, using
        // the Fock matrices and error vectors generated in this
        // and previous iterations.
        if ( do_diis ) {
            diis->Extrapolate(&(Fprime_a->pointer()[0][0]),&(Fprime_b->pointer()[0][0]));
        }

        // Now, we resume the usual SCF procedure, using the
        // extrapolated Fock matrices Fa' and Fb'.

        // Diagonalize F' to obtain C'
        Fprime_a->diagonalize(Fevec_a,epsilon_a_,ascending);
        Fprime_b->diagonalize(Fevec_b,epsilon_b_,ascending);

        // Find C = S^(-1/2)C'
        Ca_->gemm(false,false,1.0,Shalf_,Fevec_a,0.0);
        Cb_->gemm(false,false,1.0,Shalf_,Fevec_b,0.0);

        outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",iter,e_current,dele,0.5 * (gnorm_a + gnorm_b) ); 

        iter++;
        if ( iter > maxiter ) break;

    }while(fabs(dele) > e_convergence || 0.5 * (gnorm_a + gnorm_b) > d_convergence );

    if ( iter > maxiter ) {
        throw PsiException("Maximum number of iterations exceeded!",__FILE__,__LINE__);
    }

    outfile->Printf("\n");
    outfile->Printf("    SCF iterations converged!\n");
    outfile->Printf("\n");
    outfile->Printf("    * SCF total energy: %20.12lf\n",e_current);

    Process::environment.globals["SCF TOTAL ENERGY"] = e_current;
    Process::environment.globals["CURRENT ENERGY"]   = e_current;

    return e_current;

}

std::shared_ptr<Matrix> DFTSolver::OrbitalGradient(std::shared_ptr<Matrix> D, 
                                                   std::shared_ptr<Matrix> F, 
                                                   std::shared_ptr<Matrix> Shalf) {

    std::shared_ptr<Matrix> ShalfGradShalf(new Matrix("ST^{-1/2}(FDS - SDF)S^{-1/2}", nso_, nso_));

    std::shared_ptr<Matrix> FDSmSDF(new Matrix("FDS-SDF", nso_, nso_));
    std::shared_ptr<Matrix> DS(new Matrix("DS", nso_, nso_));

    DS->gemm(false,false,1.0,D,S_,0.0);
    FDSmSDF->gemm(false,false,1.0,F,DS,0.0);

    DS.reset();

    std::shared_ptr<Matrix> SDF(FDSmSDF->transpose());
    FDSmSDF->subtract(SDF);

    SDF.reset();

    std::shared_ptr<Matrix> ShalfGrad(new Matrix("ST^{-1/2}(FDS - SDF)", nso_, nso_));
    ShalfGrad->gemm(true,false,1.0,Shalf,FDSmSDF,0.0);
    FDSmSDF.reset();

    ShalfGradShalf->gemm(false,false,1.0,ShalfGrad,Shalf,0.0);

    ShalfGrad.reset();

    return ShalfGradShalf;
}

}} // End namespaces

