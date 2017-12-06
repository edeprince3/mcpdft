/*
 * @BEGIN LICENSE
 *
 * mcpdft by Psi4 Developer, a plugin to:
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

#include "psi4/libqt/qt.h"

// jk object
#include "psi4/libfock/jk.h"

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

// for potential object
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libscf_solver/hf.h"

// mcpdft 
#include "mcpdft_solver.h"

namespace psi{ namespace mcpdft {

// the MCPDFTSolver class derives from the Wavefunction class and inherits its members
MCPDFTSolver::MCPDFTSolver(std::shared_ptr<Wavefunction> reference_wavefunction,Options & options_):
    Wavefunction(options_){

    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

MCPDFTSolver::~MCPDFTSolver() {
}

// initialize members of the MCPDFTSolver class
void MCPDFTSolver::common_init() {

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
        throw PsiException("plugin mcpdft only works with symmetry c1",__FILE__,__LINE__);
    }

    // total number of symmetry orbitals
    nso_      = reference_wavefunction_->nso();

    // total number of molecular orbitals
    nmo_      = reference_wavefunction_->nmo();

    // grab the molecule from the reference wave function
    molecule_ = reference_wavefunction_->molecule();

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

    // orbital energies
    epsilon_a_= std::shared_ptr<Vector>(new Vector(nirrep_, nmopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_= std::shared_ptr<Vector>(new Vector(nirrep_, nmopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());

    // memory is from process::environment
    memory_ = Process::environment.get_memory();

    // set the wavefunction name
    name_ = "DFT";

    // restricted orbitals, unrestricted rdms
    same_a_b_orbs_ = true;
    same_a_b_dens_ = false;

    // obtain phi(r), del phi(r)

    // evaluate basis function values on a grid:

    scf::HF* scfwfn = (scf::HF*)reference_wavefunction_.get();
    std::shared_ptr<SuperFunctional> functional = scfwfn->functional();
    std::shared_ptr<VBase> potential = scfwfn->V_potential();

    // phi matrix (real-space <- AO basis mapping)
    // since grid is stored in blocks, we need to build a full phi matrix
    // from the individual blocks:
    std::shared_ptr<PointFunctions> points_func = potential->properties()[0];
    points_func->set_pointers(Da_,Db_);

    // determine number of grid points
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

    // what is the derivative level?
    deriv_   = functional->deriv();
    is_gga_  = functional->is_gga();
    is_meta_ = functional->is_meta();
    is_unpolarized_ = functional->is_unpolarized();

    super_phi_   = std::shared_ptr<Matrix>(new Matrix("SUPER PHI",phi_points_,nso_));

    if ( is_gga_ || is_meta_ ) {
        super_phi_x_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI X",phi_points_,nso_));
        super_phi_y_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Y",phi_points_,nso_));
        super_phi_z_ = std::shared_ptr<Matrix>(new Matrix("SUPER PHI Z",phi_points_,nso_));
    }

    grid_x_      = std::shared_ptr<Vector>(new Vector("GRID X",phi_points_));
    grid_y_      = std::shared_ptr<Vector>(new Vector("GRID Y",phi_points_));
    grid_z_      = std::shared_ptr<Vector>(new Vector("GRID Z",phi_points_));
    grid_w_      = std::shared_ptr<Vector>(new Vector("GRID W",phi_points_));

    // build phi matrix and derivative phi matrices
    BuildPhiMatrix(potential, points_func, "PHI",  super_phi_);

    if ( is_gga_ || is_meta_ ) {
        BuildPhiMatrix(potential, points_func, "PHI_X",super_phi_x_);
        BuildPhiMatrix(potential, points_func, "PHI_Y",super_phi_y_);
        BuildPhiMatrix(potential, points_func, "PHI_Z",super_phi_z_);
    }
}

void MCPDFTSolver::BuildPhiMatrix(std::shared_ptr<VBase> potential, std::shared_ptr<PointFunctions> points_func,
        std::string phi_type, std::shared_ptr<Matrix> myphi) {

    int nblocks = potential->nblocks();

    phi_points_ = 0;
    for (int myblock = 0; myblock < nblocks; myblock++) {
        std::shared_ptr<BlockOPoints> block = potential->get_block(myblock);
        points_func->compute_points(block);
        int npoints = block->npoints();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        double ** phi = points_func->basis_value(phi_type)->pointer();

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
                myphi->pointer()[phi_points_ + p][nug] = phi[p][nu];
            }
        }
        phi_points_ += npoints;
    }
    // grab AO->SO transformation matrix
    std::shared_ptr<Matrix> ao2so = reference_wavefunction_->aotoso();

    // transform one index of (derivative) phi matrix (AO->SO)
    std::shared_ptr<Matrix> temp (new Matrix(myphi));

    for (int p = 0; p < phi_points_; p++) {
        for (int sigma = 0; sigma < nso_; sigma++) {
            double dum = 0.0;
            for (int nu = 0; nu < nso_; nu++) {
                dum += myphi->pointer()[p][nu] * ao2so->pointer()[nu][sigma];
            }
            temp->pointer()[p][sigma] = dum;
        }
    }

    C_DCOPY(nso_*phi_points_,temp->pointer()[0],1,myphi->pointer()[0],1);

    // transform orbital index in phi matrices to MO basis
    // NOTE: this will not work if Ca and Cb are different.  it is imperative that
    // the v2rdm-casscf computation preceding mcpdft is run with reference = rhf/rohf
    for (int p = 0; p < phi_points_; p++) {
        for (int sigma = 0; sigma < nso_; sigma++) {
            double dum = 0.0;
            for (int nu = 0; nu < nso_; nu++) {
                dum += temp->pointer()[p][nu] * Ca_->pointer()[nu][sigma];
            }
            myphi->pointer()[p][sigma] = dum;
        }
    }
}

double MCPDFTSolver::compute_energy() {

    // first, let's try this without any symmetry

    // allocate memory for 1- and 2-RDM
    double * D1a  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1b  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D2aa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2bb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // read 2-RDM from disk
    ReadTPDM(D2aa,D2bb,D2ab,D1a,D1b);

    // with the 1-RDM, we can evaluate the kinetic, potential, and coulomb nergies

    // one-electron terms:

    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));

    SharedMatrix ha (new Matrix(mints->so_potential()));
    ha->add(mints->so_kinetic());
    ha->transform(Ca_);

    SharedMatrix hb (new Matrix(mints->so_potential()));
    hb->add(mints->so_kinetic());
    hb->transform(Cb_);

    double one_electron_energy = C_DDOT(nmo_*nmo_,D1a,1,ha->pointer()[0],1)
                               + C_DDOT(nmo_*nmo_,D1b,1,hb->pointer()[0],1);

    // coulomb energy should be computed using J object

    std::shared_ptr<Matrix> Ja = BuildJ(D1a,Ca_);
    std::shared_ptr<Matrix> Jb = BuildJ(D1b,Cb_);

    double caa = C_DDOT(nmo_*nmo_,D1a,1,Ja->pointer()[0],1);
    double cab = C_DDOT(nmo_*nmo_,D1a,1,Jb->pointer()[0],1);
    double cba = C_DDOT(nmo_*nmo_,D1b,1,Ja->pointer()[0],1);
    double cbb = C_DDOT(nmo_*nmo_,D1b,1,Jb->pointer()[0],1);

    double coulomb_energy = 0.5 * ( caa + cab + cba + cbb );

    // now, build rho(r), rho'(r), pi(r) and evaluate the MCPDFT xc energy

    double mcpdft_xc_energy = 0.0; 

    // build alpha- and beta-spin densities and gradients
    BuildRho(D1a,D1b);

    // build on-top pair density
    BuildPi(D2ab);

    // build R(r) = 4 * Pi(r) / rho(r)
    Build_R();
    
    // translate the alpha and beta densities and their corresponding gradients
    Translate();    
    
    // calculate the on-top energy
    mcpdft_xc_energy =  MCPDFTSolver::EX_LSDA(tr_rho_a_, tr_rho_b_, tr_zeta_) + MCPDFTSolver::EC_VWN3_RPA(tr_rho_a_, tr_rho_b_, tr_zeta_, tr_rs_);

    // print total energy and its components

    outfile->Printf("\n");
    outfile->Printf("    ==> Energetics <==\n");
    outfile->Printf("\n");

    outfile->Printf("        nuclear repulsion energy =    %20.12lf\n",molecule_->nuclear_repulsion_energy({0.0,0.0,0.0}));
    outfile->Printf("        one-electron energy =         %20.12lf\n",one_electron_energy);
    outfile->Printf("        coulomb energy =              %20.12lf\n",coulomb_energy);
    outfile->Printf("        exchange-correlation energy = %20.12lf\n",mcpdft_xc_energy);
    outfile->Printf("\n");

    double total_energy = molecule_->nuclear_repulsion_energy({0.0,0.0,0.0})+one_electron_energy+coulomb_energy+mcpdft_xc_energy;
    outfile->Printf("    * MCPDFT total energy =      %20.12lf\n",total_energy);
    outfile->Printf("\n");

    return total_energy;
}

void MCPDFTSolver::BuildPi(double * D2ab) {

    pi_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * pi_p = pi_->pointer();

    double ** phi = super_phi_->pointer();

    for (int p = 0; p < phi_points_; p++) {

        double dum = 0.0;

        // pi(r) = D(mu,nu; lambda,sigma) * phi(r,mu) * phi(r,nu) * phi(r,lambda) * phi(r,sigma)
        for (int mu = 0; mu < nmo_; mu++) {
            for (int nu = 0; nu < nmo_; nu++) {
                for (int lambda = 0; lambda < nmo_; lambda++) {
                    for (int sigma = 0; sigma < nmo_; sigma++) {

                        dum += phi[p][mu] * phi[p][lambda] * phi[p][sigma] * phi[p][nu] * D2ab[mu*nmo_*nmo_*nmo_ + nu*nmo_*nmo_ + lambda*nmo_ + sigma];

                    }
                }
            }
        }

        pi_p[p] = dum;
    }

}

void MCPDFTSolver::BuildRho(double * D1a, double * D1b) {

    rho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    zeta_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    rs_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double ** phi   = super_phi_->pointer();
    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * zeta_p = zeta_->pointer();
    double * rs_p = rs_->pointer();

    for (int p = 0; p < phi_points_; p++) {
        double duma   = 0.0;
        double dumb   = 0.0;
        for (int sigma = 0; sigma < nmo_; sigma++) {
            for (int nu = 0; nu < nmo_; nu++) {
                duma += phi[p][sigma] * phi[p][nu] * D1a[sigma*nmo_ + nu];
                dumb += phi[p][sigma] * phi[p][nu] * D1b[sigma*nmo_ + nu];
            }
        }
        rho_ap[p] = duma;
        rho_bp[p] = dumb;

        zeta_p[p] =  ( rho_ap[p] - rho_bp[p] ) / ( rho_ap[p] + rho_bp[p] );
        rs_p[p] = pow( 3.0 / ( 4.0 * M_PI * (rho_ap[p] + rho_bp[p]) ) , 1.0/3.0 );
    }

    if ( is_gga_ || is_meta_ ) {

        sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        tau_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
        tau_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));

        double ** phi_x = super_phi_x_->pointer();
        double ** phi_y = super_phi_y_->pointer();
        double ** phi_z = super_phi_z_->pointer();

        double * sigma_aap = sigma_aa_->pointer();
        double * sigma_bbp = sigma_bb_->pointer();
        double * sigma_abp = sigma_ab_->pointer();

        double * rho_a_xp = rho_a_x_->pointer();
        double * rho_b_xp = rho_b_x_->pointer();

        double * rho_a_yp = rho_a_y_->pointer();
        double * rho_b_yp = rho_b_y_->pointer();
 
        double * rho_a_zp = rho_a_z_->pointer();
        double * rho_b_zp = rho_b_z_->pointer();

        double * tau_ap = tau_a_->pointer();
        double * tau_bp = tau_b_->pointer();
        
        for (int p = 0; p < phi_points_; p++) {
            double duma_x = 0.0;
            double dumb_x = 0.0;
            double duma_y = 0.0;
            double dumb_y = 0.0;
            double duma_z = 0.0;
            double dumb_z = 0.0;
            double dumta = 0.0;
            double dumtb = 0.0;

            for (int sigma = 0; sigma < nmo_; sigma++) {
                for (int nu = 0; nu < nmo_; nu++) {
                    duma_x += ( phi_x[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_x[p][nu] ) * D1a[sigma*nmo_ + nu];
                    dumb_x += ( phi_x[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_x[p][nu] ) * D1b[sigma*nmo_ + nu];

                    duma_y += ( phi_y[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_y[p][nu] ) * D1a[sigma*nmo_ + nu];
                    dumb_y += ( phi_y[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_y[p][nu] ) * D1b[sigma*nmo_ + nu];

                    duma_z += ( phi_z[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_z[p][nu] ) * D1a[sigma*nmo_ + nu];
                    dumb_z += ( phi_z[p][sigma] * phi[p][nu] + phi[p][sigma] * phi_z[p][nu] ) * D1b[sigma*nmo_ + nu];

                    dumta += (phi_x[p][sigma] * phi_x[p][nu] + phi_y[p][sigma] * phi_y[p][nu] + phi_z[p][sigma] * phi_z[p][nu]) * D1a[sigma*nmo_ + nu];
                    dumtb += (phi_x[p][sigma] * phi_x[p][nu] + phi_y[p][sigma] * phi_y[p][nu] + phi_z[p][sigma] * phi_z[p][nu]) * D1b[sigma*nmo_ + nu];
                }
            }
            rho_a_xp[p] = duma_x;
            rho_b_xp[p] = dumb_x;

            rho_a_yp[p] = duma_y;
            rho_b_yp[p] = dumb_y;

            rho_a_zp[p] = duma_z;
            rho_b_zp[p] = dumb_z;

            sigma_aap[p] = pow(rho_a_xp[p],2.0) + pow(rho_a_yp[p],2.0) + pow(rho_a_zp[p],2.0);
            sigma_bbp[p] = pow(rho_b_xp[p],2.0) + pow(rho_b_yp[p],2.0) + pow(rho_b_zp[p],2.0);
            sigma_abp[p] = ( rho_a_xp[p] * rho_b_xp[p] ) +  ( rho_a_yp[p] * rho_b_yp[p] ) + ( rho_a_zp[p] * rho_b_zp[p] );

            tau_ap[p] = dumta;
            tau_bp[p] = dumtb;

        }
    }
}

std::shared_ptr<Matrix> MCPDFTSolver::BuildJ(double * D, std::shared_ptr<Matrix> C) {

    // get primary basis:
    std::shared_ptr<BasisSet> primary = reference_wavefunction_->get_basisset("ORBITAL");

    // get auxiliary basis:
    std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");

    // JK object (note this is hard-coded to use density fitting ...)
    std::shared_ptr<DFJK> jk = (std::shared_ptr<DFJK>)(new DFJK(primary,auxiliary));

    // memory for jk (say, 50% of what is available)
    jk->set_memory(0.5 * Process::environment.get_memory());

    // integral cutoff
    jk->set_cutoff(options_.get_double("INTS_TOLERANCE"));

    jk->set_do_J(true);
    jk->set_do_K(false);
    jk->set_do_wK(false);
    jk->set_omega(false);

    jk->initialize();

    std::vector<SharedMatrix>& C_left  = jk->C_left();
    std::vector<SharedMatrix>& C_right = jk->C_right();

    // use real and imaginary C matrices

    std::shared_ptr<Matrix> myC (new Matrix(C) );

    myC->zero();

    C_left.clear();
    C_right.clear();

    C_DCOPY(nmo_*nmo_,D,1,Da_->pointer()[0],1);
    myC->gemm('t','n',1.0,Da_,C,0.0);
    myC->transpose_this();
    C_left.push_back(myC);
    C_right.push_back(C);

    // Let jk compute for the given C_left/C_right

    jk->compute();

    std::shared_ptr<Matrix> J = jk->J()[0];
    J->transform(C);

    return J;
}

void MCPDFTSolver::Build_R(){

    R_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * pi_p = pi_->pointer();
    double * R_p = R_->pointer();
 
    for (int p = 0; p < phi_points_; p++) {
        
        double rho = rho_ap[p] + rho_bp[p];
        
        R_p[p] = 4 * pi_p[p] / (rho * rho);         
    }

} 

void MCPDFTSolver::Translate(){

    tr_rho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    tr_rho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    
    tr_zeta_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    tr_rs_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();

    double * pi_p = pi_->pointer();
    double * R_p = R_->pointer();
    double * tr_zeta_p = tr_zeta_->pointer();
    double * tr_rs_p = tr_rs_->pointer();

    double * tr_rho_ap = tr_rho_a_->pointer();
    double * tr_rho_bp = tr_rho_b_->pointer();

    for (int p = 0; p < phi_points_; p++) {

        double rho = rho_ap[p] + rho_bp[p];

        tr_rho_ap[p] = ( (R_p[p] > 1.0) ? (rho/2.0) : (rho/2.0) * ( 1.0 + sqrt(1.0 - R_p[p]) ) );
        tr_rho_bp[p] = ( (R_p[p] > 1.0) ? (rho/2.0) : (rho/2.0) * ( 1.0 - sqrt(1.0 - R_p[p]) ) );

        tr_zeta_p[p] =  ( tr_rho_ap[p] - tr_rho_bp[p] ) / ( tr_rho_ap[p] + tr_rho_bp[p] );
        tr_rs_p[p] = pow( 3.0 / ( 4.0 * M_PI * (tr_rho_ap[p] + tr_rho_bp[p]) ) , 1.0/3.0 );
    }

    if ( is_gga_ || is_meta_ ) {

       tr_rho_a_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_rho_b_x_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

       tr_rho_a_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_rho_b_y_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

       tr_rho_a_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_rho_b_z_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       
       tr_sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
       tr_sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

       double * rho_a_xp = rho_a_x_->pointer();
       double * rho_b_xp = rho_b_x_->pointer();

       double * rho_a_yp = rho_a_y_->pointer();
       double * rho_b_yp = rho_b_y_->pointer();
 
       double * rho_a_zp = rho_a_z_->pointer();
       double * rho_b_zp = rho_b_z_->pointer();

       double * tr_rho_a_xp = tr_rho_a_x_->pointer();
       double * tr_rho_b_xp = tr_rho_b_x_->pointer();

       double * tr_rho_a_yp = tr_rho_a_y_->pointer();
       double * tr_rho_b_yp = tr_rho_b_y_->pointer();
 
       double * tr_rho_a_zp = tr_rho_a_z_->pointer();
       double * tr_rho_b_zp = tr_rho_b_z_->pointer();
       
       double * tr_sigma_aap = tr_sigma_aa_->pointer();
       double * tr_sigma_abp = tr_sigma_ab_->pointer();
       double * tr_sigma_bbp = tr_sigma_bb_->pointer();

       for (int p = 0; p < phi_points_; p++) {

           tr_rho_a_xp[p] = ( (R_p[p] > 1.0) ? (rho_a_xp[p]/2.0) : (rho_a_xp[p]/2.0) * ( 1.0 + sqrt(1.0 - R_p[p]) ) );
           tr_rho_b_xp[p] = ( (R_p[p] > 1.0) ? (rho_b_xp[p]/2.0) : (rho_b_xp[p]/2.0) * ( 1.0 - sqrt(1.0 - R_p[p]) ) );
           
           tr_rho_a_yp[p] = ( (R_p[p] > 1.0) ? (rho_a_yp[p]/2.0) : (rho_a_yp[p]/2.0) * ( 1.0 + sqrt(1.0 - R_p[p]) ) );
           tr_rho_b_yp[p] = ( (R_p[p] > 1.0) ? (rho_b_yp[p]/2.0) : (rho_b_yp[p]/2.0) * ( 1.0 - sqrt(1.0 - R_p[p]) ) );
           
           tr_rho_a_zp[p] = ( (R_p[p] > 1.0) ? (rho_a_zp[p]/2.0) : (rho_a_zp[p]/2.0) * ( 1.0 + sqrt(1.0 - R_p[p]) ) );
           tr_rho_b_zp[p] = ( (R_p[p] > 1.0) ? (rho_b_zp[p]/2.0) : (rho_b_zp[p]/2.0) * ( 1.0 - sqrt(1.0 - R_p[p]) ) );

           tr_sigma_aap[p] = (tr_rho_a_xp[p] * tr_rho_a_xp[p]) + (tr_rho_a_yp[p] * tr_rho_a_yp[p]) + (tr_rho_a_zp[p] * tr_rho_a_zp[p]);  
           tr_sigma_abp[p] = (tr_rho_a_xp[p] * tr_rho_b_xp[p]) + (tr_rho_a_yp[p] * tr_rho_b_yp[p]) + (tr_rho_a_zp[p] * tr_rho_b_zp[p]);  
           tr_sigma_bbp[p] = (tr_rho_b_xp[p] * tr_rho_b_xp[p]) + (tr_rho_b_yp[p] * tr_rho_b_yp[p]) + (tr_rho_b_zp[p] * tr_rho_b_zp[p]);  
       }
    }
}

}} // end of namespaces



