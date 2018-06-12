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


#ifndef MCPDFT_SOLVER_H
#define MCPDFT_SOLVER_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "psi4/libmints/wavefunction.h"

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"

// for grid
#include "psi4/libfock/points.h"
#include "psi4/libfock/cubature.h"

#define PSIF_SO_TEI           33
#define PSIF_SO_ERF_TEI       36
#define PSIF_SO_ERFC_TEI      37
#define PSIF_MO_TEI           72
#define PSIF_DCC_QMO          268
#define PSIF_V2RDM_CHECKPOINT 269
#define PSIF_V2RDM_D2AA       270
#define PSIF_V2RDM_D2AB       271
#define PSIF_V2RDM_D2BB       272
#define PSIF_V2RDM_D3AAA      273
#define PSIF_V2RDM_D3AAB      274
#define PSIF_V2RDM_D3BBA      275
#define PSIF_V2RDM_D3BBB      276
#define PSIF_V2RDM_D1A        277
#define PSIF_V2RDM_D1B        278

namespace psi{ namespace mcpdft{

struct tpdm {
    int i;
    int j;
    int k;
    int l;
    double val;
};

struct opdm {
    int i;
    int j;
    double val;
};

class MCPDFTSolver: public Wavefunction{

  public:

    MCPDFTSolver(std::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
    ~MCPDFTSolver();
    void common_init();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return same_a_b_orbs_; }
    virtual bool same_a_b_dens() const { return same_a_b_dens_; }
    void ReadCIOPDM(double* D, const char* fileName); 
    void ReadCITPDM(double* D, const char* fileName); 
    void PrintOPDM(double* D); 
    void PrintTPDM(double* D); 

  protected:

    /// reference energy
    double reference_energy_;

    /// two-electron part of the reference energy
    double two_electron_energy_;

    /// mp2 correlation energy for double-hybrids
    double mp2_corr_energy_;

    /// list of orbital symmetries
    int * symmetry_;

    /// offset for orbitals in each irrep
    int * pitzer_offset_;

    /// level of differentiation 
    int deriv_;

    /// is gga?
    bool is_gga_;

    /// is df?
    bool is_df_;

    /// is meta?
    bool is_meta_;

    /// available memory
    long int available_memory_;

    /// is unpolarized (restricted)?
    bool is_unpolarized_;

    /// number of grid points_
    long int phi_points_;

    /// number of auxilliary basis functions
    int nQ_;

    /// phi matrix
    std::shared_ptr<Matrix> super_phi_;

    /// d phi / dx matrix
    std::shared_ptr<Matrix> super_phi_x_;

    /// d phi / dy matrix
    std::shared_ptr<Matrix> super_phi_y_;

    /// d phi / dz matrix
    std::shared_ptr<Matrix> super_phi_z_;

    /// d2 phi / dx2 matrix
    std::shared_ptr<Matrix> super_phi_xx_;

    /// d phi / dxdy matrix
    std::shared_ptr<Matrix> super_phi_xy_;

    /// d2 phi / dxdz matrix
    std::shared_ptr<Matrix> super_phi_xz_;

    /// d2 phi / dy2 matrix
    std::shared_ptr<Matrix> super_phi_yy_;

    /// d2 phi / dydz matrix
    std::shared_ptr<Matrix> super_phi_yz_;

    /// d2 phi / dz2 matrix
    std::shared_ptr<Matrix> super_phi_zz_;

    /// gamma_aa matrix
    std::shared_ptr<Matrix> super_gamma_aa_; 

    /// gamma_ab matrix
    std::shared_ptr<Matrix> super_gamma_ab_;   

    /// gamma_bb matrix
    std::shared_ptr<Matrix> super_gamma_bb_;   

    /// tau_a matrix
    std::shared_ptr<Matrix> super_tau_a_;

    /// tau_b matrix
    std::shared_ptr<Matrix> super_tau_b_;   

    /// grid x values
    std::shared_ptr<Vector> grid_x_;

    /// grid y values
    std::shared_ptr<Vector> grid_y_;

    /// grid z values
    std::shared_ptr<Vector> grid_z_;

    /// grid weights
    std::shared_ptr<Vector> grid_w_;
    
    /// inner product of alpha density gradient with itself (gamma_aa)
    std::shared_ptr<Vector> sigma_aa_;

    /// inner product of beta density gradient with itself (gamma_bb)
    std::shared_ptr<Vector> sigma_bb_;

    /// inner product of alpha density gradient with beta density gradient (gamma_ab)
    std::shared_ptr<Vector> sigma_ab_;

    /// a function to build phi/phi_x/... note: the orbital labels are in the AO basis (no symmetry)
    void BuildPhiMatrixAO(std::shared_ptr<VBase> potential, std::shared_ptr<PointFunctions> points_func,
            std::string phi_type, std::shared_ptr<Matrix> myphi);

    /// transform the orbital labels in phi/phi_x/... from the AO to the MO basis
    void TransformPhiMatrixAOMO(std::shared_ptr<Matrix> phi_in, std::shared_ptr<Matrix> phi_out);

    /// read v2RDM 2-RDM from disk
    void ReadTPDM();

    /// read v2RDM 1-RDM from disk
    void ReadOPDM();

    /// read CI 2-RDM from disk
    void ReadCITPDM(double* D, const char* fileName); 

    /// read CI 1-RDM from disk
    void ReadCIOPDM(std::shared_ptr<Matrix> D, const char* fileName); 

    /// build coulomb/exchange matrix
    std::vector< std::shared_ptr<Matrix> > BuildJK();

    /// build coulomb/exchange super matrix
    double BuildErfCoulombEnergy();

    /// alpha-spin density
    std::shared_ptr<Vector> rho_a_;

    /// beta-spin density
    std::shared_ptr<Vector> rho_b_;

    /// total density
    std::shared_ptr<Vector> rho_;
 
    /// alpha-spin density gradient x
    std::shared_ptr<Vector> rho_a_x_;

    /// beta-spin density gradient x
    std::shared_ptr<Vector> rho_b_x_;

    /// alpha-spin density gradient y
    std::shared_ptr<Vector> rho_a_y_;

    /// beta-spin density gradient y
    std::shared_ptr<Vector> rho_b_y_;

    /// alpha-spin density gradient z
    std::shared_ptr<Vector> rho_a_z_;

    /// beta-spin density gradient z
    std::shared_ptr<Vector> rho_b_z_;
    
    /// translated alpha-spin density
    std::shared_ptr<Vector> tr_rho_a_;

    /// translated beta-spin density
    std::shared_ptr<Vector> tr_rho_b_;

    /// translated alpha-spin density gradient x
    std::shared_ptr<Vector> tr_rho_a_x_;

    /// translated beta-spin density gradient x
    std::shared_ptr<Vector> tr_rho_b_x_;

    /// translated alpha-spin density gradient y
    std::shared_ptr<Vector> tr_rho_a_y_;

    /// translated beta-spin density gradient y
    std::shared_ptr<Vector> tr_rho_b_y_;

    /// translated alpha-spin density gradient z
    std::shared_ptr<Vector> tr_rho_a_z_;

    /// translated beta-spin density gradient z
    std::shared_ptr<Vector> tr_rho_b_z_;
    
    /// translated inner product of alpha density gradient with itself (gamma_aa)
    std::shared_ptr<Vector> tr_sigma_aa_;

    /// translated inner product of beta density gradient with itself (gamma_bb)
    std::shared_ptr<Vector> tr_sigma_bb_;

    /// translated inner product of alpha density gradient with beta density gradient (gamma_ab)
    std::shared_ptr<Vector> tr_sigma_ab_;
    
    /// translated inner product of total density gradient
    std::shared_ptr<Vector> tr_sigma_;

    /// R(r) = 4 * Pi(r) / rho(r) ^ 2
    std::shared_ptr<Vector> R_;

    /// zeta factor zeta = ( rho_a(r) - rho_b(r) ) / ( rho_a(r) + rho_b(r) )
    std::shared_ptr<Vector> zeta_;
   
    /// spin magnetization density
    std::shared_ptr<Vector> m_;
    
    /// translated spin magnetization density
    std::shared_ptr<Vector> tr_m_;

    /// translated zeta factor
    std::shared_ptr<Vector> tr_zeta_;

    /// effective radius of density
    std::shared_ptr<Vector> rs_;

    /// translated effective radius of density
    std::shared_ptr<Vector> tr_rs_;

    /// tau kinetic energy of alpha electrons
    std::shared_ptr<Vector> tau_a_;

    /// tau kinetic energy of beta electrons
    std::shared_ptr<Vector> tau_b_;

    /// the on-top pair density
    std::shared_ptr<Vector> pi_;

    /// x-component of the gradient of the on-top pair density
    std::shared_ptr<Vector> pi_x_;

    /// y-component of the gradient of the on-top pair density
    std::shared_ptr<Vector> pi_y_;

    /// z-component of the gradient of the on-top pair density
    std::shared_ptr<Vector> pi_z_;

    /// build spin densities and gradients
    void BuildRho();

    /// build spin densities and gradients using only non-zero elements of OPDM
    void BuildRhoFast(opdm * D1a, opdm * D1b, int na, int nb);

    /// build on-top pair density
    void BuildPi(double * D2ab);
 
    /// build on-top pair density using only non-zero elements of TPDM
    void BuildPiFast(tpdm * D2ab, int nab);
 
    /// build gradient of the on-top pair density
    void Build_Grad_Pi();

    /// build R(r) = 4 * Pi(r) / rho(r) ^ 2
    void Build_R();
    
    /// build translator function of density and its gradient
    void Translate();

    /// build full-translator function of density and its gradient
    void Fully_Translate();

    /// build G spin-interpolation formula
    double Gfunction(double r, double A, double a1, double b1, double b2, double b3, double b4, double p);

    //############################################################
    //############ Exchange functions' declarations ##############
    //############################################################

    /// build EX_LDA(rho)
    double EX_LDA(std::shared_ptr<Vector> rho_a, std::shared_ptr<Vector> rho_b);

    /// build EX_LSDA(rho_sigma)
    double EX_LSDA_Sigma(std::shared_ptr<Vector> rho_sigma);

    /// build EX_LSDA(rho_a, rho_b, zeta)
    double EX_LSDA(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> ZETA);

    /// build EX_LSDA(rho_a, rho_b)
    double EX_LSDA(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B);

    /// build EX_B86_MGC()
    double EX_B86_MGC();

    /// build EX_B88()
    double EX_B88(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB);
    double EX_B88_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB);

    /// build EX_PBE()
    double EX_PBE(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB);
    double EX_PBE_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB);

    /// build EX_wPBE()
    double EX_wPBE_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB);

    /// build EX_RPBE()
    double EX_RPBE();

    /// build EX_UPBE()
    double EX_UPBE();

    //############################################################
    //########### Correlation functions' declarations ############
    //############################################################

    /// build EC_B88()
    double EC_B88_OP(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB);

    /// build EC_LYP()
    double EC_LYP_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, 
                  std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_AB, std::shared_ptr<Vector> SIGMA_BB);

    /// build EC_PBE()
    double EC_PBE(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, 
                  std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_AB, std::shared_ptr<Vector> SIGMA_BB);
    double EC_PBE_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, 
                  std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_AB, std::shared_ptr<Vector> SIGMA_BB);

    /// build EC_VWN3_RPA()
    double EC_VWN3_RPA(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> ZETA, std::shared_ptr<Vector> RS);
    double EC_VWN3_RPA_III(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B);
    // double EC_VWN3_RPA_III(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> ZETTA);
};

}} // end of namespaces

#endif
