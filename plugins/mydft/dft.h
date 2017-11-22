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

#ifndef DFT_SOLVER_H
#define DFT_SOLVER_H

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

namespace psi{ namespace mydft{

class DFTSolver: public Wavefunction{

  public:

    DFTSolver(std::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
    ~DFTSolver();
    void common_init();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return same_a_b_orbs_; }
    virtual bool same_a_b_dens() const { return same_a_b_dens_; }

  protected:

    /// the nuclear repulsion energy
    double enuc_;

    /// evaluate the orbital gradient
    std::shared_ptr<Matrix> OrbitalGradient(std::shared_ptr<Matrix> D,
                                            std::shared_ptr<Matrix> F,
                                            std::shared_ptr<Matrix> Shalf);

    /// xc potential matrices
    std::shared_ptr<Matrix> Va_;
    std::shared_ptr<Matrix> Vb_;

    /// explicitly evaluate exchange correlation energy 
    double ExchangeCorrelationEnergy();

    /// number of grid points_
    long int phi_points_;

    /// phi matrix
    std::shared_ptr<Matrix> super_phi_;

    /// d phi / dx matrix
    std::shared_ptr<Matrix> super_phi_x_;

    /// d phi / dy matrix
    std::shared_ptr<Matrix> super_phi_y_;

    /// d phi / dz matrix
    std::shared_ptr<Matrix> super_phi_z_;

    /// grid x values
    std::shared_ptr<Vector> grid_x_;

    /// grid y values
    std::shared_ptr<Vector> grid_y_;

    /// grid z values
    std::shared_ptr<Vector> grid_z_;

    /// grid weights
    std::shared_ptr<Vector> grid_w_;

    /// a function to build phi/phi_x/...
    void BuildPhiMatrix(std::shared_ptr<VBase> potential, std::shared_ptr<PointFunctions> points_func,
            std::string phi_type, std::shared_ptr<Matrix> myphi);

    /// S^{-1/2}
    std::shared_ptr<Matrix> Shalf_;
    std::shared_ptr<Matrix> Shalf2;

    /// build G spin-interpolation formula
    double Gfunction(double r, double A, double a1, double b1, double b2, double b3, double b4, double p);

    /// zeta factor zeta = ( rho_a(r) - rho_b(r) ) / ( rho_a(r) + rho_b(r) )
    std::shared_ptr<Vector> zeta_;

    /// effective radius of density
    std::shared_ptr<Vector> rs_;

    /// tau kinetic energy of alpha electrons
    std::shared_ptr<Vector> tau_a_;
    
    /// tau kinetic energy of beta electrons
    std::shared_ptr<Vector> tau_b_;

    /// build spin densities and gradients
    void BuildRho(std::shared_ptr<Matrix> Da, std::shared_ptr<Matrix> Db);

    //############################################################
    //############ Exchange functions' declarations ##############
    //############################################################

    /// build EX_LDA(rho)
    double EX_LDA(std::shared_ptr<Vector> rho_a, std::shared_ptr<Vector> rho_b);

    /// build EX_LSDA(rho_sigma)
    double EX_LSDA_Sigma(std::shared_ptr<Vector> rho_sigma);

    /// build EX_LSDA(rho_a, rho_b, zeta)
    double EX_LSDA(std::shared_ptr<Vector> rho_a, std::shared_ptr<Vector> rho_b, std::shared_ptr<Vector> zeta);
    
    /// build EX_B86_MGC()
    double EX_B86_MGC();

    /// build EX_B88()
    double EX_B88();

    /// build EX_PBE()
    double EX_PBE();

    /// build EX_RPBE()
    double EX_RPBE();

    /// build EX_UPBE()
    double EX_UPBE();

    //############################################################
    //########### Correlation functions' declarations ############
    //############################################################

    /// build EC_B88()
    double EC_B88_OP();
   
    /// build EC_PBE()
    double EC_PBE();
    
    /// build EC_VWN3_RPA()
    double EC_VWN3_RPA();
     

};

}}


#endif

