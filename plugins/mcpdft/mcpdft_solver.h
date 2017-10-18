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

#define PSIF_DCC_QMO          268
#define PSIF_V2RDM_CHECKPOINT 269
#define PSIF_V2RDM_D2AA       270
#define PSIF_V2RDM_D2AB       271
#define PSIF_V2RDM_D2BB       272
#define PSIF_V2RDM_D3AAA      273
#define PSIF_V2RDM_D3AAB      274
#define PSIF_V2RDM_D3BBA      275
#define PSIF_V2RDM_D3BBB      276

namespace psi{ namespace mcpdft{

class MCPDFTSolver: public Wavefunction{

  public:

    MCPDFTSolver(std::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
    ~MCPDFTSolver();
    void common_init();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return same_a_b_orbs_; }
    virtual bool same_a_b_dens() const { return same_a_b_dens_; }

  protected:

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

    /// read 2-RDM from disk
    void ReadTPDM(double * D2aa, double * D2bb, double * D2ab, double * D1a, double * D1b);

    /// build coulomb matrix
    std::shared_ptr<Matrix> BuildJ(double * D, std::shared_ptr<Matrix> C);

    /// alpha-spin density
    std::shared_ptr<Vector> rho_a_;

    /// beta-spin density
    std::shared_ptr<Vector> rho_b_;

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

    /// the on-top pair density
    std::shared_ptr<Vector> pi_;

    /// build spin densities and gradients
    void BuildRho(double * D1a, double * D1b);

    /// build on-top pair density
    void BuildPi(double * D2ab);
};

}} // end of namespaces

#endif
