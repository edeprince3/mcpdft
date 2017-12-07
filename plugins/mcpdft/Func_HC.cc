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

#include "mcpdft_solver.h"

namespace psi{ namespace mcpdft {

    //###########################################################
    //# The exchange and correlation functional implementations #
    //###########################################################

double MCPDFTSolver::Gfunction(double r, double A, double a1, double b1, double b2, double b3, double b4, double p) {

    double G = -2.0 * A * (1.0 + a1 * r) * log( 1.0 + pow( 2.0 * A * ( b1 * sqrt(r) + b2 * r + b3 * pow(r,3.0/2.0) + b4 * pow(r, p+1.0 ) )  ,-1.0 ) );

    return G;

}
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++ Exchange functionals ++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


double MCPDFTSolver::EX_LDA(std::shared_ptr<Vector> rho_a, std::shared_ptr<Vector> rho_b){
    
    const double alpha = (2.0/3.0);      // Slater value
    const double Cx = (9.0/8.0) * alpha * pow(3.0/M_PI,1.0/3.0);
    
    double * rho_ap = rho_a->pointer();
    double * rho_bp = rho_b->pointer();

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rho = rho_ap[p] + rho_bp[p];

        exc += -Cx * pow( rho, 4.0/3.0) * grid_w_->pointer()[p]; 
    }
    return exc;
}

double MCPDFTSolver::EX_LSDA(std::shared_ptr<Vector> rho_a, std::shared_ptr<Vector> rho_b){
    
    const double alpha = (2.0/3.0);      // Slater value
    const double Cx = (9.0/8.0) * alpha * pow(3.0/M_PI,1.0/3.0);

    double * rho_ap = rho_a->pointer();
    double * rho_bp = rho_b->pointer();
    
    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double exa = pow(2.0,1.0/3.0) * Cx * pow( rho_ap[p], 4.0/3.0) ;
        double exb = pow(2.0,1.0/3.0) * Cx * pow( rho_bp[p], 4.0/3.0) ;
        double ex_LSDA = exa + exb;
        exc += - ex_LSDA * grid_w_->pointer()[p]; 
    }
    return exc;
}

double MCPDFTSolver::EX_LSDA(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> ZETA){
    
    const double alpha = (2.0/3.0);      // Slater value
    const double Cx = (9.0/8.0) * alpha * pow(3.0/M_PI,1.0/3.0);

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();
    double * zeta_p = ZETA->pointer();
    
    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {
   
        // build f(zeta) weight factor where f(0) = 0 and f(1) = 1
        double fZet =  ( pow( (1.0 + zeta_p[p]) ,4.0/3.0 ) + pow( (1.0 - zeta_p[p]) ,4.0/3.0) - 2.0) / ( 2.0 * pow(2,1.0/3.0) - 2.0 );

        double ex0 = Cx * pow( (rho_ap[p] + rho_bp[p]), 1.0/3.0) ;
        double ex1 = pow(2.0,1.0/3.0) * ex0;
        double ex_LSDA = ex0 + (ex1 - ex0) * fZet;
        exc += - ex_LSDA * (rho_ap[p] + rho_bp[p]) * grid_w_->pointer()[p]; 
    }
    return exc;
}

double MCPDFTSolver::EX_B86_MGC(){
    
    const double Cx = 0.73855876638202240586; 
    const double beta = 0.00375;
    const double c = pow(2.0,1.0/3.0) * Cx;
   
    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * sigma_aap = sigma_aa_->pointer();
    double * sigma_bbp = sigma_bb_->pointer();
 
    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {
        
        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rhoa_43 = pow( rhoa, 4.0/3.0); 
        double rhob_43 = pow( rhob, 4.0/3.0); 
        double Xa = sqrt(sigma_aap[p]) / rhoa_43;
        double Xb = sqrt(sigma_bbp[p]) / rhob_43;
        double Xa_2 = Xa * Xa;
        double Xb_2 = Xb * Xb;
         
        exc += ( -c * rhoa_43 - (beta * Xa_2 * rhoa_43) / pow(1.0 + 0.007 * Xa_2,4.0/5.0) ) * grid_w_->pointer()[p]; 
        exc += ( -c * rhob_43 - (beta * Xb_2 * rhob_43) / pow(1.0 + 0.007 * Xb_2,4.0/5.0) ) * grid_w_->pointer()[p]; 
    }
    return exc;
}

double MCPDFTSolver::EX_B88(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB){
    
    const double Cx = 0.73855876638202240586; 
    const double beta = 0.0042;
    const double c = pow(2.0,1.0/3.0) * Cx;
 
    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();
    double * sigma_aap = SIGMA_AA->pointer();
    double * sigma_bbp = SIGMA_BB->pointer();
 
    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {
        
        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rhoa_43 = pow( rhoa, 4.0/3.0); 
        double rhob_43 = pow( rhob, 4.0/3.0); 
        double Xa = sqrt(sigma_aap[p]) / rhoa_43;
        double Xb = sqrt(sigma_bbp[p]) / rhob_43;
        double Xa_2 = Xa * Xa;
        double Xb_2 = Xb * Xb;

        exc += -rhoa_43 * ( c + (beta * Xa_2) / (1.0 + 6.0 * beta * Xa * asinh(Xa)) ) * grid_w_->pointer()[p]; 
        exc += -rhob_43 * ( c + (beta * Xb_2) / (1.0 + 6.0 * beta * Xb * asinh(Xb)) ) * grid_w_->pointer()[p]; 
    }
    return exc;
}

double MCPDFTSolver::EX_PBE(){
    
    const double delta = 0.06672455060314922;
    const double MU = (1.0/3.0) * delta * M_PI * M_PI;
    const double KAPPA = 0.804;

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * sigma_aap = sigma_aa_->pointer();
    double * sigma_bbp = sigma_bb_->pointer();
    
    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {
   
        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rhoa_43 = pow( rhoa, 4.0/3.0); 
        double rhob_43 = pow( rhob, 4.0/3.0); 
        double Xa = sqrt(sigma_aap[p]) / rhoa_43;
        double Xb = sqrt(sigma_bbp[p]) / rhob_43;
        
        double Sa = (Xa * pow(6.0, 2.0/3.0)) / (12.0 * pow(M_PI, 2.0/3.0));
        double Sb = (Xb * pow(6.0, 2.0/3.0)) / (12.0 * pow(M_PI, 2.0/3.0));
        
        double Fsa = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(Sa,2.0)) / KAPPA ), -1.0 );
        double Fsb = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(Sb,2.0)) / KAPPA ), -1.0 );
       
        auto E = [](double rhos, double Fss) -> double{

                 double temp = -0.75 * pow(3.0, 1.0/3.0) * pow(M_PI, 2.0/3.0) * pow(rhos,4.0/3.0) * Fss / M_PI;
                 return temp;
        };
 
        double EX_GGAa = 0.5 * E(2.0*rhoa,Fsa);
        double EX_GGAb = 0.5 * E(2.0*rhob,Fsb);
       
        exc += ( EX_GGAa + EX_GGAb ) * grid_w_->pointer()[p]; 
    }
    return exc;
}

// double MCPDFTSolver::EX_PBE(){
// 
//     const double alpha = (2.0/3.0);      // Slater value
//     const double Cx = (9.0/8.0) * alpha * pow(3.0/M_PI,1.0/3.0);
//     const double MU = 0.2195149727645171;
//     const double KAPPA = 0.804;
// 
//     double * rho_ap = rho_a_->pointer();
//     double * rho_bp = rho_b_->pointer();
// 
//     double * rho_a_xp = rho_a_x_->pointer();
//     double * rho_b_xp = rho_b_x_->pointer();
// 
//     double * rho_a_yp = rho_a_y_->pointer();
//     double * rho_b_yp = rho_b_y_->pointer();
// 
//     double * rho_a_zp = rho_a_z_->pointer();
//     double * rho_b_zp = rho_b_z_->pointer();
// 
//     double * zeta_p = zeta_->pointer();
//     double * sigma_aap = sigma_aa_->pointer();
//     double * sigma_abp = sigma_ab_->pointer();
//     double * sigma_bbp = sigma_bb_->pointer();
// 
//     double exc = 0.0;
//     for (int p = 0; p < phi_points_; p++) {
// 
//         double sig = sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p];
//         // double sig = sigma_aap[p]; 
//         double rhoa = rho_ap[p];
//         double rhob = rho_bp[p];
//         double rho = rhoa + rhob;
//         // local fermi wave vector
//         double kf = pow( ( 3.0 * pow(M_PI,2.0) * rho), 1.0/3.0);
//         // double kfa = pow( ( 3.0 * pow(M_PI,2.0) * rhoa ) , 1.0/3.0);
//         // double kfb = pow( ( 3.0 * pow(M_PI,2.0) * rhob ) , 1.0/3.0);
// 
//         // double EXa = -(3.0 * kfa) / (4.0 * M_PI); 
//         // double EXb = -(3.0 * kfb) / (4.0 * M_PI);
//         // double absDelRho = sqrt( ( pow( (rho_a_xp[p] + rho_b_xp[p]) ,2.0) + pow( (rho_a_yp[p] + rho_b_yp[p]) ,2.0) + pow( (rho_a_zp[p] + rho_b_zp[p]) ,2.0) ) );
//         double absDelRho = sqrt( sig );
// 
//         // double Sa = sqrt(4.0 * grada_2) / (2.0 * kfa * rhoa); 
//         // double Sb = sqrt(4.0 * gradb_2) / (2.0 * kfb * rhob);
// 
//         double s = absDelRho / ( 2.0 * kf * rho);
//         // double sa = absDelRho / ( 2.0 * kfa *2.0* rho_ap[p] );
//         // double sb = absDelRho / ( 2.0 * kfb *2.0* rho_bp[p] );
// 
//         double Fs = pow( (1.0 + 1.296 * pow(s,2.0) + 14.0 * pow(s,4.0) + 0.2 * pow(s,6.0) ) , 1.0/15.0 ) ;
//         // double Fsa = pow( (1.0 + 1.296 * pow(sa,2.0) + 14.0 * pow(sa,4.0) + 0.2 * pow(sa,6.0) ) , 1.0/15.0 );
//         // double Fsb = pow( (1.0 + 1.296 * pow(sb,2.0) + 14.0 * pow(sb,4.0) + 0.2 * pow(sb,6.0) ) , 1.0/15.0 );
//         // double Fs = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(s,2.0)) / KAPPA ), -1.0 );
//         // double Fsa = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(Sa,2.0)) / KAPPA ), -1.0 );
//         // double Fsb = 1.0 + KAPPA - KAPPA * pow( (1.0 + (MU * pow(Sb,2.0)) / KAPPA ), -1.0 );
// 
//         // double EX_GGAa = rhoa * EXa * Fsa;
//         // double EX_GGAb = rhob * EXb * Fsb;
//         exc += -Cx * pow( rho, 4.0/3.0) * Fs * grid_w_->pointer()[p];
//         // exc += -pow(2.0,1.0/3.0) * Cx * ( pow( rho_ap[p], 4.0/3.0) + pow( rho_bp[p], 4.0/3.0) ) * Fs * grid_w_->pointer()[p]; 
//         // exc += -pow(2.0,1.0/3.0) * Cx * ( pow( rho_ap[p], 4.0/3.0) * Fsa + pow( rho_bp[p], 4.0/3.0) * Fsb ) * grid_w_->pointer()[p]; 
//         // exc += -0.5 * ( -3.0 *2.0* kfa/(4*M_PI) * rho_ap[p] * Fsa - 3.0 *2.0* kfb / (4*M_PI) * rho_bp[p] * Fsb ) * grid_w_->pointer()[p]; 
//         // exc += 0.5 * ( EX_GGAa + EX_GGAb ) * grid_w_->pointer()[p]; 
//     }
//     return exc;
// }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++ Correlation Functionals +++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Note: This correlation functional depends on B86MGC exchange functional
// with empirical atomic parameters t and u defined below.
// From T. Tsuneda, J. Chem. Phys. 110, 10664 (1999).

double MCPDFTSolver::EC_B88_OP(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB){

   const double beta = 0.0042;
 
   double * rho_ap = RHO_A->pointer();
   double * rho_bp = RHO_B->pointer();
   double * sigma_aap = SIGMA_AA->pointer();
   double * sigma_bbp = SIGMA_BB->pointer();

   double exc = 0.0;
   for (int p = 0; p < phi_points_; p++) {

       double rhoa = rho_ap[p];
       double rhob = rho_bp[p];
       double rhoa_43 = pow( rhoa, 4.0/3.0);
       double rhob_43 = pow( rhob, 4.0/3.0);
       double rhoa_13 = pow( rhoa, 1.0/3.0);
       double rhob_13 = pow( rhob, 1.0/3.0);
       double sigmaaa = sigma_aap[p];
       double sigmabb = sigma_bbp[p];
       double Xa = sqrt(sigmaaa) / rhoa_43;
       double Xb = sqrt(sigmabb) / rhob_43;
       double Xa_2 = Xa * Xa;
       double Xb_2 = Xb * Xb;
      
       auto Ks = [=](double Xs, double Xs_2) -> double{
               
                 double temp3 = 3.0 * pow(3.0/(4.0 * M_PI),1.0/3.0) + 2.0 * (beta * Xs_2) / (1.0 + 6.0 * beta * Xs * asinh(Xs));
                 return temp3; 
       };

       double Ka = Ks(Xa, Xa_2);
       double Kb = Ks(Xb, Xb_2);

       double BETA_ab = 2.3670 * (rhoa_13 * rhob_13 * Ka * Kb) / (rhoa_13 * Ka + rhob_13 * Kb);
       
       exc += -rhoa * rhob * ( (1.5214 * BETA_ab + 0.5764) / ( pow(BETA_ab,4.0) + 1.1284 * pow(BETA_ab,3.0) + 0.3183 * pow(BETA_ab,2.0) ) ) * grid_w_->pointer()[p];

   }
   return exc;
}

// double MCPDFTSolver::EC_B88(){
// 
//    const double Cx = 0.73855876638202240586;
//    const double c = pow(2.0,1.0/3.0) * Cx;
//    const double t = 0.63;
//    const double u = 0.96;
//    const double beta = 0.00375;
//    // const double beta = 0.0042;
//    const double lambda = 0.007;
//  
//    double * rho_ap = rho_a_->pointer();
//    double * rho_bp = rho_b_->pointer();
//    double * sigma_aap = sigma_aa_->pointer();
//    double * sigma_bbp = sigma_bb_->pointer();
//    double * tau_ap = tau_a_->pointer();
//    double * tau_bp = tau_b_->pointer();
// 
//    double exc = 0.0;
//    for (int p = 0; p < phi_points_; p++) {
// 
//        double rhoa = rho_ap[p];
//        double rhob = rho_bp[p];
//        double rho = rhoa + rhob;
//        double rhoa_43 = pow( rhoa, 4.0/3.0);
//        double rhob_43 = pow( rhob, 4.0/3.0);
//        double rhoa_13 = pow( rhoa, 1.0/3.0);
//        double rhob_13 = pow( rhob, 1.0/3.0);
//        double taua = tau_ap[p];
//        double taub = tau_bp[p];
//        double sigmaaa = sigma_aap[p];
//        double sigmabb = sigma_bbp[p];
//        double Xa = sqrt(sigmaaa) / rhoa_43;
//        double Xb = sqrt(sigmabb) / rhob_43;
//        double Xa_2 = Xa * Xa;
//        double Xb_2 = Xb * Xb;
// 
//        auto xyfunc = [=](double rhos_13, double Xs_2) -> double{
//                 
//                      double temp = 0.5 * pow(c * rhos_13 + (beta * Xs_2 * rhos_13) / pow(1.0 + lambda * Xs_2 ,4.0/5.0 ) ,-1.0);
//                      return temp;
//        };
//        
//        auto rfunc = [=](double rhos, double rhos_43, double Xs_2) -> double{
//        
//                 double temp1 =  0.5 * rhos * pow( c * rhos_43 + (beta * Xs_2 * rhos_43)/ pow(1.0 + lambda * Xs_2 ,4.0/5.0) ,-1.0);
//                 return temp1;
//        };
// 
//        auto dfunc = [=](double rhos, double taus, double sigmass) -> double {
// 
//                 double temp2 = taus - 0.25 * sigmass / rhos;
//                 return temp2;
//        };
// 
//        double x = xyfunc(rhoa_13,Xa_2);
//        double y = xyfunc(rhob_13,Xb_2);
//        double q = t * (x + y);
//        double q_2 = q * q;
//        double r_a = rfunc(rhoa, rhoa_43, Xa_2);
//        double r_b = rfunc(rhob, rhob_43, Xb_2);
//        double z_a = 2.0 * u * r_a;
//        double z_b = 2.0 * u * r_b;
//        double d_a = dfunc(rhoa, taua, sigmaaa);
//        double d_b = dfunc(rhob, taub, sigmabb);
//        
//        double f = -0.8 * rhoa * rhob * q_2 * (1.0 - log(1.0 + q) / q);
//        double g_a = -0.01 * rhoa * d_a * pow(z_a,4.0) * ( 1.0 - 2.0 * log(1.0 + 0.5 * z_a) / z_a );
//        double g_b = -0.01 * rhob * d_b * pow(z_b,4.0) * ( 1.0 - 2.0 * log(1.0 + 0.5 * z_b) / z_b );
//       
//        auto Ks = [=](double Xs, double Xs_2) -> double{
//                
//                  double temp3 = 3.0 * pow(3.0/(4.0 * M_PI),1.0/3.0) + 2.0 * (beta * Xs_2) / (1.0 + 6.0 * beta * Xs * asinh(Xs));
//                  return temp3; 
//        };
// 
//        double Ka = Ks(Xa, Xa_2);
//        double Kb = Ks(Xb, Xb_2);
// 
//        exc += ( f + g_a + g_b ) * grid_w_->pointer()[p];
//        double BETA_ab = 2.3670 * (rhoa_13 * rhob_13 * Ka * Kb) / (rhoa_13 * Ka + rhob_13 * Kb);
//        // exc += -rhoa * rhob * ( (1.5214 * BETA_ab + 0.5764) / ( pow(BETA_ab,4.0) + 1.1284 * pow(BETA_ab,3.0) + 0.3183 * pow(BETA_ab,2.0) ) ) * grid_w_->pointer()[p];
// 
//    }
//    return exc;
// }

double MCPDFTSolver::EC_PBE(){

    const double P = 1.0;
    const double T1 = 0.031091;
    const double T2 = 0.015545;
    const double T3 = 0.016887;
    const double U1 = 0.21370;
    const double U2 = 0.20548;
    const double U3 = 0.11125;
    const double V1 = 7.5957;
    const double V2 = 14.1189;
    const double V3 = 10.357;
    const double W1 = 3.5876;
    const double W2 = 6.1977;
    const double W3 = 3.6231;
    const double X1 = 1.6382;
    const double X2 = 3.3662;
    const double X3 = 0.88026;
    const double Y1 = 0.49294;
    const double Y2 = 0.62517;
    const double Y3 = 0.49671;

    const double KSI = 23.266;
    const double PHII = 0.007389;
    const double LAMBDA = 8.723;
    const double UPSILON = 0.472;
    const double c = 1.709921;
    const double t = 0.0716;
    const double v = (16.0/M_PI) * pow(3.0 * M_PI * M_PI ,1.0/3.0);
    const double k = 0.004235;
    const double Z = -0.001667;
    const double lambda = v * k;
   
    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    
    double * sigma_aap = sigma_aa_->pointer();
    double * sigma_abp = sigma_ab_->pointer();
    double * sigma_bbp = sigma_bb_->pointer();

    double * zeta_p = zeta_->pointer();
    double * rs_p = rs_->pointer();

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double zeta = zeta_p[p];
        double sigmaaa = sigma_aap[p];
        double sigmaab = sigma_abp[p];
        double sigmabb = sigma_bbp[p];
        double rs = rs_p[p];
        
        double rho = rhoa + rhob;
        double sigma = sqrt(sigmaaa + sigmabb + 2.0 * sigmaab);

        // build f(zeta) weight factor where f(0) = 0 and f(1) = 1
        double omega =  (pow((1.0 + zeta) ,4.0/3.0) + pow((1.0 - zeta) ,4.0/3.0) - 2.0) / (2.0 * pow(2,1.0/3.0) - 2.0);

        // build spin scaling factor
        double u = 0.5 * ( pow( (1.0 + zeta) ,2.0/3.0 ) + pow( (1.0 - zeta) ,2.0/3.0) );
  
        double d = ( sqrt(sigma) / 12.0 * u ) * pow( pow(3.0,5.0) * M_PI / pow(rho,7.0) ,1.0/6.0);
     
        auto e = [](double r, double T, double U, double V, double W, double X, double Y, int P) -> double {

                 double dum = -2.0 * T * (1.0 + U * r) * log (1.0 + 0.5 / ( T * (V * sqrt(r) + W * r + X * pow(r,3.0/2.0) + Y * pow(r,P+1)))); 
        };
              
        double eps = e(rs,T1,U1,V1,W1,X1,Y1,P) - ( e(rs,T3,U3,V3,W3,X3,Y3,P) * omega * (1.0 - pow(zeta ,4.0)) ) / c 
                   + ( e(rs,T2,U2,V2,W2,X2,Y2,P) - e(rs,T1,U1,V1,W1,X1,Y1,P)) * omega * pow(zeta ,4.0);
                 
        double A = 2.0 * t * pow(lambda ,-1.0) * pow( exp( (-2.0 * t * eps) / (pow(u,3.0) * pow(lambda ,2.0)) ) - 1.0 , -1.0);
       
        double H = 0.5 * pow(u,3.0) * pow(lambda ,2.0) * log(1.0 + ( 2.0 * t * (pow(d ,2.0) + A * pow(d ,4.0)) )
                             / ( lambda * (1.0 + A * pow(d ,2.0) + pow(A,2.0) * pow(d ,4.0)) )) * pow(t ,-1.0);

        auto Q = [=](double sigmass) -> double{

                 double temp = sqrt(sigmass) * pow(2.0, 1.0/3.0) * pow( pow(3.0,5.0) * M_PI ,1.0/6.0) / ( 12.0 * pow(rho,7.0/6.0) );
                 return temp;
        };
    
        auto phi = [=](double r) -> double {

                   double theta = 0.001 * (2.568 + KSI * r + PHII * pow(r,2.0)) / (1.0 + LAMBDA * r + UPSILON * pow(r,2.0) + 10.0 * PHII * pow(r,3.0));
                   double temp = theta - Z; 
                   return temp;
        };
       
        exc += rho * (eps + H) * grid_w_->pointer()[p];
        exc += rho * (d + u + omega + sigma) * grid_w_->pointer()[p];
    }
    return exc;
}


double MCPDFTSolver::EC_VWN3_RPA(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B, std::shared_ptr<Vector> ZETA, std::shared_ptr<Vector> RS){
  
    const double k1 = 0.0310907;
    const double k2 = 0.01554535;
    const double l1 = -0.409286;
    const double l2 = -0.743294;
    const double m1 = 13.0720;
    const double m2 = 20.1231;
    const double n1 = 42.7198;
    const double n2 = 101.578;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();
    
    double * zeta_p = ZETA->pointer();
    double * rs_p = RS->pointer();

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {
        
        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;
        double zeta = zeta_p[p];
        double rs = rs_p[p];
        double x = sqrt(rs);
        
        double y = (9.0/8.0) * pow(1.0 + zeta, 4.0/3.0) + (9.0/8.0) * pow(1.0 - zeta, 4.0/3.0) - (9.0/4.0);
        double z = (4.0 * y) / (9.0 * pow(2.0,1.0/3.0) - 9.0);

        auto X = [](double i, double c, double d) -> double{
                 
                 double temp = pow(i,2.0) + c * i + d;
                 return temp;
        };
        
        auto Q = [](double c, double d) -> double{
        
                 double temp1 = sqrt( 4 * d - pow(c,2.0) );
                 return temp1;
        };

        auto q = [=](double A, double p, double c, double d) -> double{
            
                 double dum1 = A * ( log( pow(x,2.0) / X(x,c,d) ) + 2.0 * c * atan( Q(c,d)/(2.0*x + c) ) * pow(Q(c,d),-1.0)    
                           - c * p * ( log( pow(x-p,2.0) / X(x,c,d) ) + 2.0 * (c + 2.0 * p) * atan( Q(c,d)/(2.0*x + c) ) * pow(Q(c,d),-1.0) ) * pow(X(p,c,d),-1.0) ); 
                 return dum1;
        };

        double Lambda = q(k1, l1, m1, n1);
        double lambda = q(k2, l2, m2, n2);
       
        double e =  Lambda + z * (lambda - Lambda); 
 
        exc += e * rho * grid_w_->pointer()[p]; 
    }
    return exc;
    
}

// double MCPDFTSolver::EC_VWN3_RPA(){
// 
//             // build alpha_c(rs) factor
//             double alphac = Gfunction(rs,Aa_,a1a_,b1a_,b2a_,b3a_,b4a_,pa_);
// 
//             // build ec(rs,zeta) at (rs,0)
//             double ec_rs0 = Gfunction(rs,c0p_,a1p_,b1p_,b2p_,b3p_,b4p_,pe_);
// 
//             // build ec(rs,zeta) at (rs,1)        
//             double ec_rs1 = Gfunction(rs,c0f_,a1f_,b1f_,b2f_,b3f_,b4f_,pe_);

}} // End namespaces
