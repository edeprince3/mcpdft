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

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @ The four major conventions regarding to the principal variables used for                                              @
// @ building exchange-correlation (XC) functionals are as follows:                                                        @
// @                                                                                                                       @
// @ Convention (I)   spin densities rho_a, rho_b and their gradients:                   EXC[rho_a, rho_b, rho_a', rho_b'] @
// @ Convention (II)  total density, spin-magnetization density m and their gradients:   EXC[rho, m, rho', m']             @
// @ Convention (III) total density, spin polarization and their gradients:              EXC[rho, zeta, rho', zeta']       @
// @ Convention (IV)  singlet and triplet charge densities                                                                 @
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

namespace psi{ namespace RDMinoles {

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++ Exchange functionals ++++++++++++++++++
    // wPBE, wB97, Lh-BLYP, wB88                                +
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double MCPDFTSolver::EX_wPBE_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B,
                             std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB) {

    const double OMEGA = options_.get_double("MCPDFT_OMEGA");

    const double A_bar = 0.757211;
    const double B = -0.106364;
    const double C = -0.118649;
    const double D = 0.609650;
    const double E = -0.0477963;

    const double a2 = 0.0159941;
    const double a3 = 0.0852995;
    const double a4 = -0.160368;
    const double a5 = 0.152645;
    const double a6 = -0.0971263;
    const double a7 = 0.0422061;
    const double b1 = 5.33319;
    const double b2 = -12.4780;
    const double b3 = 11.0988;
    const double b4 = -5.11013;
    const double b5 = 1.71468;
    const double b6 = -0.610380;
    const double b7 = 0.307555;
    const double b8 = -0.0770547;
    const double b9 = 0.0334840;

    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();

    double * sigma_aap = SIGMA_AA->pointer();
    double * sigma_bbp = SIGMA_BB->pointer();

    auto kF = [](double RHO) -> double {

              double dum = pow(3.0 * M_PI * M_PI * RHO ,1.0/3.0);
              return dum;
    };

    auto S = [=](double RHO, double SIGMA) -> double {

             double temp = sqrt(SIGMA) / (2.0 * kF(RHO) * RHO);
             return temp;
    };

    auto H = [=](double RHO, double SIGMA) -> double {

             double s1 = S(RHO,SIGMA);
             double s2 = s1 * s1;
             double s3 = s2 * s1;
             double s4 = s3 * s1;
             double s5 = s4 * s1;
             double s6 = s5 * s1;
             double s7 = s6 * s1;
             double s8 = s7 * s1;
             double s9 = s8 * s1;

             double numerator   = a2 * s2 + a3 * s3 + a4 * s4 + a5 * s5 + a6 * s6 + a7 * s7;
             double denominator = 1.0 + b1 * s1 + b2 * s2 + b3 * s3 + b4 * s4 + b5 * s5 + b6 * s6 + b7 * s7 + b8 * s8 + b9 * s9;
             double dum = numerator/denominator;

             return dum;
    };

    auto Nu = [=](double RHO) -> double {

              double dum = OMEGA / kF(RHO);
              return dum;
    };

    auto ZETA = [=](double RHO, double SIGMA) -> double {

             double dum = pow(S(RHO,SIGMA), 2.0) * H(RHO, SIGMA);
             return dum;

    };

    auto ETA = [=](double RHO, double SIGMA) -> double {

             double temp = A_bar + ZETA(RHO, SIGMA);
             return temp;

    };

    auto LAMBDA = [=](double RHO, double SIGMA) -> double {

             double temp = D + ZETA(RHO, SIGMA);
             return temp;

    };

    auto Chi = [=](double RHO, double SIGMA) -> double {

              double dum = Nu(RHO) / sqrt(LAMBDA(RHO,SIGMA) + pow(Nu(RHO),2.0));
              return dum;
    };

    auto BG_LAMBDA = [=](double RHO, double SIGMA) -> double {

             double temp = LAMBDA(RHO,SIGMA) / (1.0 - Chi(RHO,SIGMA));
             return temp;

    };

    auto F_bar = [=](double RHO, double SIGMA) -> double {

                 double s0 = 2.0;
                 double s1 = S(RHO,SIGMA);
                 double s2 = s1 * s1;
                 double ze_v = ZETA(RHO,SIGMA);

                 double dum = 1.0 - (1.0/(27.0*C)) * (s2/(1.0 + (s2/pow(s0,2.0)))) - (1.0/(2.0*C)) * ze_v;

                 return dum;

    };

    auto G_bar = [=](double RHO, double SIGMA) -> double {

                 double ze_v   = ZETA(RHO,SIGMA);
                 double et_v   = ETA(RHO,SIGMA);
                 double lam_v  = LAMBDA(RHO,SIGMA);
                 double lam_v2 = lam_v * lam_v;
                 double lam_v3 = lam_v2 * lam_v;
                 double lam_v72 = pow(lam_v,7.0/2.0);
                 double sq_ze = sqrt(ze_v);
                 double sq_et = sqrt(et_v);

                 double dum = -(2.0/5.0) * C * F_bar(RHO,SIGMA) * lam_v - (4.0/15.0) * B * lam_v2 - (6.0/5.0) * A_bar * lam_v3
                            - (4.0/5.0) * sqrt(M_PI) * lam_v72 - (12.0/5.0) * lam_v72 * (sq_ze - sq_et);

                 dum *= (1.0/E);

                 return dum;

    };

    auto eX = [=](double RHO) -> double {

              double temp = -(3.0 * kF(RHO)) / (4.0 * M_PI);
              return temp;
    };

    auto FX = [=](double RHO, double SIGMA) -> double {

              double chi_v  = Chi(RHO,SIGMA);
              double chi_v2 = chi_v * chi_v;
              double chi_v3 = chi_v2 * chi_v;
              double chi_v5 = chi_v3 * chi_v2;
              double nu_v   = Nu(RHO);
              double ze_v   = ZETA(RHO,SIGMA);
              double et_v   = ETA(RHO,SIGMA);
              double bg_lam_v  = BG_LAMBDA(RHO,SIGMA);
              double bg_lam_v2 = bg_lam_v * bg_lam_v;
              double bg_lam_v3 = bg_lam_v2 * bg_lam_v;
              double lam_v  = LAMBDA(RHO,SIGMA);
              double lam_v2 = lam_v * lam_v;
              double lam_v3 = lam_v2 * lam_v;
              double sq_ze_nu = sqrt(ze_v  + nu_v*nu_v);
              double sq_et_nu = sqrt(et_v  + nu_v*nu_v);
              double sq_la_nu = sqrt(lam_v + nu_v*nu_v);

              double temp = A_bar * ( ze_v/( (sq_ze_nu + sq_et_nu) * (sq_ze_nu + nu_v) ) + et_v/( (sq_ze_nu + sq_et_nu) * (sq_et_nu + nu_v) ) )
                          - (4.0/9.0) * (B/bg_lam_v) - (4.0/9.0) * (C * F_bar(RHO,SIGMA) / bg_lam_v2) * (1.0 + (1.0/2.0) * chi_v )
                          - (8.0/9.0) * (E * G_bar(RHO,SIGMA) / bg_lam_v3) * (1.0 + (9.0/8.0) * chi_v + (3.0/8.0) * chi_v2)
                          + 2.0 * ze_v * log( 1.0 - (lam_v - ze_v) / ( (nu_v + sq_la_nu) * (sq_la_nu + sq_ze_nu) ) )
                          - 2.0 * et_v * log( 1.0 - (lam_v - et_v) / ( (nu_v + sq_la_nu) * (sq_la_nu + sq_et_nu) ) );

              return temp;
    };

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;

        double sigmaaa = sigma_aap[p];
        double sigmabb = sigma_bbp[p];
        double sigma = 0.0;

        double tol = 1.0e-20;
        if ( rho > tol ) {
           if ( rhoa < tol ){

              rho = rhob;
              sigmabb = std::max(0.0,sigmabb);
              sigma = sigmabb;

              double zk = eX(2.0 * rho) * FX(2.0 * rho, 4.0 * sigma);
              exc += rho * zk * grid_w_->pointer()[p];

           }else if ( rhob < tol ){

                    rho = rhoa;
                    sigmaaa = std::max(0.0,sigmaaa);
                    sigma = sigmaaa;

                    double zk = eX(2.0 * rho) * FX(2.0 * rho, 4.0 * sigma);
                    exc += rho * zk * grid_w_->pointer()[p];
           }else {

                 double zka = rhoa * eX(2.0 * rhoa) * FX(2.0 * rhoa, 4.0 * sigmaaa);
                 double zkb = rhob * eX(2.0 * rhob) * FX(2.0 * rhob, 4.0 * sigmabb);
                 double zk = zka + zkb;
                 exc += zk * grid_w_->pointer()[p];
           }
        }else{
                //double zk = 0.0;        
                exc += 0.0;
             }
    }
    return exc;
}

double MCPDFTSolver::Lh_EX_B88_I(std::shared_ptr<Vector> RHO_A, std::shared_ptr<Vector> RHO_B,
                                 std::shared_ptr<Vector> SIGMA_AA, std::shared_ptr<Vector> SIGMA_BB){

    double tol = 1.0e-20;

    const double Cx = 0.73855876638202240586;
    const double beta = 0.0042;
    const double c = pow(2.0,1.0/3.0) * Cx;

    double * ex_exact_p = ex_exact_->pointer();
    double * lmf_p  = lmf_->pointer();
    double * rho_ap = RHO_A->pointer();
    double * rho_bp = RHO_B->pointer();
    double * sigma_aap = SIGMA_AA->pointer();
    double * sigma_bbp = SIGMA_BB->pointer();

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;
        double sigmaaa = sigma_aap[p];
        double sigmabb = sigma_bbp[p];
        double rhoa_43 = 0.0;
        double rhob_43 = 0.0;
        double Xa   = 0.0;
        double Xb   = 0.0;
        double Xa_2 = 0.0;
        double Xb_2 = 0.0;

        if ( rho > tol ) {
           if ( rhoa < tol ){

              rhoa = 0.0;
              sigmabb = std::max(0.0,sigma_bbp[p]);
              rhob_43 = pow(rhob, 4.0/3.0);
              Xb = sqrt(sigma_bbp[p]) / rhob_43;
              Xb_2 = Xb * Xb;

           }else if ( rhob < tol ){

                    rhob = 0.0;
                    sigmaaa = std::max(0.0,sigma_aap[p]);
                    rhoa_43 = pow( rhoa, 4.0/3.0);
                    Xa = sqrt(sigma_aap[p]) / rhoa_43;
                    Xa_2 = Xa * Xa;
           }else{

                sigmaaa = std::max(0.0,sigma_aap[p]);
                sigmabb = std::max(0.0,sigma_bbp[p]);
                rhoa_43 = pow( rhoa, 4.0/3.0);
                rhob_43 = pow( rhob, 4.0/3.0);
                Xa = sqrt(sigma_aap[p]) / rhoa_43;
                Xb = sqrt(sigma_bbp[p]) / rhob_43;
                Xa_2 = Xa * Xa;
                Xb_2 = Xb * Xb;

           }
           exc += -rhoa_43 * ( c + (beta * Xa_2) / (1.0 + 6.0 * beta * Xa * asinh(Xa)) ) * (1.0 - lmf_p[p]) * grid_w_->pointer()[p];
           exc += -rhob_43 * ( c + (beta * Xb_2) / (1.0 + 6.0 * beta * Xb * asinh(Xb)) ) * (1.0 - lmf_p[p]) * grid_w_->pointer()[p];
	   exc += pow(rho,1.0/3.0) * ex_exact_p[p] * lmf_p[p] * grid_w_->pointer()[p];// TODO:Check the spin-polarized version as well

        }else{

             exc += 0.0;
        }
    }
    return exc;
}

double DFTSolver::EX_wB88_I(){
    
    const double OMEGA = options_.get_double("RS_OMEGA");
    
    const double A_bar = 0.757211;
    const double B = -0.106364;
    const double C = -0.118649;
    const double D = 0.609650;
    const double E = -0.0477963;

    const double a2 = 0.0253933;
    const double a3 = -0.0673075;
    const double a4 = 0.0891476;
    const double a5 = -0.0454168;
    const double a6 = -0.0076581;
    const double a7 = 0.0142506;
    const double b1 = -2.65060;
    const double b2 = 3.91108;
    const double b3 = -3.31509;
    const double b4 = 1.54485;
    const double b5 = -0.198386;
    const double b6 = -0.136112;
    const double b7 = 0.0647862;
    const double b8 = 0.0159586;
    const double b9 = -0.000245066;

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
   
    double * sigma_aap = sigma_aa_->pointer();
    double * sigma_bbp = sigma_bb_->pointer();
    
    auto kF = [](double RHO) -> double {

              double dum = pow(3.0 * M_PI * M_PI * RHO ,1.0/3.0);
              return dum;
    };

    auto SS = [=](double RHO, double SIGMA) -> double {

             double temp = sqrt(SIGMA) / (2.0 * kF(RHO) * RHO); 
             return temp;
    };

    auto S = [=](double RHO, double SIGMA) -> double {

             double xi = 1.0 / (exp(20.0) - 1.0);
             double temp = -log( ( exp(-SS(RHO,SIGMA)) + xi ) / (1.0 + xi) );
             return temp;
    };

    auto H = [=](double RHO, double SIGMA) -> double {

             double s1 = S(RHO,SIGMA);
             double s2 = s1 * s1;
             double s3 = s2 * s1;
             double s4 = s3 * s1;
             double s5 = s4 * s1;
             double s6 = s5 * s1;
             double s7 = s6 * s1;
             double s8 = s7 * s1;
             double s9 = s8 * s1;
  
             double numerator   = a2 * s2 + a3 * s3 + a4 * s4 + a5 * s5 + a6 * s6 + a7 * s7;
             double denominator = 1.0 + b1 * s1 + b2 * s2 + b3 * s3 + b4 * s4 + b5 * s5 + b6 * s6 + b7 * s7 + b8 * s8 + b9 * s9;
             double dum = numerator/denominator;

             return dum;

    };

    auto Nu = [=](double RHO) -> double {

              double dum = OMEGA / kF(RHO);
              return dum;
    };
 
    auto ZETA = [=](double RHO, double SIGMA) -> double {
  
             double dum = pow(S(RHO,SIGMA), 2.0) * H(RHO, SIGMA);
             return dum;

    };

    auto ETA = [=](double RHO, double SIGMA) -> double {
  
             double temp = A_bar + ZETA(RHO, SIGMA);
             return temp;

    };

    auto LAMBDA = [=](double RHO, double SIGMA) -> double {
  
             double temp = D + ZETA(RHO, SIGMA);
             return temp;

    };

    auto Chi = [=](double RHO, double SIGMA) -> double {

              double dum = Nu(RHO) / sqrt(LAMBDA(RHO,SIGMA) + pow(Nu(RHO),2.0));
              return dum;
    };

    auto BG_LAMBDA = [=](double RHO, double SIGMA) -> double {
  
             double temp = LAMBDA(RHO,SIGMA) / (1.0 - Chi(RHO,SIGMA));
             return temp;

    };

    auto F_bar = [=](double RHO, double SIGMA) -> double {

                 double s0 = 2.0;
                 double s1 = S(RHO,SIGMA);
                 double s2 = s1 * s1;
                 double ze_v = ZETA(RHO,SIGMA);
  
                 double dum = 1.0 - (1.0/(27.0*C)) * (s2/(1.0 + (s2/pow(s0,2.0)))) - (1.0/(2.0*C)) * ze_v;

                 return dum;

    };
    
    auto G_bar = [=](double RHO, double SIGMA) -> double {

                 double ze_v   = ZETA(RHO,SIGMA);
                 double et_v   = ETA(RHO,SIGMA);
                 double lam_v  = LAMBDA(RHO,SIGMA);
                 double lam_v2 = lam_v * lam_v;
                 double lam_v3 = lam_v2 * lam_v;
                 double lam_v72 = pow(lam_v,7.0/2.0);
                 double sq_ze = sqrt(ze_v);
                 double sq_et = sqrt(et_v);
  
                 double dum = -(2.0/5.0) * C * F_bar(RHO,SIGMA) * lam_v - (4.0/15.0) * B * lam_v2 - (6.0/5.0) * A_bar * lam_v3
                            - (4.0/5.0) * sqrt(M_PI) * lam_v72 - (12.0/5.0) * lam_v72 * (sq_ze - sq_et);

                 dum *= (1.0/E);

                 return dum;

    };

    auto eX = [=](double RHO) -> double {
    
              double temp = -(3.0 * kF(RHO)) / (4.0 * M_PI);
              return temp;
    };

    auto FX = [=](double RHO, double SIGMA) -> double {

              double chi_v  = Chi(RHO,SIGMA);
              double chi_v2 = chi_v * chi_v;
              double chi_v3 = chi_v2 * chi_v;
              double chi_v5 = chi_v3 * chi_v2;
              double nu_v   = Nu(RHO);
              double ze_v   = ZETA(RHO,SIGMA);
              double et_v   = ETA(RHO,SIGMA);
              double bg_lam_v  = BG_LAMBDA(RHO,SIGMA);
              double bg_lam_v2 = bg_lam_v * bg_lam_v;
              double bg_lam_v3 = bg_lam_v2 * bg_lam_v;
              double lam_v  = LAMBDA(RHO,SIGMA);
              double lam_v2 = lam_v * lam_v;
              double lam_v3 = lam_v2 * lam_v;
              double sq_ze_nu = sqrt(ze_v  + nu_v*nu_v);
              double sq_et_nu = sqrt(et_v  + nu_v*nu_v);
              double sq_la_nu = sqrt(lam_v + nu_v*nu_v);

              double temp = A_bar * ( ze_v/( (sq_ze_nu + sq_et_nu) * (sq_ze_nu + nu_v) ) + et_v/( (sq_ze_nu + sq_et_nu) * (sq_et_nu + nu_v) ) )
                          - (4.0/9.0) * (B/bg_lam_v) - (4.0/9.0) * (C * F_bar(RHO,SIGMA) / bg_lam_v2) * (1.0 + (1.0/2.0) * chi_v )
                          - (8.0/9.0) * (E * G_bar(RHO,SIGMA) / bg_lam_v3) * (1.0 + (9.0/8.0) * chi_v + (3.0/8.0) * chi_v2)
                          + 2.0 * ze_v * log( 1.0 - (lam_v - ze_v) / ( (nu_v + sq_la_nu) * (sq_la_nu + sq_ze_nu) ) )
                          - 2.0 * et_v * log( 1.0 - (lam_v - et_v) / ( (nu_v + sq_la_nu) * (sq_la_nu + sq_et_nu) ) );

              return temp;
    };

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {
   
        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;

        double sigmaaa = sigma_aap[p];
        double sigmabb = sigma_bbp[p];
        double sigma = 0.0;
        
        double tol = 1.0e-20;     
        if ( rho > tol ) {
           if ( rhoa < tol ){
              
              rho = rhob;
              sigmabb = std::max(0.0,sigmabb);
              sigma = sigmabb;   
              
              double zk = eX(2.0 * rho) * FX(2.0 * rho, 4.0 * sigma);
              exc += rho * zk * grid_w_->pointer()[p];

           }else if ( rhob < tol ){
                                  
                    rho = rhoa;
                    sigmaaa = std::max(0.0,sigmaaa);
                    sigma = sigmaaa;   
                    
                    double zk = eX(2.0 * rho) * FX(2.0 * rho, 4.0 * sigma);
                    exc += rho * zk * grid_w_->pointer()[p];
           }else {
           
                 double zka = rhoa * eX(2.0 * rhoa) * FX(2.0 * rhoa, 4.0 * sigmaaa);
                 double zkb = rhob * eX(2.0 * rhob) * FX(2.0 * rhob, 4.0 * sigmabb);
                 double zk = zka + zkb;
                 exc += zk * grid_w_->pointer()[p];
           }
        }else{
                //double zk = 0.0;        
                exc += 0.0;
             }
    }
    return exc;
}

double MCDFTSolver::EX_SR_B97() {
    const double cx_0 = 1.00000;
    const double cx_1 = 1.13116;
    const double cx_2 = -2.74915;
    const double cx_3 = 12.0900;
    const double cx_4 = -5.71642;
    const double alpha = (2.0/3.0);      // Slater value
    const double Cx = (9.0/8.0) * alpha * pow(3.0/M_PI,1.0/3.0);
    const double OMEGA = 0.4;
    const double GAMMA_X = 0.004;
    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * sigma_aap = sigma_aa_->pointer();
    double * sigma_bbp = sigma_bb_->pointer();

    double ex_SR_B97 = 0.0;
    for (int p = 0; p < phi_points_; p++) {
        double s2_a = sigma_aap[p] / pow(rho_ap[p], 8.0/3.0);
        double s2_b = sigma_bbp[p] / pow(rho_bp[p], 8.0/3.0);
        double ux_a = (GAMMA_X * s2_a) / (1.0 + GAMMA_X * s2_a);
        double ux_b = (GAMMA_X * s2_b) / (1.0 + GAMMA_X * s2_b);
        double gx_a = cx_0 + cx_1 * ux_a + cx_2 * pow(ux_a,2.0) + cx_3 * pow(ux_a,3.0) + cx_4 * pow(ux_a,4.0);
        double gx_b = cx_0 + cx_1 * ux_b + cx_2 * pow(ux_b,2.0) + cx_3 * pow(ux_b,3.0) + cx_4 * pow(ux_b,4.0);
        double kf_a = pow(6.0 * M_PI * M_PI* rho_ap[p], 1.0/3.0);
        double kf_b = pow(6.0 * M_PI * M_PI* rho_bp[p], 1.0/3.0);
        double a_a = OMEGA / (2.0 * kf_a);
        double a_b = OMEGA / (2.0 * kf_b);
        double Fa_a = 1.0 - (8.0/3.0) * a_a * ( sqrt(M_PI) * erf(1.0 / (2.0*a_a)) - 3.0 * a_a + 4.0 * pow(a_a,3.0)
                    + (2.0 * a_a - 4.0 * pow(a_a, 3.0)) * exp(-(1.0/(4.0 * a_a * a_a))) );
        double Fa_b = 1.0 - (8.0/3.0) * a_b * ( sqrt(M_PI) * erf(1.0 / (2.0*a_b)) - 3.0 * a_b + 4.0 * pow(a_b,3.0)
                    + (2.0 * a_b - 4.0 * pow(a_b, 3.0)) * exp(-(1.0/(4.0 * a_b * a_b))) );
        double exa_SR_LSDA = pow(2.0,1.0/3.0) * Cx * pow( rho_ap[p], 4.0/3.0) * Fa_a;
        double exb_SR_LSDA = pow(2.0,1.0/3.0) * Cx * pow( rho_bp[p], 4.0/3.0) * Fa_b;
        double exa_SR_B97 = exa_SR_LSDA * gx_a * grid_w_->pointer()[p];
        double exb_SR_B97 = exb_SR_LSDA * gx_b * grid_w_->pointer()[p];
        ex_SR_B97 += exa_SR_B97 + exb_SR_B97; 
    }
    return ex_SR_B97;

}
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++ Correlation Functionals +++++++++++++++++
    // WB97 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double MCPDFTSolver::EC_B97(){

    const double pa = 1.0;
    const double Aa = 0.0168869;
    const double a1a = 0.11125;
    const double b1a = 10.357;
    const double b2a = 3.6231;
    const double b3a = 0.88026;
    const double b4a = 0.49671;
    const double pe = 1.0;
    const double c0p = 0.0310907;
    const double a1p = 0.21370;
    const double b1p = 7.5957;
    const double b2p = 3.5876;
    const double b3p = 1.6382;
    const double b4p = 0.49294;
    const double c0f = 0.01554535;
    const double a1f = 0.20548;
    const double b1f = 14.1189;
    const double b2f = 6.1977;
    const double b3f = 3.3662;
    const double b4f = 0.62517;
    const double d2Fz = 1.7099209341613656173;

    const double c0_ss = 1.00000;
    const double c0_ab = 1.00000;
    const double c1_ss = -2.55352;
    const double c1_ab = 3.99051;
    const double c2_ss = 11.8926;
    const double c2_ab = -17.0066;
    const double c3_ss = -26.9452;
    const double c3_ab = 1.07292;
    const double c4_ss = 17.0927;
    const double c4_ab = 8.88211;
    const double GAMMA_SS = 0.2;
    const double GAMMA_AB = 0.006;

    double tol = 1.0e-20;

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * sigma_aap = sigma_aa_->pointer();
    double * sigma_bbp = sigma_bb_->pointer();

    auto Fz = [](double ZETA) -> double {

              double dum = (pow((1.0 + ZETA) ,4.0/3.0) + pow((1.0 - ZETA) ,4.0/3.0) - 2.0) / (2.0 * pow(2.0,1.0/3.0) - 2.0);
              return dum;
    };

    auto G = [](double r, double T, double a1, double b1, double b2, double b3, double b4, double p) -> double {

             double dum = -2.0 * T * (1.0 + a1 * r) * log(1.0 + 0.5 * pow(T * (b1 * sqrt(r) + b2 * r + b3 * pow(r,3.0/2.0) + b4 * pow(r, p+1.0)) ,-1.0));
             return dum;

    };

    auto Ac = [=](double r) -> double {

              double temp = -G(r,Aa,a1a,b1a,b2a,b3a,b4a,pa);
              return temp;
    };

    auto EcP = [=](double r) -> double {

               double dum = G(r,c0p,a1p,b1p,b2p,b3p,b4p,pe);
               return dum;
    };

    auto EcF = [=](double r) -> double {

               double dumm = G(r,c0f,a1f,b1f,b2f,b3f,b4f,pe);
               return dumm;
    };

    auto Ec = [=](double r, double ZETA) -> double {

              double dum = EcP(r) + ( Ac(r) * Fz(ZETA) * (1.0 - pow(ZETA ,4.0)) ) / d2Fz + ( EcF(r) - EcP(r) ) * Fz(ZETA) * pow(ZETA ,4.0);
              return dum;
    };

    std::shared_ptr<Vector> ec_LSDA (new Vector(phi_points_) );
    std::shared_ptr<Vector> ec_ab   (new Vector(phi_points_) );
    std::shared_ptr<Vector> ec_aa   (new Vector(phi_points_) );
    std::shared_ptr<Vector> ec_bb   (new Vector(phi_points_) );

    ec_LSDA->zero();
    ec_ab->zero();
    ec_aa->zero();
    ec_bb->zero();

    double * ec_LSDAp = ec_LSDA->pointer();
    double * ec_abp = ec_ab->pointer();
    double * ec_aap = ec_aa->pointer();
    double * ec_bbp = ec_bb->pointer();

    double exc = 0.0;
    for (int p = 0; p < phi_points_; p++) {

        double rhoa = rho_ap[p];
        double rhob = rho_bp[p];
        double rho = rhoa + rhob;
        double zeta = (rhoa - rhob) / rho;
        double rs =  pow( 3.0 / ( 4.0 * M_PI * rho) , 1.0/3.0 );

        if ( rho > tol ) {
           if ( rhoa < tol ){

              rho = rhob;
              zeta = 1.0;
              ec_bbp[p] = Ec(rs,zeta) * rho;
              ec_aap[p] = 0.0;

           }else if ( rhob < tol ){

                    rho = rhoa;
                    zeta = 1.0;
                    ec_aap[p] = Ec(rs,zeta) * rho;
                    ec_bbp[p] = 0.0;

           }
           ec_LSDAp[p] = Ec(rs,zeta) * rho;
        }else{
             ec_LSDAp[p] =  0.0;
             ec_aap[p] =  0.0;
             ec_bbp[p] =  0.0;
        }
        ec_abp[p] = ec_LSDAp[p] - ec_aap[p] - ec_bbp[p];
        double s2_a = sigma_aap[p] / pow(rho_ap[p], 8.0/3.0);
        double s2_b = sigma_bbp[p] / pow(rho_bp[p], 8.0/3.0);
        double s2_av = 0.5 * (s2_a + s2_b);
        double uc_aa = (GAMMA_SS * s2_a) / (1.0 + GAMMA_SS * s2_a);
        double uc_bb = (GAMMA_SS * s2_b) / (1.0 + GAMMA_SS * s2_b);
        double uc_ab = (GAMMA_AB * s2_av) / (1.0 + GAMMA_AB * s2_av);
        double gc_aa = c0_ss + c1_ss * uc_aa + c2_ss * pow(uc_aa,2.0) + c3_ss * pow(uc_aa,3.0) + c4_ss * pow(uc_aa,4.0);
        double gc_bb = c0_ss + c1_ss * uc_bb + c2_ss * pow(uc_bb,2.0) + c3_ss * pow(uc_bb,3.0) + c4_ss * pow(uc_bb,4.0);
        double gc_ab = c0_ab + c1_ab * uc_ab + c2_ab * pow(uc_ab,2.0) + c3_ab * pow(uc_ab,3.0) + c4_ab * pow(uc_ab,4.0);
        
        double EC_AA = ec_aap[p] * grid_w_->pointer()[p];
        double EC_BB = ec_bbp[p] * grid_w_->pointer()[p];
        double EC_AB = ec_abp[p] * grid_w_->pointer()[p];

        exc += EC_AA + EC_BB + EC_AB;
    }
    return exc;
}

}} // End namespaces
