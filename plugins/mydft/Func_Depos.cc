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

#include "dft.h"

namespace psi{ namespace mydft {

    //###########################################################
    //# The exchange and correlation functional implementations #
    //###########################################################

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++ Exchange functionals ++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double DFTSolver::EX_LDA(std::shared_ptr<Vector> rho_a, std::shared_ptr<Vector> rho_b){

    vrho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_a2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * rho_ap = rho_a->pointer();
    double * rho_bp = rho_b->pointer();

    double * vrho_ap = vrho_a_->pointer();
    double * v2rho_a2p = v2rho_a2_->pointer();

    const double Cx = 0.73855876638202240586;
    double exc = 0.0;

    if (deriv_ == 0) {

       for (int p = 0; p < phi_points_; p++) {

           exc += -Cx * pow( (rho_ap[p] + rho_bp[p]), 4.0/3.0) * grid_w_->pointer()[p];
       }
    }else if (deriv_ == 1) {

       for (int p = 0; p < phi_points_; p++) {

           double rho13 = pow( (rho_ap[p] + rho_bp[p]), 1.0/3.0) ;

           exc += -Cx * rho13 * ( rho_ap[p] + rho_bp[p] )  * grid_w_->pointer()[p];

           vrho_ap[p] = -0.98474502184269654114 * rho13;
       }
    }else if (deriv_ == 2) {

       for (int p = 0; p < phi_points_; p++) {

           double rho13 = pow( (rho_ap[p] + rho_bp[p]), 1.0/3.0) ;
           double rho_23 = pow( rho13, -2.0);

           exc += -Cx * rho13 * ( rho_ap[p] + rho_bp[p] ) * grid_w_->pointer()[p];

           vrho_ap[p] = -0.98474502184269654114 * rho13;
           v2rho_a2p[p] = -0.65649668122846436077 * rho_23;

       }
    }else {

       throw PsiException("The differentiation order should be 0, 1 or 2!",__FILE__,__LINE__);
    }
    return exc;
}

double DFTSolver::EX_LSDA(std::shared_ptr<Vector> rho_a, std::shared_ptr<Vector> rho_b, std::shared_ptr<Vector> zeta){

    vrho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    vrho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_a2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_b2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * rho_ap = rho_a->pointer();
    double * rho_bp = rho_b->pointer();

    double * vrho_ap = vrho_a_->pointer();
    double * vrho_bp = vrho_b_->pointer();
    double * v2rho_a2p = v2rho_a2_->pointer();
    double * v2rho_abp = v2rho_ab_->pointer();
    double * v2rho_b2p = v2rho_b2_->pointer();


    const double Cx = 0.73855876638202240586;
    double exc = 0.0;

    if (deriv_ == 0) {

       for (int p = 0; p < phi_points_; p++) {

           exc += -pow(2.0,1.0/3.0) * Cx * ( pow( rho_ap[p], 4.0/3.0) + pow( rho_bp[p], 4.0/3.0) ) * grid_w_->pointer()[p];
       }
    }else if (deriv_ == 1) {

       for (int p = 0; p < phi_points_; p++) {

           double rho13a = pow( rho_ap[p], 1.0/3.0);
           double rho13b = pow( rho_bp[p], 1.0/3.0);

           exc += -pow(2.0,1.0/3.0) * Cx * ( pow( rho_ap[p], 4.0/3.0) + pow( rho_bp[p], 4.0/3.0) ) * grid_w_->pointer()[p];

           vrho_ap[p] = -1.2407009817988000333 * rho13a;
           vrho_bp[p] = -1.2407009817988000333 * rho13b;
       }
    }else if (deriv_ == 2) {

       for (int p = 0; p < phi_points_; p++) {

           double rho13a = pow( rho_ap[p], 1.0/3.0);
           double rho13b = pow( rho_bp[p], 1.0/3.0);
           double rho_23a = pow( rho13a, -2.0);
           double rho_23b = pow( rho13b, -2.0);

           exc += -pow(2.0,1.0/3.0) * Cx * ( pow( rho_ap[p], 4.0/3.0) + pow( rho_bp[p], 4.0/3.0) ) * grid_w_->pointer()[p];

           vrho_ap[p] = -1.2407009817988000333 * rho13a;
           vrho_bp[p] = -1.2407009817988000333 * rho13b;
           v2rho_a2p[p] = -0.41356699393293334443 * rho_23a;
           v2rho_b2p[p] = -0.41356699393293334443 * rho_23b;
           v2rho_abp[p] = 0.0;
       }
    }else {

       throw PsiException("The differentiation order should be 0, 1 or 2!",__FILE__,__LINE__);
    }
    return exc;
}

double DFTSolver::EX_B88(){
    
    vrho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    vrho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_a2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    vsigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    vsigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    vsigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    
    v2rho_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_b2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * sigma_aap = sigma_aa_->pointer();
    double * sigma_abp = sigma_ab_->pointer();
    double * sigma_bbp = sigma_bb_->pointer();
    
    double * vrho_ap = vrho_a_->pointer();
    double * vrho_bp = vrho_b_->pointer();
    double * vsigma_aap = vsigma_aa_->pointer();
    double * vsigma_abp = vsigma_ab_->pointer();
    double * vsigma_bbp = vsigma_bb_->pointer();

    double * v2rho_a2p = v2rho_a2_->pointer();
    double * v2rho_abp = v2rho_ab_->pointer();
    double * v2rho_b2p = v2rho_b2_->pointer();

    double ex_LSDA = DFTSolver::EX_LSDA(rho_a_, rho_b_, zeta_);
    double ex_b88 = 0.0;
    const double Cx = 0.73855876638202240586; 
    const double beta = 0.0042;
    double exc = 0.0;

    if (deriv_ == 0) {
       
       for (int p = 0; p < phi_points_; p++) {
   
           double sig = sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p]; 
           double rho = rho_ap[p] + rho_bp[p];
           double rhoa = rho_ap[p];
           double rhob = rho_bp[p];
           double rhoa_13 = pow( rhoa, 1.0/3.0); 
           double rhob_13 = pow( rhob, 1.0/3.0); 
           double rhoa_23 = pow( rhoa, 2.0/3.0); 
           double rhob_23 = pow( rhob, 2.0/3.0); 
           double rhoa_43 = pow( rhoa, 4.0/3.0); 
           double rhob_43 = pow( rhob, 4.0/3.0); 
           double rhoa_2  = pow( rhoa, 2.0); 
           double rhob_2  = pow( rhob, 2.0); 
           double Xa = sqrt(sigma_aap[p]) / rhoa_43;
           double Xb = sqrt(sigma_bbp[p]) / rhob_43;
  
           exc += -pow(2.0,1.0/3.0) * Cx * ( pow( rho_ap[p], 4.0/3.0) + pow( rho_bp[p], 4.0/3.0) ) * grid_w_->pointer()[p]; 
       }
    }else if (deriv_ == 1) {
    
       for (int p = 0; p < phi_points_; p++) {

           double sig = sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p]; 
           double rho = rho_ap[p] + rho_bp[p];
           double rhoa = rho_ap[p];
           double rhob = rho_bp[p];
           double rhoa_13 = pow( rhoa, 1.0/3.0); 
           double rhob_13 = pow( rhob, 1.0/3.0); 
           double rhoa_23 = pow( rhoa, 2.0/3.0); 
           double rhob_23 = pow( rhob, 2.0/3.0); 
           double rhoa_43 = pow( rhoa, 4.0/3.0); 
           double rhob_43 = pow( rhob, 4.0/3.0); 
           double rhoa_2  = pow( rhoa, 2.0); 
           double rhob_2  = pow( rhob, 2.0); 
           double rhoa_m43 = pow( rhoa, -4.0/3.0); 
           double rhob_m43 = pow( rhob, -4.0/3.0); 
           double rhoa_m73 = pow( rhoa, -7.0/3.0); 
           double rhob_m73 = pow( rhob, -7.0/3.0); 
           double rhoa_m83 = pow( rhoa, -8.0/3.0); 
           double rhob_m83 = pow( rhob, -8.0/3.0); 
          
           double grada = sqrt(sigma_aap[p]);
           double gradb = sqrt(sigma_bbp[p]);
           double Xa = grada / rhoa_43;
           double Xb = gradb / rhob_43;
           double Ta = sigma_aap[p] / rhoa_43;
           double Tb = sigma_bbp[p] / rhob_43;
           double fa = log(Xa + sqrt(1.0 + pow(Xa,2.0) ) );
           double fb = log(Xb + sqrt(1.0 + pow(Xb,2.0) ) );
           double ga = 1.0 + 0.0252 * Xa * fa;
           double gb = 1.0 + 0.0252 * Xb * fb;
           double ha = sqrt( 1.0 + sigma_aap[p] * rhoa_m83 );
           double hb = sqrt( 1.0 + sigma_bbp[p] * rhob_m83 );
           double inv_ga = pow( ga, -1.0);
           double inv_gb = pow( gb, -1.0);
           double inv2_ga = pow( ga, -2.0);
           double inv2_gb = pow( gb, -2.0);
           double inv_ha = pow( ha, -1.0);      
           double inv_hb = pow( hb, -1.0);      
     
           exc += -pow(2.0,1.0/3.0) * Cx * ( pow( rho_ap[p], 4.0/3.0) + pow( rho_bp[p], 4.0/3.0) ) * grid_w_->pointer()[p]; 
           exc += -beta * ( Ta * inv_ga + Tb * inv_gb ) * grid_w_->pointer()[p]; 
           
           vrho_ap[p] = -1.2407009817988000333 * rhoa_13 + (0.0056 * rhoa_m73 * sigma_aap[p] * inv_ga) + (0.0042 * Ta * inv2_ga *
                        ( -0.0336 * grada * rhoa_m73 * fa -0.0336 * sigma_aap[p] / rhoa_23 / rhoa_2 / rhoa * inv_ha ));
           vrho_bp[p] = -1.2407009817988000333 * rhob_13 + (0.0056 * rhob_m73 * sigma_bbp[p] * inv_gb) + (0.0042 * Tb * inv2_gb *
                        ( -0.0336 * gradb * rhob_m73 * fb -0.0336 * sigma_bbp[p] / rhob_23 / rhob_2 / rhob * inv_hb ));
           vsigma_aap[p] = -0.0042 * rhoa_m43 * inv_ga + 0.0042 * Ta * inv2_ga * ( 0.0126 / grada * rhoa_m43 * fa + 0.0126 * rhoa_m83 * inv_ha );
           vsigma_bbp[p] = -0.0042 * rhob_m43 * inv_gb + 0.0042 * Tb * inv2_gb * ( 0.0126 / gradb * rhob_m43 * fb + 0.0126 * rhob_m83 * inv_hb );
           vsigma_abp[p] = 0.0;
       }
    }else if (deriv_ == 2) {

       for (int p = 0; p < phi_points_; p++) {

           double sig = sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p]; 
           double rho = rho_ap[p] + rho_bp[p];
           double rhoa = rho_ap[p];
           double rhob = rho_bp[p];
           double rhoa_13 = pow( rhoa, 1.0/3.0); 
           double rhob_13 = pow( rhob, 1.0/3.0); 
           double rhoa_23 = pow( rhoa, 2.0/3.0); 
           double rhob_23 = pow( rhob, 2.0/3.0); 
           double rhoa_43 = pow( rhoa, 4.0/3.0); 
           double rhob_43 = pow( rhob, 4.0/3.0); 
           double rhoa_2  = pow( rhoa, 2.0); 
           double rhob_2  = pow( rhob, 2.0); 
           double Xa = sqrt(sigma_aap[p]) / rhoa_43;
           double Xb = sqrt(sigma_bbp[p]) / rhob_43;

           exc += -pow(2.0,1.0/3.0) * Cx * ( pow( rho_ap[p], 4.0/3.0) + pow( rho_bp[p], 4.0/3.0) ) * grid_w_->pointer()[p]; 

           vrho_ap[p] = -1.2407009817988000333 * rhoa_13;
           vrho_bp[p] = -1.2407009817988000333 * rhob_13;
           v2rho_a2p[p] = -0.41356699393293334443 * rhoa_23;
           v2rho_b2p[p] = -0.41356699393293334443 * rhob_23;
           v2rho_abp[p] = 0.0;
       }
    }else {

       throw PsiException("The differentiation order should be 0, 1 or 2!",__FILE__,__LINE__);
    }
    ex_b88 = exc;   
    return ex_b88;
}

double DFTSolver::EX_RPBE(){

    const double MU = 0.2195149727645171;
    const double KAPPA = 0.804;

    vrho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    vsigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_a2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2sigma_aa2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_a_sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * sigma_aap = sigma_aa_->pointer();
    double * sigma_abp = sigma_ab_->pointer();
    double * sigma_bbp = sigma_bb_->pointer();

    double * vrho_ap = vrho_a_->pointer();
    double * vsigma_aap = vsigma_aa_->pointer();

    double * v2rho_a2p = v2rho_a2_->pointer();
    double * v2sigma_aa2p = v2sigma_aa2_->pointer();
    double * v2rho_a_sigma_aap = v2rho_a_sigma_aa_->pointer();

    const double Cx = 0.73855876638202240586;
    double exc = 0.0;

    if (deriv_ == 0) {

       for (int p = 0; p < phi_points_; p++) {

           // double sig = std::max( 0.0, sigma_aap[p]);
           double sig = sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p];
           double rho = rho_ap[p] + rho_bp[p];
           double rho13 = pow( rho, 1.0/3.0);
           double rho23 = pow( rho, 2.0/3.0);
           double rho43 = pow( rho, 4.0/3.0);
           double rho2  = pow( rho, 2.0);

           exc += -Cx * rho43 * (1.0 + KAPPA - KAPPA/( 1.0 + 0.71318265876004893556e-2 * sig / rho23 / rho2 ) ) * grid_w_->pointer()[p];
       }
    }else if (deriv_ == 1) {

       for (int p = 0; p < phi_points_; p++) {

           // double sig = std::max( 0.0, sigma_aap[p]);
           double sig = sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p];
           double rho = rho_ap[p] + rho_bp[p];
           double rho13 = pow( rho, 1.0/3.0);
           double rho23 = pow( rho, 2.0/3.0);
           double rho43 = pow( rho, 4.0/3.0);
           double rho2  = pow( rho, 2.0);

           double f =  ( 1.0 + 0.71318265876004893556e-2 * sig / rho23 / rho2 );
           double g = pow(f,-2.0);
           double F_PBE = 1.0 + KAPPA - KAPPA/f ;

           exc += -Cx * rho43 * F_PBE * grid_w_->pointer()[p];

           vrho_ap[p] = -0.9847450218426965 * rho13 * F_PBE + 0.01129303341188623/ rho13 / rho2 * g * sig;
           vsigma_aap[p] = -0.01693955011782934 / rho43 * g;
       }
    }else if (deriv_ == 2) {

       for (int p = 0; p < phi_points_; p++) {

           double sig = sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p];
           double rho = rho_ap[p] + rho_bp[p];
           double rho13 = pow( rho, 1.0/3.0);
           double rho23 = pow( rho, 2.0/3.0);
           double rho43 = pow( rho, 4.0/3.0);
           double rho2  = pow( rho, 2.0);
           double rho4  = pow( rho, 4.0);

           double f =  ( 1.0 + 0.71318265876004893556e-2 * sig / rho23 / rho2 );
           double g = pow(f,-2.0);
           double h = 1.0 / rho13 / rho2 * g;
           double i = g / f;
           double sig2 = pow( sigma_aap[p], 2.0);

           double F_PBE = 1.0 + KAPPA - KAPPA/f ;

           exc += -Cx * rho43 * F_PBE * grid_w_->pointer()[p];

           v2rho_a2p[p] = ( -0.6564966812284644 / rho23 * F_PBE ) - ( 0.02258606682377246 / rho13 / rho2 / rho * g * sig ) +
           0.8590928633765426e-3 / rho4 / rho2 * i * sig2;
           v2rho_a_sigma_aap[p] = 0.02258606682377246 * h - 0.644319647532407e-3 / rho4 / rho * i * sig;
           v2sigma_aa2p[p] = 0.9664794712986105e-3 / rho4 * i;
       }
    }else {

       throw PsiException("The differentiation order should be 0, 1 or 2!",__FILE__,__LINE__);
    }
    return exc;
}

double DFTSolver::EX_UPBE(){

    const double MU = 0.2195149727645171;
    const double KAPPA = 0.804;

    vrho_a_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    vrho_b_   = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    vsigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    vsigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    vsigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    v2rho_a2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_b2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2sigma_aa2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2sigma_bb2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_a_sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_b_sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_a_sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_b_sigma_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_a_sigma_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2rho_b_sigma_aa_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2sigma_aa_ab_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2sigma_aa_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2sigma_ab2_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));
    v2sigma_ab_bb_ = (std::shared_ptr<Vector>)(new Vector(phi_points_));

    double * rho_ap = rho_a_->pointer();
    double * rho_bp = rho_b_->pointer();
    double * sigma_aap = sigma_aa_->pointer();
    double * sigma_abp = sigma_ab_->pointer();
    double * sigma_bbp = sigma_bb_->pointer();

    double * vrho_ap = vrho_a_->pointer();
    double * vrho_bp = vrho_b_->pointer();
    double * vsigma_aap = vsigma_aa_->pointer();
    double * vsigma_abp = vsigma_ab_->pointer();
    double * vsigma_bbp = vsigma_bb_->pointer();
    double * v2rho_a2p = v2rho_a2_->pointer();
    double * v2rho_abp = v2rho_ab_->pointer();
    double * v2rho_b2p = v2rho_b2_->pointer();
    double * v2sigma_aa2p = v2sigma_aa2_->pointer();
    double * v2sigma_bb2p = v2sigma_bb2_->pointer();
    double * v2rho_a_sigma_aap = v2rho_a_sigma_aa_->pointer();
    double * v2rho_b_sigma_bbp = v2rho_b_sigma_bb_->pointer();
    double * v2rho_a_sigma_abp = v2rho_a_sigma_ab_->pointer();
    double * v2rho_b_sigma_abp = v2rho_b_sigma_ab_->pointer();
    double * v2rho_a_sigma_bbp = v2rho_a_sigma_bb_->pointer();
    double * v2rho_b_sigma_aap = v2rho_b_sigma_aa_->pointer();
    double * v2sigma_aa_abp = v2sigma_aa_ab_->pointer();
    double * v2sigma_aa_bbp = v2sigma_aa_bb_->pointer();
    double * v2sigma_ab2p = v2sigma_ab2_->pointer();
    double * v2sigma_ab_bbp = v2sigma_ab_bb_->pointer();

    const double Cx = 0.73855876638202240586;
    double exc = 0.0;

    if (deriv_ == 0) {

       for (int p = 0; p < phi_points_; p++) {

           // double sig = std::max( 0.0, sigma_aap[p]);
           double sig = sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p];
           double rhoa = rho_ap[p];
           double rhob = rho_bp[p];
           double rhoa_13 = pow( rhoa, 1.0/3.0);
           double rhob_13 = pow( rhob, 1.0/3.0);
           double rhoa_23 = pow( rhoa, 2.0/3.0);
           double rhob_23 = pow( rhob, 2.0/3.0);
           double rhoa_43 = pow( rhoa, 4.0/3.0);
           double rhob_43 = pow( rhob, 4.0/3.0);
           double rhoa_2  = pow( rhoa, 2.0);
           double rhob_2  = pow( rhob, 2.0);

           // exc += -Cx * rho43 * (1.0 + KAPPA - KAPPA/( 1.0 + 0.71318265876004893556e-2 * sig/rho23/rho2 ) ) * grid_w_->pointer()[p]; 
           exc += ( -pow( 2.0, 1.0/3.0) * Cx * rhoa_43 * (1.0 + KAPPA - KAPPA/( 1.0 + 0.449276922095889e-2 * sigma_aap[p]/rhoa_23/rhoa_2 ) )
                    -pow( 2.0, 1.0/3.0) * Cx * rhob_43 * (1.0 + KAPPA - KAPPA/( 1.0 + 0.449276922095889e-2 * sigma_bbp[p]/rhob_23/rhob_2 ) ) ) * grid_w_->pointer()[p];
       }
   }else if (deriv_ == 1) {

       for (int p = 0; p < phi_points_; p++) {

           double sig = sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p];
           double rhoa = rho_ap[p];
           double rhob = rho_bp[p];
           double rhoa_13 = pow( rhoa, 1.0/3.0);
           double rhob_13 = pow( rhob, 1.0/3.0);
           double rhoa_23 = pow( rhoa, 2.0/3.0);
           double rhob_23 = pow( rhob, 2.0/3.0);
           double rhoa_43 = pow( rhoa, 4.0/3.0);
           double rhob_43 = pow( rhob, 4.0/3.0);
           double rhoa_2  = pow( rhoa, 2.0);
           double rhob_2  = pow( rhob, 2.0);

           double fa =  ( 1.0 + 0.449276922095889e-2 * sigma_aap[p] / rhoa_23 / rhoa_2 );
           double fb =  ( 1.0 + 0.449276922095889e-2 * sigma_bbp[p] / rhob_23 / rhob_2 );
           double ga = pow(fa,-2.0);
           double gb = pow(fb,-2.0);
           double F_PBEa = 1.0 + KAPPA - KAPPA/fa ;
           double F_PBEb = 1.0 + KAPPA - KAPPA/fb ;

           // exc += -Cx * rho43 * (1.0 + KAPPA - KAPPA/( 1.0 + 0.71318265876004893556e-2 * sig/rho23/rho2 ) ) * grid_w_->pointer()[p]; 
           exc += -pow( 2.0, 1.0/3.0) * Cx * ( rhoa_43 * F_PBEa + rhob_43 * F_PBEb ) * grid_w_->pointer()[p];

           vrho_ap[p] = -1.2407009817988000333 * rhoa_13 * F_PBEa + 0.8963286558970112e-2 / rhoa_13 / rhoa_2 * ga * sigma_aap[p];
           vrho_bp[p] = -1.2407009817988000333 * rhob_13 * F_PBEb + 0.8963286558970112e-2 / rhob_13 / rhob_2 * gb * sigma_bbp[p];
           vsigma_aap[p] = -0.3361232459613792e-2 / rhoa_43 * ga;
           vsigma_bbp[p] = -0.3361232459613792e-2 / rhob_43 * gb;
           vsigma_abp[p] = 0.0;
       }
    }else if (deriv_ == 2) {

       for (int p = 0; p < phi_points_; p++) {

           double sig = sigma_aap[p] + 2.0 * sigma_abp[p] + sigma_bbp[p];
           double rhoa = rho_ap[p];
           double rhob = rho_bp[p];
           double rhoa_13 = pow( rhoa, 1.0/3.0);
           double rhob_13 = pow( rhob, 1.0/3.0);
           double rhoa_23 = pow( rhoa, 2.0/3.0);
           double rhob_23 = pow( rhob, 2.0/3.0);
           double rhoa_43 = pow( rhoa, 4.0/3.0);
           double rhob_43 = pow( rhob, 4.0/3.0);
           double rhoa_2  = pow( rhoa, 2.0);
           double rhob_2  = pow( rhob, 2.0);
           double rhoa_4  = pow( rhoa, 4.0);
           double rhob_4  = pow( rhob, 4.0);

           double fa =  ( 1.0 + 0.449276922095889e-2 * sigma_aap[p] / rhoa_23 / rhoa_2 );
           double fb =  ( 1.0 + 0.449276922095889e-2 * sigma_bbp[p] / rhob_23 / rhob_2 );
           double ga = pow(fa,-2.0);
           double gb = pow(fb,-2.0);
           double ha = 1.0 / rhoa_13 / rhoa_2 * ga;
           double hb = 1.0 / rhob_13 / rhob_2 * gb;
           double ia = ga / fa;
           double ib = gb / fb;
           double siga_2 = pow( sigma_aap[p], 2.0);
           double sigb_2 = pow( sigma_bbp[p], 2.0);
           double F_PBEa = 1.0 + KAPPA - KAPPA/fa ;
           double F_PBEb = 1.0 + KAPPA - KAPPA/fb ;

           // exc += -Cx * rho43 * (1.0 + KAPPA - KAPPA/( 1.0 + 0.71318265876004893556e-2 * sig/rho23/rho2 ) ) * grid_w_->pointer()[p]; 
           exc += -pow( 2.0, 1.0/3.0) * Cx * ( rhoa_43 * F_PBEa + rhob_43 * F_PBEb ) * grid_w_->pointer()[p];

           vrho_ap[p] = -1.2407009817988000333 * rhoa_13 * F_PBEa + 0.8963286558970112e-2 * ha * sigma_aap[p];
           vrho_bp[p] = -1.2407009817988000333 * rhob_13 * F_PBEb + 0.8963286558970112e-2 * hb * sigma_bbp[p];
           vsigma_aap[p] = -0.3361232459613792e-2 / rhoa_43 * ga;
           vsigma_bbp[p] = -0.3361232459613792e-2 / rhob_43 * gb;
           vsigma_abp[p] = 0.0;

           v2rho_a2p[p] = ( -0.4135669939329333 / rhoa_23 * F_PBEa ) - ( 0.8963286558970112e-2 / rhoa_13 / rhoa_2 / rhoa * ga * sigma_aap[p] ) +
           0.2147732158441357e-3 / rhoa_4 / rhoa_2 * ia * siga_2;
           v2rho_b2p[p] = ( -0.4135669939329333 / rhob_23 * F_PBEb ) - ( 0.8963286558970112e-2 / rhob_13 / rhob_2 / rhob * gb * sigma_bbp[p] ) +
           0.2147732158441357e-3 / rhob_4 / rhob_2 * ib * sigb_2;
           v2rho_abp[p] = 0.0;
           v2rho_a_sigma_aap[p] = 0.4481643279485056e-2 * ha - 0.8053995594155087e-4 / rhoa_4 / rhoa * ia * sigma_aap[p];
           v2rho_b_sigma_bbp[p] = 0.4481643279485056e-2 * hb - 0.8053995594155087e-4 / rhob_4 / rhob * ib * sigma_bbp[p];
           v2sigma_aa2p[p] = 0.3020248347808158e-4 / rhoa_4 * ia;
           v2sigma_bb2p[p] = 0.3020248347808158e-4 / rhob_4 * ib;
           v2rho_a_sigma_abp[p] = 0.0;
           v2rho_b_sigma_abp[p] = 0.0;
           v2rho_a_sigma_bbp[p] = 0.0;
           v2rho_b_sigma_aap[p] = 0.0;
           v2sigma_aa_abp[p] = 0.0;
           v2sigma_aa_bbp[p] = 0.0;
           v2sigma_ab2p[p] = 0.0;
           v2sigma_ab_bbp[p] = 0.0;
       }
    }else {

       throw PsiException("The differentiation order should be 0, 1 or 2!",__FILE__,__LINE__);
    }
    return exc;
}


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++ Correlation Functionals +++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

}} // End namespaces
