/*
 *@BEGIN LICENSE
 *
 * mcpdft, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
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
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 * 
 *@END LICENSE
 *
 */

#include "psi4/psi4-dec.h"
#include <psi4/psifiles.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>

#include "mcpdft_solver.h"

using namespace psi;

namespace psi{namespace mcpdft{

struct tpdm {
    int i;
    int j;
    int k;
    int l;
    double val;
};

void MCPDFTSolver::ReadTPDM(double * D2aa, double * D2bb, double * D2ab, double * D1a, double * D1b){

    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    std::shared_ptr<PSIO> psio (new PSIO());

    if ( !psio->exists(PSIF_V2RDM_D2AB) ) throw PsiException("No D2ab on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D2AA) ) throw PsiException("No D2aa on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D2BB) ) throw PsiException("No D2bb on disk",__FILE__,__LINE__);

    //Ca_->print();

    psio_address addr_aa = PSIO_ZERO;
    psio_address addr_bb = PSIO_ZERO;
    psio_address addr_ab = PSIO_ZERO;

    // ab
    psio->open(PSIF_V2RDM_D2AB,PSIO_OPEN_OLD);

    long int nab;
    psio->read_entry(PSIF_V2RDM_D2AB,"length",(char*)&nab,sizeof(long int));

    for (int n = 0; n < nab; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2AB,"D2ab",(char*)&d2,sizeof(tpdm),addr_ab,&addr_ab);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2ab[id] = d2.val;
    }
    psio->close(PSIF_V2RDM_D2AB,1);

    // aa
    psio->open(PSIF_V2RDM_D2AA,PSIO_OPEN_OLD);

    long int naa;
    psio->read_entry(PSIF_V2RDM_D2AA,"length",(char*)&naa,sizeof(long int));

    for (int n = 0; n < naa; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2AA,"D2aa",(char*)&d2,sizeof(tpdm),addr_aa,&addr_aa);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2aa[id] = d2.val;
    }
    psio->close(PSIF_V2RDM_D2AA,1);

    // bb
    psio->open(PSIF_V2RDM_D2BB,PSIO_OPEN_OLD);

    long int nbb;
    psio->read_entry(PSIF_V2RDM_D2BB,"length",(char*)&nbb,sizeof(long int));

    for (int n = 0; n < nbb; n++) {
        tpdm d2;
        psio->read(PSIF_V2RDM_D2BB,"D2bb",(char*)&d2,sizeof(tpdm),addr_bb,&addr_bb);
        int i = d2.i;
        int j = d2.j;
        int k = d2.k;
        int l = d2.l;
        long int id = i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l;
        D2bb[id] = d2.val;
    }
    psio->close(PSIF_V2RDM_D2BB,1);

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            for (int k = 0; k < nmo_; k++) {
                for (int l = 0; l < nmo_; l++) {
                    int ik = i*nmo_+k;
                    int jl = j*nmo_+l;
                    //double dum = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]-D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];

                    double aa_1 = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                    double aa_2 = D2aa[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l];
                    double aa_3 = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+k];
                    double aa_4 = D2aa[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+l*nmo_+k];

                    double ab_1 = D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l];
                    double ab_2 = D2ab[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+l];
                    double ab_3 = D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+k];
                    double ab_4 = D2ab[j*nmo_*nmo_*nmo_+i*nmo_*nmo_+l*nmo_+k];
                    double val = 0.5 * (ab_1-ab_2-ab_3+ab_4);
                    double dum = D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]-val;

                    // print here check against ci
                    //printf("%5i %5i %20.12lf\n",ik,jl,D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]);
                    //printf("%5i %5i %5i %5i %20.12lf %20.12lf %20.12lf %20.12lf\n",i,j,k,l,ab_1,ab_2,ab_3,ab_4);
                    //if ( fabs(dum) > 1e-6 ) {
                    //    //printf("%5i %5i %20.12lf %20.12lf\n",ik,jl,D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l],D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l]);
                    //    printf("%5i %5i %5i %5i %20.12lf %20.12lf\n",i,j,k,l,D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+k*nmo_+l],val);
                    //}
                }
            }
        }
    }

    // check traces:
    double traa = 0.0;
    double trbb = 0.0;
    double trab = 0.0;
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            traa += D2aa[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
            trbb += D2bb[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
            trab += D2ab[i*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+j];
        }
    }
    //printf("  tr(d2aa) = %20.12lf\n",traa);
    //printf("  tr(d2bb) = %20.12lf\n",trbb);
    //printf("  tr(d2ab) = %20.12lf\n",trab);

    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));

    double tra = 0.0;
    double trb = 0.0;

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {

            double duma = 0.0;
            double dumb = 0.0;
            for (int k = 0; k < nmo_; k++) {
                duma += D2ab[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];
                duma += D2aa[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];

                dumb += D2ab[k*nmo_*nmo_*nmo_+i*nmo_*nmo_+k*nmo_+j];
                dumb += D2bb[i*nmo_*nmo_*nmo_+k*nmo_*nmo_+j*nmo_+k];
            }
            D1a[i*nmo_+j] = 1.0/(nalpha_+nbeta_-1.0) * duma;
            D1b[i*nmo_+j] = 1.0/(nalpha_+nbeta_-1.0) * dumb;

        }

        tra += D1a[i*nmo_+i];
        trb += D1b[i*nmo_+i];

    }

    //printf("  tr(da) = %20.12lf\n",tra);
    //printf("  tr(db) = %20.12lf\n",trb);

}

}} //end namespaces


