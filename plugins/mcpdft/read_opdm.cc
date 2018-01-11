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

#include <fstream>
#include <iostream>
#include <iomanip>

#include "mcpdft_solver.h"

using namespace psi;

namespace psi{namespace mcpdft{

void MCPDFTSolver::ReadOPDM(double * D1a, double * D1b){

    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));

    std::shared_ptr<PSIO> psio (new PSIO());

    if ( !psio->exists(PSIF_V2RDM_D1A) ) throw PsiException("No D1a on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D1B) ) throw PsiException("No D1b on disk",__FILE__,__LINE__);

    // D1a

    psio->open(PSIF_V2RDM_D1A,PSIO_OPEN_OLD);

    long int na;
    psio->read_entry(PSIF_V2RDM_D1A,"length",(char*)&na,sizeof(long int));

    opdm * opdm_a = (opdm *)malloc(na * sizeof(opdm));
    psio->read_entry(PSIF_V2RDM_D1A,"D1a",(char*)opdm_a,na * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1A,1);

    for (int n = 0; n < na; n++) {
        int i = opdm_a[n].i;
        int j = opdm_a[n].j;
        D1a[i*nmo_+j] = opdm_a[n].val;
    }

    // D1b

    psio->open(PSIF_V2RDM_D1B,PSIO_OPEN_OLD);

    long int nb;
    psio->read_entry(PSIF_V2RDM_D1B,"length",(char*)&nb,sizeof(long int));

    opdm * opdm_b = (opdm *)malloc(nb * sizeof(opdm));
    psio->read_entry(PSIF_V2RDM_D1B,"D1b",(char*)opdm_b,nb * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1B,1);

   for (int n = 0; n < nb; n++) {
        int i = opdm_b[n].i;
        int j = opdm_b[n].j;
        D1b[i*nmo_+j] = opdm_b[n].val;
    }

    BuildRhoFast(opdm_a,opdm_b,na,nb);

    free(opdm_a);
    free(opdm_b);

}

void MCPDFTSolver::ReadCIOPDM(double* D, const char* fileName) {

    std::ifstream dataIn;

    dataIn.open(fileName);

    if (!dataIn)
       std::cout << "Error opening file.\n";
    else {
         for (int i = 0; i < nmo_; i++)
             for (int j = 0; j < nmo_; j++) {

                 dataIn >> D[i*nmo_+j];
                 if (D[i*nmo_+j] < 1.e-20)
                     D[i*nmo_+j] = 0.0;
             }
        dataIn.close();
    }
}

}} //end namespaces


