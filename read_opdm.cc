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
#include <psi4/libmints/matrix.h>
#include "psi4/libmints/vector.h"
#include <psi4/libpsi4util/PsiOutStream.h>

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libqt/qt.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "mcpdft_solver.h"

using namespace psi;
using namespace std;

namespace psi{namespace mcpdft{

void MCPDFTSolver::ReadOPDM() {

    std::shared_ptr<PSIO> psio (new PSIO());

    // TODO: should be added back when reading-in the density PSIOH files 
    // psio->set_pid("18332");

    if ( !psio->exists(PSIF_V2RDM_D1A) ) throw PsiException("No D1a on disk",__FILE__,__LINE__);
    if ( !psio->exists(PSIF_V2RDM_D1B) ) throw PsiException("No D1b on disk",__FILE__,__LINE__);

    // D1a

    psio->open(PSIF_V2RDM_D1A,PSIO_OPEN_OLD);

    long int na;
    psio->read_entry(PSIF_V2RDM_D1A,"length",(char*)&na,sizeof(long int));

    opdm_a_ = (opdm *)malloc(na * sizeof(opdm));
    psio->read_entry(PSIF_V2RDM_D1A,"D1a",(char*)opdm_a_,na * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1A,1);

    for (int n = 0; n < na; n++) {

        int i = opdm_a_[n].i;
        int j = opdm_a_[n].j;

        int hi = symmetry_[i];
        int hj = symmetry_[j];

        if ( hi != hj ) {
            throw PsiException("error: something is wrong with the symmetry of the alpha OPDM",__FILE__,__LINE__);
        }

        int ii = i - pitzer_offset_[hi];
        int jj = j - pitzer_offset_[hi];

        Da_->pointer(hi)[ii][jj] = opdm_a_[n].val;

    }

    // D1b

    psio->open(PSIF_V2RDM_D1B,PSIO_OPEN_OLD);

    long int nb;
    psio->read_entry(PSIF_V2RDM_D1B,"length",(char*)&nb,sizeof(long int));

    opdm_b_ = (opdm *)malloc(nb * sizeof(opdm));
    psio->read_entry(PSIF_V2RDM_D1B,"D1b",(char*)opdm_b_,nb * sizeof(opdm));
    psio->close(PSIF_V2RDM_D1B,1);

   for (int n = 0; n < nb; n++) {

        int i = opdm_b_[n].i;
        int j = opdm_b_[n].j;

        int hi = symmetry_[i];
        int hj = symmetry_[j];

        if ( hi != hj ) {
            throw PsiException("error: something is wrong with the symmetry of the beta OPDM",__FILE__,__LINE__);
        }

        int ii = i - pitzer_offset_[hi];
        int jj = j - pitzer_offset_[hi];

        Db_->pointer(hi)[ii][jj] = opdm_b_[n].val;

    }

    if ( !is_low_memory_ ) {
        outfile->Printf("\n");
        outfile->Printf("    ==> Build Rho's ...\n");
        BuildRhoFast(na,nb);
        outfile->Printf("    ... Done. <==\n\n");

        free(opdm_a_);
        free(opdm_b_);
    }

}

void MCPDFTSolver::ReadCIOPDM(std::shared_ptr<Matrix> D, const char* fileName) {

    std::ifstream dataIn;

    dataIn.open(fileName);

    if (!dataIn) throw PsiException("No D1 on disk",__FILE__,__LINE__);
    else {
        double ** dp = D->pointer();
        for (int i = 0; i < nmo_; i++) {
            for (int j = 0; j < nmo_; j++) {
                dataIn >> dp[i][j];
                if (dp[i][j] < 1.0e-20)
                    dp[i][j] = 0.0;
            }
        }
        dataIn.close();
    }
}


void MCPDFTSolver::WriteQTAIM(std::shared_ptr<Matrix> Ca,
                              std::shared_ptr<Matrix> Cb,
                              std::shared_ptr<Vector> Ea,
                              std::shared_ptr<Vector> Eb,
                              std::shared_ptr<Vector> Na,
                              std::shared_ptr<Vector> Nb,
                              const char* fileName) {

    FILE * outfile;
    outfile = fopen(fileName,"w");


    /* ==============
     * Title card
     * ============== */
    // std::string filename = get_writer_file_prefix(reference_wavefunction_->molecule()->name());
    fprintf(outfile," Title Card Required\n");

    /* ==============================
     * Nuclear coordinates & charges
     * ============================== */
    fprintf(outfile,"GAUSSIAN               %d MOL ORBITALS      %d PRIMITIVES         %d NUCLEI\n");
    
    // std::shared_ptr<Matrix> geom = std::shared_ptr<Matrix>(new Matrix(molecule_->full_geometry()));
    // geom->print();
    int nAtoms = molecule_->nunique();
    double geomX, geomY, geomZ;
    double charge;
    std::string labelAtom;
    for ( int i = 0; i < nAtoms; i++ ) {
        labelAtom = molecule_->fsymbol(i);
        geomX = molecule_->xyz(i,0);
        geomY = molecule_->xyz(i,1);
        geomZ = molecule_->xyz(i,2);
        charge = molecule_->fcharge(i);

        fprintf(outfile,"  %s    %d    (CENTER  %d)   ",labelAtom.c_str(),i+1,i+1);
        fprintf(outfile,"% 10.8lf  % 10.8lf  % 10.8lf  CHARGE =   %-5.2f\n",geomX,geomY,geomZ,charge);
    }

    /* ==============================================================
     * Note: Coefficients of basis fxns in each MO belong to a column 
     * of the C matrix. Each column of C matrix therefore should be 
     * arranged in a row of five coefficients in the wfn file as 
     * shown below. 
     * ============================================================== */
    int nRows = 5;
    
    double ** Cap = Ca->pointer();
    double ** Cbp = Cb->pointer();
    double *  Eap = Ea->pointer();
    double *  Ebp = Eb->pointer();
    double *  Nap = Na->pointer();
    double *  Nbp = Nb->pointer();

    for (int i = 1; i < nmo_+1; i++) {
        fprintf(outfile,"MO    %d                    OCC NO =     %15.8lf  ",i,Nap[i-1]);
        fprintf(outfile,"ORB. ENERGY =    %15.8lf\n",Eap[i-1]);
        fflush(stdout);

        for (int j = 1; j < nmo_+1; j++) {
            fprintf(outfile,"%15.8le ",Cap[j-1][i-1]);

            if(j % nRows == 0 )
              fprintf(outfile,"\n");
        }
    }
    fclose(outfile);
}

}} //end namespaces


