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

#include "dft.h"

namespace psi{ namespace mydft {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "MYDFT"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        /*- Do use DIIS extrapolation? -*/
        options.add_bool("DIIS", true);
        /*- DFT functional -*/
        options.add_str("FUNCTIONAL", "B3LYP");
        /*- DFT functional -*/
        options.add_bool("IP_FITTING", false);
        /*- DFT functional -*/
        options.add_double("DFT_OMEGA", 0.0);
    }

    return true;
}

extern "C"
SharedWavefunction mydft(SharedWavefunction ref_wfn, Options& options)
{

    outfile->Printf("\n\n");
    outfile->Printf( "        *******************************************************\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    mydft                                            *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    An unrestricted Kohn-Sham DFT plugin to Psi4     *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *    Eugene DePrince                                  *\n");
    outfile->Printf( "        *                                                     *\n");
    outfile->Printf( "        *******************************************************\n");

    std::shared_ptr<DFTSolver> dft (new DFTSolver(ref_wfn,options));
    dft->compute_energy();

    return dft;

}

}}
