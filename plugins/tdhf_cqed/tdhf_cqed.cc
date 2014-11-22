#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <../bin/fnocc/frozen_natural_orbitals.h>

#include "tdhf.h"

INIT_PLUGIN

using namespace boost;
using namespace psi;
using namespace fnocc;

namespace psi{ namespace tdhf_cqed {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "TDHF_CQED"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        /*- total time in a.u. -*/
        options.add_double("TOTAL_TIME", 100.0);
        /*- time step in a.u. -*/
        options.add_double("TIME_STEP", 0.2);
        /*- pulse shape -*/
        options.add_str("LASER_SHAPE", "SIN_SQUARED", "SIN_SQUARED TRAPEZOID PI_PULSE CONTINUOUS");
        /*- transition dipole moment for pi pulse -*/
        options.add_double("LASER_TDPM", -0.415638122584);
        /*- amplitude of pulse in a.u.-*/
        options.add_double("LASER_AMP", 0.05);
        /*- frequency of pulse in a.u. (default is the 1 fs pulse) -*/
        options.add_double("LASER_FREQ", 0.1519829846);
        /*- width of pulse in a.u. (1 fs default) -*/
        options.add_double("LASER_TIME", 41.3414);
        /*- flag for linear-response absorption instead of
            real-time interaction with some external field -*/
        options.add_bool("LINEAR_RESPONSE",false);
        /*- polarization (default x+y+z). -*/
        options.add("POLARIZATION",new ArrayType());

    }

    return true;
}

extern "C" 
PsiReturnType tdhf_cqed(Options& options)
{
    int print = options.get_int("PRINT");

    // 3-index integrals (generated/read by fno class)
    boost::shared_ptr<DFFrozenNO> fno(new DFFrozenNO(Process::environment.wavefunction(),options));
    fno->ThreeIndexIntegrals();
    fno.reset();
    //boost::shared_ptr<TDHF> tdhf ( new TDHF(Process::environment.wavefunction(),options) );
    MyTDHF = (boost::shared_ptr<TDHF>) (new TDHF(Process::environment.wavefunction(),options) );
    MyTDHF->compute_energy();

    /* Your code goes here */

    return Success;
}

}} // End namespaces

