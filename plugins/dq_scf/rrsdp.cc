#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "sdp.h"

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace dq_scf {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "DQ_SCF"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C" 
PsiReturnType dq_scf(Options& options)
{

    boost::shared_ptr<SDPSolver> sdp ( new SDPSolver(Process::environment.wavefunction(),options) );
    //Process::environment.set_wavefunction(sdp);
    sdp->compute_energy();

    return Success;
}

}} // End namespaces

