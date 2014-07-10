#include <sbpsdp.h>
#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace bpsdp {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "BPSDP"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C" 
PsiReturnType bpsdp(Options& options)
{
  boost::shared_ptr<BPSDPsolver > bpsdp (new BPSDPsolver(Process::environment.wavefunction(),options));
  double energy = bpsdp->compute_energy();
    Process::environment.globals["CURRENT ENERGY"] = energy;
  //  int print = options.get_int("PRINT");

    /* Your code goes here */

    return Success;
}

}} // End namespaces

