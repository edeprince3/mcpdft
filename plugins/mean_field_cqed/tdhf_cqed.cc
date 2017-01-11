#include "psi4/libplugin/plugin.h"
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
//#include "psi4/libmints/mints.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"

#include "frozen_natural_orbitals.h"
#include "tdhf.h"

INIT_PLUGIN

using namespace std;
using namespace psi;

namespace psi{ namespace tdhf_cqed {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "MEAN_FIELD_CQED"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        /*- total time in a.u. -*/
        options.add_double("TOTAL_TIME", 100.0);
        /*- time step in a.u. -*/
        options.add_double("TIME_STEP", 0.2);
        /*- pulse shape -*/
        options.add_str("LASER_SHAPE", "SIN_SQUARED", "SIN_SQUARED TRAPEZOID PI_PULSE CONTINUOUS GAUSSIAN");
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
        /*- get the spectrum for the molecule only -*/
        options.add_bool("MOLECULE_ONLY",false);
        /*- get the spectrum for the plasmon only -*/
        options.add_bool("PLASMON_ONLY",false);
        /*- polarization (default x+y+z). -*/
        options.add("POLARIZATION",new ArrayType());
        /*- plasmonic-molecule distance (bohr)-*/
        options.add("PLASMON_COORDINATES",new ArrayType());
        /*- plasmonic states -*/
        options.add_int("N_PLASMON_STATES", 1);
        /*- plasmon excitation energy (a.u.) -*/
        options.add("PLASMON_E", new ArrayType());
        /*- coupling energy (a.u.) -*/
        options.add_double("COUPLING_E", 10.8e-3/27.21138);
        /*- plasmon transition dipole moment (a.u.) -*/
        options.add("PLASMON_TDM",/* 2990.0/2.54175*/new ArrayType());
        /*- plasmon damping rate (a.u.) -*/
        options.add_double("PLASMON_DR", 150e-3/27.21138);
        /*- electron damping rate (a.u.) -*/
        options.add_double("ELECTRON_DR", 1e-6);
        /*- dielectric constant of the medium -*/
        options.add_double("EPSILON_M", 0.079577472);
    }

    return true;
}

extern "C" 
SharedWavefunction mean_field_cqed(SharedWavefunction ref_wfn, Options& options)
{
    int print = options.get_int("PRINT");

    tstart();

    // 3-index integrals (generated/read by fno class)
    std::shared_ptr<DFFrozenNO> fno(new DFFrozenNO(ref_wfn,options));
    fno->ThreeIndexIntegrals();
    fno.reset();
    //std::shared_ptr<TDHF> tdhf ( new TDHF(Process::environment.wavefunction(),options) );

    std::shared_ptr<TDHF> tdhf (new TDHF(ref_wfn,options));
    tdhf->compute_energy();

    tstop();

    return fno;
}

}} // End namespaces

