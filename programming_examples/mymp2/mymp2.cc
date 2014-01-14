#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace mymp2 {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "MYMP2"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C" 
PsiReturnType mymp2(Options& options)
{
    int print = options.get_int("PRINT");

    /* Your code goes here */

    // what do we need?  

    // orbital energies!

    boost::shared_ptr<Wavefunction> ref = Process::environment.wavefunction();
    boost::shared_ptr<Vector> eps = ref->epsilon_a();
    int nmo = ref->nmo();
    int ndocc = ref->doccpi()[0];
    int nvirt = nmo - ndocc;
    printf("\n");
    printf("  no. occupied orbitals: %5i\n",ndocc);
    printf("  no. virtual orbitals:  %5i\n",nvirt);
    printf("\n");
    printf("  orbital energies:\n");
    for (int i = 0; i < nmo; i++) {
        printf("%20.12lf\n",eps->pointer()[i]);
    }
    printf("\n");
    
    // two-electron integrals!
    boost::shared_ptr<MintsHelper> mints (new MintsHelper());
    boost::shared_ptr<Matrix> V = mints->mo_eri(ref->Ca(),ref->Ca());
    double ** Vp = V->pointer();
    double * epsp = eps->pointer();

    int o = ndocc; 
    int v = nvirt;
    double emp2_aa = 0.0;
    double emp2_bb = 0.0;
    double emp2_ab = 0.0;

    // alpha
    for (int a = 0; a < o; a++) {
        // alpha
        for (int b = 0; b < o; b++) {
            // alpha
            for (int r = 0; r < v; r++) {
                // alpha
                for (int s = 0; s < v; s++) {
                    double denom = epsp[a] + epsp[b] - epsp[r+o] - epsp[s+o];
                    double numer = Vp[a*nmo+(r+o)][b*nmo+(s+o)] - Vp[a*nmo+(s+o)][b*nmo+(r+o)];
                    emp2_aa += 0.25 * numer * numer / denom;
                    // closed-shell (beta = alpha)
                    emp2_bb += 0.25 * numer * numer / denom;
                }
            }
        }
    }
    // alpha
    for (int a = 0; a < o; a++) {
        // beta
        for (int b = 0; b < o; b++) {
            // alpha
            for (int r = 0; r < v; r++) {
                // beta
                for (int s = 0; s < v; s++) {
                    double denom = epsp[a] + epsp[b] - epsp[r+o] - epsp[s+o];
                    double numer = Vp[a*nmo+(r+o)][b*nmo+(s+o)];
                    emp2_ab += numer * numer / denom;
                }
            }
        }
    }
    printf("  alpha-alpha contribution: %20.12lf\n",emp2_aa);
    printf("  beta-beta contribution:   %20.12lf\n",emp2_bb);
    printf("  alpha-beta contribution:  %20.12lf\n",emp2_ab);
    printf("\n");
    printf("  mp2 correlation energy:         %20.12lf\n",emp2_aa+emp2_bb+emp2_ab);
    double escf = Process::environment.globals["SCF TOTAL ENERGY"];
    printf("  scf energy:                     %20.12lf\n",escf);
    printf("  total mp2 energy:               %20.12lf\n",emp2_aa+emp2_bb+emp2_ab+escf);
    printf("\n");

    return Success;
}

}} // End namespaces

