/*                                                                      
 *@BEGIN LICENSE
 *                                                                                                                             
 * PSI4: an ab initio quantum chemistry software package
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                                                                                                                     
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *   
 *@END LICENSE
 */
// AED:
// this code performs a boundary point semidefinite programing  on reduced 
//density functional method (PRL 106,083001, 2011)

#ifndef BPSDP_H
#define BPSDP_H
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>
#include<libmints/vector.h>
namespace boost {
  template<class T> class shared_ptr;
}
namespace psi{ namespace bpsdp{

    class BPSDPsolver: public Wavefunction{
    public: 
      BPSDPsolver(boost::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
      ~BPSDPsolver();
      void set_init_param();
      double compute_energy();
      virtual bool same_a_b_orbs() const { return true; }
      virtual bool same_a_b_dens() const { return true; } //true b/c closed shell systems
      long int maxiter;
    protected:
      long int memory;
      long int ndocc,nso,nmo,nvirt,ndoccact,nfrzc,nfrzv,nconstraints;
      double escf, tr2dm, enuc, *eps;
      // Matrices 
      SharedMatrix sMat; //Overlap matrix
      SharedMatrix tMat;  // kinetic energy matrix
      SharedMatrix vMat;  // electron- nuclear potential
      SharedMatrix K2;    // two electron repulsion integral
      SharedMatrix D2;    // two electron reduced density
      SharedMatrix K1;     // one electron integral in MO
      SharedMatrix H;    // one electron integral in AO
      // methods
      void Guess();
      void K();
      void BuildD2();
      double K2D2();
      double TraceD2();
      //Character table
      int *symmetry;
      int Totalsymetry(int i,int j,int k,int l);
      int sym_pdt(int i, int j);
      //container to hold geminal indices
      std::vector <std::vector <std::pair<int,int> > > gems;

    };
  }
}

#endif	
      
