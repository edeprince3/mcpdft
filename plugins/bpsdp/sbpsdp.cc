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

#include "sbpsdp.h"
#include<psifiles.h>
#include<libpsio/psio.hpp>
#include<libmints/mints.h>
#include<libmints/vector.h>
#include<libmints/matrix.h>
#include<../bin/fnocc/blas.h>
#include<time.h>
#include<libmints/basisset_parser.h>

using namespace psi;
using namespace bpsdp;

namespace psi{ namespace bpsdp{
    BPSDPsolver::BPSDPsolver(boost::shared_ptr<Wavefunction> reference_wavefunction,Options & options):
      Wavefunction(options,_default_psio_lib_){
      reference_wavefunction_ = reference_wavefunction;
      set_init_param();
    }
    BPSDPsolver::~BPSDPsolver()
    {
    }

    double BPSDPsolver::compute_energy() {
      K();
      return 0.0;
    }
    void  BPSDPsolver::set_init_param(){
      // This function is a copy of common_init() from SDPSolver
      escf    = reference_wavefunction_->reference_energy();
      nalpha_ = reference_wavefunction_->nalpha();
      nbeta_ = reference_wavefunction_->nbeta();
      nalphapi_ = reference_wavefunction_->nalphapi();
      nbetapi_ = reference_wavefunction_->nbetapi();
      doccpi_ = reference_wavefunction_->doccpi();
      soccpi_ = reference_wavefunction_->soccpi();
      frzcpi_ = reference_wavefunction_->frzcpi();
      nmopi_ = reference_wavefunction_->nmopi();
      nirrep_ = reference_wavefunction_->nirrep();
      nso_ = reference_wavefunction_->nso();
      nsopi_ = reference_wavefunction_->nsopi();
      
      S_ = SharedMatrix(reference_wavefunction_->S());
      Da_ = SharedMatrix(reference_wavefunction_->Da());
      Ca_ = SharedMatrix(reference_wavefunction_->Ca());
      Fa_ = SharedMatrix(reference_wavefunction_->Fa());
      Db_ = SharedMatrix(reference_wavefunction_->Db());
      Cb_ = SharedMatrix(reference_wavefunction_->Cb());
      Fb_ = SharedMatrix(reference_wavefunction_->Fb());
      
      boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
      boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
      boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser,molecule,"BASIS");
      
      //Creat a n integral factory later use to create integral objects 
      boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(aoBasis,aoBasis,aoBasis,aoBasis));
      boost::shared_ptr<SOBasisSet> soBasis(new SOBasisSet(aoBasis,integral));
      const Dimension dimension=soBasis->dimension();

      boost::shared_ptr<MatrixFactory> factory(new MatrixFactory);
      factory->init_with(dimension,dimension);
      sMat = factory->create_shared_matrix("Overlap");
      tMat = factory->create_shared_matrix("Kinetic");
      vMat = factory->create_shared_matrix("potential");
      H = factory->create_shared_matrix("one electron integral in AO");
      K1 = factory->create_shared_matrix("one electron integral in MO");
      K2 = factory->create_shared_matrix("Two electron integral in MO");

      boost::shared_ptr<OneBodySOInt> sOBI(integral->so_overlap());
      boost::shared_ptr<OneBodySOInt> tOBI(integral->so_kinetic());
      boost::shared_ptr<OneBodySOInt> vOBI(integral->so_potential());
     
      //compute the overlap integral and store in sMat
      sOBI->compute(sMat); 
      tOBI->compute(tMat); 
      vOBI->compute(vMat); 

      H->zero();
      K1->zero();
      K2->zero();
      H->add(vMat);
      H->add(tMat);

      epsilon_a_= boost::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
      epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
      epsilon_b_= boost::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
      epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());
      
      nso = nmo = ndocc = nvirt = nfrzc = nfrzv = 0;
      for (int h=0; h<nirrep_; h++){
        nfrzc   += frzcpi_[h];
        nfrzv   += frzvpi_[h];
        nso    += nsopi_[h];
        nmo    += nmopi_[h]-frzcpi_[h]-frzvpi_[h];
	ndocc  += doccpi_[h];
      }
      ndoccact = ndocc - nfrzc;
      nvirt    = nmo - ndoccact;
      if (nfrzc > 0) {
        throw PsiException("bpsdp does not yet work with frozen core",__FILE__,__LINE__);
      }
      if (nmo != nso) {
        throw PsiException("bpsdp does not yet work when nmo!=nso",__FILE__,__LINE__);
    }
      // memory is from process::environment        
      memory = Process::environment.get_memory();
      // set the wavefunction name
      name_ = "BPSDP";
      // qt                                       
      eps = (double*)malloc(2*nmo*sizeof(double));
      int count=0;
      for (int h=0; h<nirrep_; h++){
        for (int norb = frzcpi_[h]; norb<doccpi_[h]; norb++){
	  eps[count++] = epsilon_a_->get(h,norb);
        }
      }
      for (int h=0; h<nirrep_; h++){
        for (int norb = doccpi_[h]; norb<nmopi_[h]-frzvpi_[h]; norb++){
	  eps[count++] = epsilon_a_->get(h,norb);
	}
      }
      for (int h=0; h<nirrep_; h++){
	for (int norb = doccpi_[h]; norb<nmopi_[h]-frzvpi_[h]; norb++){
	  eps[count++] = epsilon_a_->get(h,norb);
        }
      }
      // Initialize more parameter from here as we move on
      free(eps);
    }
    //build one and two electron integral
     void BPSDPsolver::K(){
      //Compute the dot product 1RDM 12 int and 2RDM 2e int
       sMat->print();       
       tMat->print();       
       vMat->print();       
       H->print();       

       
       /*  for (long int i=0 ; i <nmo ; i++) {
	for (long int j = 0; j <nmo ; j++) {
	  double sum = 0.0;
	  for (long int ii=0 ; ii <nmo ; ii++) {
	    for (long int jj = 0; jj <nmo ; jj++) {
	      sum += Ca_->get(ii,i)*Ca_->get(jj,j)*H->get(ii,jj);
	    }
	  }
	  K1->set(i,j,sum);
	}
      }
      K1->print();*/
      //compute and hold the two electron integral in a buffer
      boost::shared_ptr<MintsHelper> mints (new MintsHelper());
      K2 = mints->mo_eri(Ca_,Ca_);
      K2->print();
      enuc = Process::environment.molecule()->nuclear_repulsion_energy();
      //      printf("nuclear rep energy %20.9lf\n", enuc);
    }
  }
}




