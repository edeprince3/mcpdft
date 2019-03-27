<p align="center">
<img src="logo.png" style='height: 30%; width: 50%; object-fit: contain'/> 
<br>
<a href="https://opensource.org/licenses/GPL-2.0"><img src="https://img.shields.io/github/license/edeprince3/v2rdm_casscf.svg" /></a>
<!-- <img src="logo.png" alt="RDM-inoles" align="middle"> -->
</p>

# Multi-configurational Pair-density Functional Theory (MCPDFT)
A plugin to Psi4

## OVERVIEW

MCPDFT provides an efficient way of recovering both static and dynamical correlation effects with reasonable cost. It becomes most useful in strongly-correlated systems where the number of active electrons and orbitals in the active space is large. In principle, the multi-configurational reference 1-electron (1-RDMs) and 2-electron (2-RDMs) reduced-density matrices can be provided by any method which is able to calculate them. Both translated and fully-translated versions of Slater and Vosko-Wilk-Nusair random-phase approximation expression III (SVWN3), Perdew-Burke-Ernzerhof (PBE), revised PBE (revPBE), Becke88 exchange and one-parameter correlation functional (BOP) and Becke and Lee-Yang-Parr (BLYP) on-top pair-density exchange-correlation functionals are available at the moment. In addition, the global-, double- and range-separated hybrid multi-configurational on-top pair density functionals such as wPBE and LRC-wPBE have also been implemented. However, this part of the project also is under the ongoing developement.


## INSTALLATION

To run the Psi4 plugin MCPDFT:

* Download Psi4 (1.2rc2 or later) from github.com: https://github.com/psi4/psi4, and follow the installation instructions given here: http://psicode.org/psi4manual/master/build_planning.html . Make sure to keep the name of the plugin directory mcpdft .

*  Configure with CMake to generate a Makefile. Run `psi4 --plugin-compile` to get a CMake command. Modify it as needed with `-D` for compiler, libraries, and options.

* Note that, if you configured Psi4 with a fortran compiler, you shouldn't have to specify these things here. If the configure shows no errors, compile the plugin:

  > make

* We recommend the intel compilers (icc, icpc and ifort) for the above procedures.

## INPUT OPTIONS

* **MCPDFT_METHOD** (string):

    The type of the multi-configurational on-top pair-density functional theory adopted for the calculation.
    The legitimate values include MCPDFT, 1H_MCPDFT, 1DH_MCPDFT, RS_MCPDFT, RS1H_MCPDFT, RS1DH_MCPDFT, and LS1DH_MCPDFT. The default
    value is MCPDFT.

* **MCPDFT_REFERENCE** (string):

    The type of reference 1- and 2-RDMs provided for MCPDFT calculations. The allowed values are V2RDM and CI. The
    default is v2RDM.

* **MCPDFT_TYPE** (string):

    The algorithm used for the MCPDFT computation. The valid options are DF and PK. Default is DF.    

* **MCPDFT_FUNCTIONAL** (string):

    The on-top pair-density exchange-correlation functionals used for MCPDFT. The valid options are SVWN, PBE, REVPBE
    BOP, BLYP, WPBE, LRC_WPBE. The default value is SVWN.

* **MCPDFT_TRANSLATION** (string):

    The density transformation scheme for MCPDFT. The valid options are REGULAR and FULL. The default is REGULAR.

* **MCPDFT_OMEGA** (double):

    The range-separation parameter. The default value is 0.0.

* **MCPDFT_LAMBDA** (double):

    The global-hybrid coupling parameter for hybrid functionals. The default value is 0.0.

* **WRITE_QTAIM_WFN** (bool):

    Writes the QTAIM .wfn file (AIMPAC and its successors). The default is false.

## REFERENCES

[1] M. Mostafanejad and A. E. DePrince III, J. Chem. Theory Comput. 15, 290-302 (2019). "Combining Pair-Density Functional Theory and Variational Two-Electron Reduced-Density Matrix Methods"
