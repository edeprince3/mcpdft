#
# @BEGIN LICENSE
#
# mydft by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util
from psi4.driver.procrouting import proc

def run_mydft(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    mydft can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('mydft')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here

    func = psi4.core.get_option('MYDFT','FUNCTIONAL')

    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    #print('Attention! This SCF may be density-fitted.')

    # TODO: this is where empty reference wave function should be constructed
    # but something changed recently in psi4 and this piece of code is broken

    #scf_molecule = kwargs.get('molecule', psi4.core.get_active_molecule())
    #base_wfn = psi4.core.Wavefunction.build(scf_molecule, psi4.core.get_global_option('BASIS'))
    #scf_wfn = proc.scf_wavefunction_factory(psi4.core.get_option('SCF', 'REFERENCE'), base_wfn, func)


    # TODO: scf_wfn should be set above, not here.
    scf_wfn = kwargs.get('ref_wfn', None)
    #if ref_wfn is None:
    #    ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # Ensure IWL files have been written when not using DF/CD
    proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), scf_wfn)

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    mydft_wfn = psi4.core.plugin('mydft.so', scf_wfn)

    return mydft_wfn


# Integration with driver routines
psi4.driver.procedures['energy']['mydft'] = run_mydft


def exampleFN():
    # Your Python code goes here
    pass
