#
# @BEGIN LICENSE
#
# myplugin by Psi4 Developer, a plugin to:
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
from psi4.driver.procedures import proc_util

def run_mean_field_cqed(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    mean_field_cqed can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('mean_field_cqed')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    # Your plugin's psi4 run sequence goes here
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    mean_field_cqed_wfn = psi4.core.plugin('mean_field_cqed.so',ref_wfn)

    optstash.restore()

    return mean_field_cqed_wfn


# Integration with driver routines
psi4.driver.procedures['energy']['mean_field_cqed'] = run_mean_field_cqed
psi4.driver.procedures['energy']['cqed'] = run_mean_field_cqed

def exampleFN():
    # Your Python code goes here
    pass
