import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util


def run_dq_scf(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    dq_scf can be called via :py:func:`~driver.energy`.

    >>> energy('dq_scf')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # override symmetry:
    molecule = psi4.get_active_molecule()
    molecule.update_geometry()
    molecule.reset_point_group('c1')
    molecule.fix_orientation(1)
    molecule.update_geometry()

    # Your plugin's psi4 run sequence goes here
    #scf_helper(name, **kwargs)
    returnvalue = psi4.plugin('dq_scf.so')

    return returnvalue


# Integration with driver routines
procedures['energy']['dq_scf'] = run_dq_scf


def exampleFN():
    # Your Python code goes here
    pass
